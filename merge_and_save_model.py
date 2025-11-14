#!/usr/bin/env python3
"""
Merge base ProstT5 model with fine-tuned LoRA weights and save as HuggingFace format.
This creates a full model that can then be converted to GGUF if needed.

Examples:
  # For 5k_1k dataset:
  python merge_and_save_model.py \
      --lora_path test_prostt5_lora_5k_1k/prostt5_finetuned_model_5k_1k.pth \
      --output_dir merged_prostt5_model_5k_1k \
      --num_labels 20 \
      --model_dir model 2>&1 | tee merge_and_save_model_5k_1k.log

  # For full dataset:
  python merge_and_save_model.py \
      --lora_path test_prostt5_lora/prostt5_finetuned_model.pth \
      --output_dir merged_prostt5_model \
      --num_labels 20 \
      --model_dir model 2>&1 | tee merge_and_save_model.log
"""

import argparse
import torch
import copy
import re
import numpy as np
from pathlib import Path
import sys

# Add the ProstT5 scripts directory to path
sys.path.insert(0, str(Path(__file__).parent / "ProstT5" / "scripts"))

# Import necessary components
from transformers import T5EncoderModel, T5Tokenizer, T5Config, T5PreTrainedModel
from transformers.models.t5.modeling_t5 import T5Stack
from transformers.modeling_outputs import TokenClassifierOutput
import torch.nn as nn
import torch.nn.functional as F

# LoRA components (copied from finetune script since they're defined inside main())
class LoRAConfig:
    def __init__(self):
        self.lora_rank = 4
        self.lora_init_scale = 0.01
        self.lora_modules = ".*SelfAttention|.*EncDecAttention"
        self.lora_layers = "q|k|v|o"
        self.trainable_param_names = ".*layer_norm.*|.*lora_[ab].*"
        self.lora_scaling_rank = 1

class LoRALinear(nn.Module):
    def __init__(self, linear_layer, rank, scaling_rank, init_scale):
        super().__init__()
        self.in_features = linear_layer.in_features
        self.out_features = linear_layer.out_features
        self.rank = rank
        self.scaling_rank = scaling_rank
        self.weight = linear_layer.weight
        self.bias = linear_layer.bias
        if self.rank > 0:
            self.lora_a = nn.Parameter(
                torch.randn(rank, linear_layer.in_features) * init_scale
            )
            if init_scale < 0:
                self.lora_b = nn.Parameter(
                    torch.randn(linear_layer.out_features, rank) * init_scale
                )
            else:
                self.lora_b = nn.Parameter(
                    torch.zeros(linear_layer.out_features, rank)
                )
        if self.scaling_rank:
            self.multi_lora_a = nn.Parameter(
                torch.ones(self.scaling_rank, linear_layer.in_features)
                + torch.randn(self.scaling_rank, linear_layer.in_features) * init_scale
            )
            if init_scale < 0:
                self.multi_lora_b = nn.Parameter(
                    torch.ones(linear_layer.out_features, self.scaling_rank)
                    + torch.randn(linear_layer.out_features, self.scaling_rank) * init_scale
                )
            else:
                self.multi_lora_b = nn.Parameter(
                    torch.ones(linear_layer.out_features, self.scaling_rank)
                )

    def forward(self, input):
        if self.scaling_rank == 1 and self.rank == 0:
            if self.multi_lora_a.requires_grad:
                hidden = F.linear(
                    (input * self.multi_lora_a.flatten()), self.weight, self.bias
                )
            else:
                hidden = F.linear(input, self.weight, self.bias)
            if self.multi_lora_b.requires_grad:
                hidden = hidden * self.multi_lora_b.flatten()
            return hidden
        else:
            weight = self.weight
            if self.scaling_rank:
                weight = (
                    weight * torch.matmul(self.multi_lora_b, self.multi_lora_a) / self.scaling_rank
                )
            if self.rank:
                weight = weight + torch.matmul(self.lora_b, self.lora_a) / self.rank
            input = input.to(torch.half)
            weight = weight.to(torch.half)
            return F.linear(input, weight, self.bias)

def modify_with_lora(transformer, config):
    """Modify transformer by adding LoRA layers."""
    for m_name, module in dict(transformer.named_modules()).items():
        if re.fullmatch(config.lora_modules, m_name):
            for c_name, layer in dict(module.named_children()).items():
                if re.fullmatch(config.lora_layers, c_name):
                    assert isinstance(
                        layer, nn.Linear
                    ), f"LoRA can only be applied to torch.nn.Linear, but {layer} is {type(layer)}."
                    setattr(
                        module,
                        c_name,
                        LoRALinear(
                            layer,
                            config.lora_rank,
                            config.lora_scaling_rank,
                            config.lora_init_scale,
                        ),
                    )
    return transformer


# Define the model classes (copied from finetune script)
class ClassConfig:
    def __init__(self, dropout=0.2, num_labels=3):
        self.dropout_rate = dropout
        self.num_labels = num_labels

class T5EncoderForTokenClassification(T5PreTrainedModel):
    def __init__(self, config: T5Config, class_config):
        super().__init__(config)
        self.num_labels = class_config.num_labels
        self.config = config

        self.shared = nn.Embedding(config.vocab_size, config.d_model)

        encoder_config = copy.deepcopy(config)
        encoder_config.use_cache = False
        encoder_config.is_encoder_decoder = False
        self.encoder = T5Stack(encoder_config, self.shared)

        self.dropout = nn.Dropout(class_config.dropout_rate)
        self.classifier = nn.Linear(config.hidden_size, class_config.num_labels)

        # Initialize weights and apply final processing
        self.post_init()

    def forward(
        self,
        input_ids=None,
        attention_mask=None,
        head_mask=None,
        inputs_embeds=None,
        labels=None,
        output_attentions=None,
        output_hidden_states=None,
        return_dict=None,
    ):
        return_dict = (
            return_dict if return_dict is not None else self.config.use_return_dict
        )

        outputs = self.encoder(
            input_ids=input_ids,
            attention_mask=attention_mask,
            inputs_embeds=inputs_embeds,
            head_mask=head_mask,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict,
        )

        sequence_output = outputs[0]
        sequence_output = self.dropout(sequence_output)
        sequence_output_fl = sequence_output.to(torch.float32)
        logits = self.classifier(sequence_output_fl)

        if not return_dict:
            output = (logits,) + outputs[2:]
            return output

        return TokenClassifierOutput(
            loss=None,
            logits=logits,
            hidden_states=outputs.hidden_states,
            attentions=outputs.attentions,
        )


def load_prostt5_classification_model(num_labels, model_dir, half_precision=False):
    """Load ProstT5 model with classification head (replicates PT5_classification_model logic)."""
    # Load PT5 and tokenizer
    if not half_precision:
        model = T5EncoderModel.from_pretrained(
            f"Rostlab/ProstT5", cache_dir=f"{model_dir}/", use_safetensors=True
        )
        tokenizer = T5Tokenizer.from_pretrained(
            f"Rostlab/ProstT5", cache_dir=f"{model_dir}/"
        )
    else:
        tokenizer = T5Tokenizer.from_pretrained(
            f"Rostlab/ProstT5_fp16", do_lower_case=False, cache_dir=f"{model_dir}/"
        )
        model = T5EncoderModel.from_pretrained(
            f"Rostlab/ProstT5_fp16",
            torch_dtype=torch.float16,
            cache_dir=f"{model_dir}/",
            use_safetensors=True,
        ).to(torch.device("cuda"))

    # Create new Classifier model with PT5 dimensions
    class_config = ClassConfig(num_labels=num_labels)
    class_model = T5EncoderForTokenClassification(model.config, class_config)

    # Set encoder and embedding weights to checkpoint weights
    class_model.shared = model.shared
    class_model.encoder = model.encoder

    # Delete the checkpoint model
    model = class_model
    del class_model

    # Add LoRA layers (but we'll load the trained LoRA weights later)
    config = LoRAConfig()
    model = modify_with_lora(model, config)

    return model, tokenizer


def merge_lora_with_base(
    lora_path,
    output_dir,
    num_labels=20,
    model_dir="model",
    device="cuda"):
    """Merge LoRA weights with base model and save as HuggingFace format.
    
    Args:
        lora_path: Path to LoRA weights .pth file
        output_dir: Directory to save merged model (HuggingFace format)
        num_labels: Number of classification labels
        model_dir: Directory where base ProstT5 model is cached
        device: Device to load model on
    """
    print(f"Loading base ProstT5 model from {model_dir}...")
    
    # Load base model
    model, tokenizer = load_prostt5_classification_model(
        num_labels=num_labels,
        model_dir=model_dir,
        half_precision=False
    )
    
    print(f"Loading LoRA weights from {lora_path}...")
    # Load LoRA weights
    lora_weights = torch.load(lora_path, map_location=device)
    
    # Apply LoRA weights to model
    print("Merging LoRA weights with base model...")
    merged_count = 0
    for param_name, param in model.named_parameters():
        if param_name in lora_weights:
            param.data = lora_weights[param_name].data
            merged_count += 1
    print(f"  Merged {merged_count} parameter groups from LoRA weights")
    
    model = model.to(device)
    model.eval()
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\nSaving merged model to {output_dir}...")
    
    # Save tokenizer first (creates tokenizer_config.json, vocab files, etc.)
    print("  Saving tokenizer...")
    tokenizer.save_pretrained(output_dir)
    
    # Save model using HuggingFace's save_pretrained (proper way)
    # This saves config.json, pytorch_model.bin, and other required files
    print("  Saving model...")
    try:
        model.save_pretrained(output_dir, safe_serialization=False)  # Use .bin format for compatibility
        print("  ✓ Model saved using save_pretrained()")
    except Exception as e:
        print(f"  Warning: save_pretrained() failed: {e}")
        print("  Falling back to manual save...")
        # Fallback: manual save
        torch.save(model.state_dict(), output_dir / "pytorch_model.bin")
        
        # Save config manually with proper structure
        config = model.config if hasattr(model, 'config') else None
        if config:
            import json
            config_dict = config.to_dict()
            # Ensure model_type is set (needed for HuggingFace recognition)
            if 'model_type' not in config_dict:
                config_dict['model_type'] = 't5'
            # Add architecture info if missing
            # Note: Using T5EncoderModel instead of T5EncoderForTokenClassification
            # for compatibility with llama.cpp/convert_hf_to_gguf.py (though it still won't work)
            # The actual model is T5EncoderForTokenClassification, but this helps with some tools
            if 'architectures' not in config_dict:
                config_dict['architectures'] = ['T5EncoderModel']  # Base architecture for compatibility
            # Add custom classifier info
            if hasattr(model, 'num_labels'):
                config_dict['num_labels'] = model.num_labels
            elif 'num_labels' not in config_dict:
                config_dict['num_labels'] = num_labels  # Use function parameter
            
            with open(output_dir / "config.json", "w") as f:
                json.dump(config_dict, f, indent=2)
            print("  ✓ Config saved manually")
        else:
            print("  ✗ Warning: Could not save model config")
    
    # Verify required files exist
    print("\n  Verifying HuggingFace model structure...")
    required_files = ["config.json", "pytorch_model.bin"]
    tokenizer_files = ["tokenizer_config.json"]
    
    all_present = True
    for file in required_files:
        if (output_dir / file).exists():
            print(f"    ✓ {file}")
        else:
            print(f"    ✗ {file} MISSING")
            all_present = False
    
    for file in tokenizer_files:
        if (output_dir / file).exists():
            print(f"    ✓ {file}")
        else:
            print(f"    ⚠ {file} (optional)")
    
    # List all files for debugging
    print(f"\n  All files in {output_dir}:")
    for f in sorted(output_dir.iterdir()):
        size = f.stat().st_size / (1024 * 1024)  # Size in MB
        print(f"    - {f.name} ({size:.2f} MB)")
    
    if not all_present:
        print("\n  ⚠ WARNING: Some required files are missing!")
        print("  The model may not be recognized by HuggingFace/llama.cpp")
    
    print(f"\n✓ Merged model saved to {output_dir}")
    print("\nNote: The official llama.cpp/convert_hf_to_gguf.py does not support T5EncoderForTokenClassification.")
    print("Use convert_gguf.py (custom converter) or inference_finetuned_prostt5.py for inference.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge LoRA weights with base model and save as HuggingFace format"
    )
    parser.add_argument(
        "--lora_path", type=str, required=True,
        help="Path to LoRA weights .pth file"
    )
    parser.add_argument(
        "--output_dir", type=str, required=True,
        help="Output directory for merged model (HuggingFace format)"
    )
    parser.add_argument(
        "--num_labels", type=int, default=20,
        help="Number of classification labels (default: 20)"
    )
    parser.add_argument(
        "--model_dir", type=str, default="model",
        help="Directory where base ProstT5 model is cached"
    )
    parser.add_argument(
        "--device", type=str, default="cuda",
        help="Device to use (default: cuda)"
    )
    
    args = parser.parse_args()
    
    if not torch.cuda.is_available() and args.device == "cuda":
        print("Warning: CUDA not available, using CPU")
        args.device = "cpu"
    
    merge_lora_with_base(
        args.lora_path,
        args.output_dir,
        num_labels=args.num_labels,
        model_dir=args.model_dir,
        device=args.device
    )

    print(f"Finished merging and saving model: {Path(__file__).resolve()}")
