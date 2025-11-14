#!/usr/bin/env python3
"""
Inference script for fine-tuned ProstT5 model.
Loads the fine-tuned LoRA weights and predicts 3Di sequences from amino acid sequences.
"""
import argparse
import torch
from Bio import SeqIO
from Bio.Seq import Seq

# Import model loading functions from merge_and_save_model
# This module has the extracted model classes and functions
from merge_and_save_model import load_prostt5_classification_model


def load_finetuned_model(
    model_path,
    num_labels=20,
    model_dir="model",
    device="cuda"):
    """Load the fine-tuned model with LoRA weights.
    
    Args:
        model_path: Path to the saved .pth file with LoRA weights
        num_labels: Number of classification labels (20 for 3Di)
        model_dir: Directory where base ProstT5 model is cached
        device: Device to load model on
        
    Returns:
        model, tokenizer
    """
    print(f"Loading base ProstT5 model from {model_dir}...")
    
    # Load base model (this will download if not cached)
    model, tokenizer = load_prostt5_classification_model(
        num_labels=num_labels,
        model_dir=model_dir,
        half_precision=False
    )
    
    print(f"Loading fine-tuned weights from {model_path}...")
    # Load the fine-tuned parameters
    non_frozen_params = torch.load(model_path, map_location=device)
    
    # Assign the fine-tuned parameters to the model
    for param_name, param in model.named_parameters():
        if param_name in non_frozen_params:
            param.data = non_frozen_params[param_name].data
            print(f"  Loaded: {param_name}")
    
    model = model.to(device)
    model.eval()
    
    print("Model loaded successfully!")
    return model, tokenizer


def predict_3di(model, tokenizer, aa_sequence, device="cuda"):
    """Predict 3Di sequence from amino acid sequence.
    
    Args:
        model: Fine-tuned ProstT5 model
        tokenizer: T5 tokenizer
        aa_sequence: Amino acid sequence string
        device: Device to run inference on
        
    Returns:
        Predicted 3Di sequence string
    """
    # Replace uncommon AAs with X
    aa_sequence = aa_sequence.replace("O", "X").replace("B", "X").replace("U", "X").replace("Z", "X")
    
    # Add spaces between amino acids for T5 tokenization
    aa_sequence_spaced = " ".join(list(aa_sequence))
    
    # Tokenize (matching training script: max_length=1024, padding=False, truncation=True)
    tokenized = tokenizer(
        aa_sequence_spaced,
        max_length=1024,
        padding=False,
        truncation=True,
        return_tensors="pt"
    )
    
    # Move to device
    input_ids = tokenized["input_ids"].to(device)
    attention_mask = tokenized.get("attention_mask", None)
    if attention_mask is not None:
        attention_mask = attention_mask.to(device)
    
    # Predict
    with torch.no_grad():
        if attention_mask is not None:
            outputs = model(input_ids=input_ids, attention_mask=attention_mask)
        else:
            outputs = model(input_ids=input_ids)
        
        logits = outputs.logits  # Shape: (batch, seq_len, num_labels)
        
        # Get predicted class for each position (argmax)
        predictions = torch.argmax(logits, dim=-1)  # Shape: (batch, seq_len)
        
        # Convert to numpy
        predictions = predictions.cpu().numpy()[0]
        
        # 3Di label mapping (from training script)
        ss_mapping = {
            0: "A", 1: "C", 2: "D", 3: "E", 4: "F", 5: "G", 6: "H", 7: "I",
            8: "K", 9: "L", 10: "M", 11: "N", 12: "P", 13: "Q", 14: "R",
            15: "S", 16: "T", 17: "V", 18: "W", 19: "Y"
        }
        
        # Extract predictions for actual sequence tokens
        # T5 tokenizer adds special tokens at the beginning
        # Each amino acid becomes one token (because we space them)
        # So predictions[1:seq_len+1] should correspond to the sequence
        # But we need to account for the actual tokenized length
        
        # Get actual sequence length (number of tokens, excluding special tokens)
        # The tokenized sequence has: [special_token] + [AA tokens] + [optional special tokens]
        input_ids_np = input_ids.cpu().numpy()[0]
        actual_seq_len = len(aa_sequence)  # Original sequence length
        
        # Skip first token (special token) and take predictions for sequence
        # Match training: labels correspond to positions 1 to seq_len (0-indexed: 1 to seq_len)
        if len(predictions) > 1:
            # Take predictions starting from index 1 (skip special token)
            seq_predictions = predictions[1:actual_seq_len+1]
            # Ensure we don't exceed available predictions
            seq_predictions = seq_predictions[:actual_seq_len]
        else:
            seq_predictions = []
        
        # Convert to 3Di string
        di3_sequence = "".join([ss_mapping.get(int(pred), "X") for pred in seq_predictions])
    
    return di3_sequence


def process_fasta(
    input_fasta,
    output_fasta,
    model_path,
    model_dir="model",
    device="cuda"):
    """Process a FASTA file and generate 3Di predictions.
    
    Args:
        input_fasta: Input FASTA file with amino acid sequences
        output_fasta: Output FASTA file for 3Di sequences
        model_path: Path to fine-tuned model .pth file
        model_dir: Directory where base ProstT5 model is cached
        device: Device to run inference on
    """
    print(f"Loading model...")
    model, tokenizer = load_finetuned_model(model_path, num_labels=20, model_dir=model_dir, device=device)
    
    print(f"Processing sequences from {input_fasta}...")
    
    records = list(SeqIO.parse(input_fasta, "fasta"))
    output_records = []
    
    for i, record in enumerate(records):
        print(f"Processing {i+1}/{len(records)}: {record.id}")
        
        aa_seq = str(record.seq)
        di3_seq = predict_3di(model, tokenizer, aa_seq, device=device)
        
        # Create new record with 3Di sequence (must be Seq object, not string)
        output_record = record
        output_record.seq = Seq(di3_seq)
        output_records.append(output_record)
    
    # Write output
    print(f"Writing predictions to {output_fasta}...")
    with open(output_fasta, "w") as f:
        SeqIO.write(output_records, f, "fasta")
    
    print(f"Done! Generated {len(output_records)} 3Di sequences.")


def main():
    parser = argparse.ArgumentParser(
        description="Predict 3Di sequences from AA sequences using fine-tuned ProstT5"
    )
    parser.add_argument(
        "--input", type=str, required=True,
        help="Input FASTA file with amino acid sequences"
    )
    parser.add_argument(
        "--output", type=str, required=True,
        help="Output FASTA file for 3Di sequences"
    )
    parser.add_argument(
        "--model", type=str, required=True,
        help="Path to fine-tuned model .pth file"
    )
    parser.add_argument(
        "--model_dir", type=str, default="model",
        help="Directory where base ProstT5 model is cached (default: model)"
    )
    parser.add_argument(
        "--device", type=str, default="cuda",
        help="Device to use (default: cuda)"
    )
    
    args = parser.parse_args()
    
    if not torch.cuda.is_available() and args.device == "cuda":
        print("Warning: CUDA not available, using CPU")
        args.device = "cpu"
    
    process_fasta(
        args.input,
        args.output,
        args.model,
        model_dir=args.model_dir,
        device=args.device
    )


if __name__ == "__main__":
    main()
