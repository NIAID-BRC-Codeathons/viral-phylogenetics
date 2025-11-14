# Viral Phylogenetics - ProstT5 Fine-tuning and Inference

Complete workflow for fine-tuning ProstT5 models and running inference on new sequences.

## Overview

This repository contains scripts for:
1. Splitting full datasets into train/validation sets
2. Creating smaller subsets for debugging
3. Fine-tuning ProstT5 with LoRA
4. Merging LoRA weights with base model
5. Running inference on new sequences

## Important: Where Files Are Saved

**After fine-tuning (`finetune_prostt5*.sh`):**
- **LoRA weights** (`.pth` file): Saved in output directory (e.g., `test_prostt5_lora_5k_1k/prostt5_finetuned_model_5k_1k.pth`)
- **Base model**: Downloaded to `model/` directory (HuggingFace cache, shared across all fine-tuning runs)
- **Training logs**: Saved in output directory

**After merging (`merge_and_save_model.py`):**
- **Complete model**: Saved in `merged_prostt5_model_*/` directory (complete HuggingFace model ready for inference)

**Key point:** The `model/` directory contains the **base ProstT5 model** (downloaded once). The fine-tuning output directories contain **only LoRA weights**. You need to run `merge_and_save_model.py` to combine them into a complete model.

---

## Complete Workflow

### Step 1: Split Full Dataset into Train/Validation

**Script:** `split_fastas.py`

**Purpose:** Split paired AA and 3Di FASTA files into train/validation sets.

**Usage:**
```bash
python split_fastas.py \
    --aa_fasta_path /path_to/bfvd.fasta \
    --di_fasta_path /path_to/bfvd_ss.fasta \
    --train_ratio 0.8 \
    --output_path . \
    --random_seed 42
```

**Output:**
- `train_aa.faa` - Training amino acid sequences
- `valid_aa.faa` - Validation amino acid sequences
- `train_3di.faa` - Training 3Di sequences
- `valid_3di.faa` - Validation 3Di sequences

**Parameters:**
- `--aa_fasta_path`: Input amino acid FASTA file
- `--di_fasta_path`: Input 3Di FASTA file
- `--train_ratio`: Train/validation split ratio (default: 0.8)
- `--output_path`: Output directory (default: current directory)
- `--random_seed`: Random seed for reproducibility

---

### Step 2: Create Smaller Dataset Subset (Optional)

**Script:** `create_small_dataset.py` or `create_dataset_5k_1k.sh`

**Purpose:** Create subsets from the full train/valid splits for debugging.

**Usage (Python script):**
```bash
python create_small_dataset.py \
    --train_aa train_aa.faa \
    --train_3di train_3di.faa \
    --valid_aa valid_aa.faa \
    --valid_3di valid_3di.faa \
    --train_size 5000 \
    --valid_size 1000 \
    --seed 42
```

**Usage (Bash script for 5k/1k):**
```bash
./create_dataset_5k_1k.sh
```

**Output files (with sample sizes encoded in filename):**
- `train_aa_5k_1k.faa` - 5000 training samples
- `train_3di_5k_1k.faa`
- `valid_aa_5k_1k.faa` - 1000 validation samples
- `valid_3di_5k_1k.faa`

**Note:** Filenames encode both train and validation sizes (e.g., `_5k_1k` = 5k train, 1k valid)

---

### Step 3: Fine-tune ProstT5 Model

**Script:** `finetune_prostt5_5k_1k.sh`, `finetune_prostt5_9k_1k.sh`, or `finetune_prostt5.sh`

**Purpose:** Fine-tune ProstT5 with LoRA on your dataset.

**Usage:**
```bash
./finetune_prostt5_5k_1k.sh
```

**What it does:**
- Fine-tunes ProstT5 with LoRA (Low-Rank Adaptation)
- Downloads base ProstT5 model to `model/` directory (HuggingFace cache)
- Saves only trainable parameters (LoRA weights + classifier head) to `.pth` file

**Output files:**
- **LoRA weights:** Saved in the output directory specified by `-o` flag:
  - `finetune_prostt5.sh` → `test_prostt5_lora/prostt5_finetuned_model.pth`
  - `finetune_prostt5_5k_1k.sh` → `test_prostt5_lora_5k_1k/prostt5_finetuned_model_5k_1k.pth`
  - `finetune_prostt5_9k_1k.sh` → `test_prostt5_lora_9k_1k/prostt5_finetuned_model_9k_1k.pth`
- **Base model:** Automatically downloaded to `model/` directory (HuggingFace cache)
- **Training logs:** Also saved in the output directory (`.log`, `loss_plot.png`, etc.)

**Important:** 
- The `.pth` file contains **only LoRA weights**, NOT a complete model
- The `model/` directory contains the **base ProstT5 model** (downloaded from HuggingFace)
- You need **both** to create a complete model (see Step 4)

---

### Step 4: Merge LoRA Weights with Base Model

**Script:** `merge_and_save_model.py`

**Purpose:** Merge fine-tuned LoRA weights with base ProstT5 model and save as complete HuggingFace model.

**Usage:**
```bash
# For 5k_1k dataset:
python merge_and_save_model.py \
    --lora_path test_prostt5_lora_5k_1k/prostt5_finetuned_model_5k_1k.pth \
    --output_dir merged_prostt5_model_5k_1k \
    --num_labels 20 \
    --model_dir model 2>&1 | tee merge_and_save_model_5k_1k.log

# For 9k_1k dataset:
python merge_and_save_model.py \
    --lora_path test_prostt5_lora_9k_1k/prostt5_finetuned_model_9k_1k.pth \
    --output_dir merged_prostt5_model_9k_1k \
    --num_labels 20 \
    --model_dir model 2>&1 | tee merge_and_save_model_9k_1k.log

# For full dataset:
python merge_and_save_model.py \
    --lora_path test_prostt5_lora/prostt5_finetuned_model.pth \
    --output_dir merged_prostt5_model \
    --num_labels 20 \
    --model_dir model 2>&1 | tee merge_and_save_model.log
```

**What it does:**
1. Loads base ProstT5 model from `model/` directory (HuggingFace cache)
2. Loads LoRA weights from the `.pth` file (from Step 3)
3. Merges them together
4. Saves complete model in HuggingFace format

**Output:**
- `merged_prostt5_model_*/` directory containing:
  - `pytorch_model.bin` - Complete model weights
  - `config.json` - Model configuration
  - Tokenizer files (`tokenizer_config.json`, `spiece.model`, etc.)

**Why it's needed:**
- Fine-tuning only saves LoRA weights (incomplete model)
- This step creates the complete model by combining base + LoRA weights
- The merged model can be used for inference or further processing

---

### Step 5: Run Inference on New Sequences

**Script:** `run_inference.sh` or `inference_finetuned_prostt5.py`

**Purpose:** Predict 3Di sequences from amino acid sequences using the fine-tuned model.

**Usage (bash script):**
```bash
# Use default model (5k_1k)
./run_inference.sh

# Use specific model
./run_inference.sh test_prostt5_lora_9k_1k/prostt5_finetuned_model_9k_1k.pth
```

**Usage (Python script directly):**
```bash
python inference_finetuned_prostt5.py \
    --input input_sequences.fasta \
    --output predicted_3di.fasta \
    --model test_prostt5_lora_5k_1k/prostt5_finetuned_model_5k_1k.pth \
    --model_dir model \
    --device cuda
```

**What it does:**
- Loads base ProstT5 model and applies fine-tuned LoRA weights
- Predicts 3Di sequences from input AA sequences
- Outputs predictions in FASTA format

**Note:** You can use the LoRA weights (`.pth` file) directly - no need to merge first for inference.

---

## Quick Reference

### Complete Workflow (Full Dataset)

```bash
# 1. Split full dataset
python split_fastas.py \
    --aa_fasta_path full_aa.faa \
    --di_fasta_path full_3di.faa \
    --train_ratio 0.8

# 2. Fine-tune
./finetune_prostt5.sh

# 3. Merge LoRA weights
python merge_and_save_model.py \
    --lora_path test_prostt5_lora/prostt5_finetuned_model.pth \
    --output_dir merged_prostt5_model \
    --num_labels 20 \
    --model_dir model

# 4. Run inference on new sequences
./run_inference.sh test_prostt5_lora/prostt5_finetuned_model.pth
```

### Quick Workflow (Small Dataset for Debugging)

```bash
# 1. Create small dataset (5k train, 1k valid)
./create_dataset_5k_1k.sh

# 2. Fine-tune on small dataset
./finetune_prostt5_5k_1k.sh

# 3. Merge LoRA weights
python merge_and_save_model.py \
    --lora_path test_prostt5_lora_5k_1k/prostt5_finetuned_model_5k_1k.pth \
    --output_dir merged_prostt5_model_5k_1k \
    --num_labels 20 \
    --model_dir model

# 4. Run inference on new sequences
./run_inference.sh test_prostt5_lora_5k_1k/prostt5_finetuned_model_5k_1k.pth
```

---

## File Naming Conventions

- **Full dataset splits:** `train_aa.faa`, `valid_aa.faa`, `train_3di.faa`, `valid_3di.faa`
- **Small subsets:** `train_aa_5k_1k.faa` (encodes train_size and valid_size)
- **Fine-tuned weights:** `prostt5_finetuned_model_small.pth` (LoRA weights only)
- **Merged model:** `merged_prostt5_model/` (complete HuggingFace model)

---

## Scripts Overview

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `split_fastas.py` | Split full dataset | Full AA/3Di FASTAs | train/valid FASTAs |
| `create_dataset_5k_1k.sh` | Create smaller subsets | train/valid FASTAs | Smaller train/valid FASTAs |
| `finetune_prostt5_5k_1k.sh` | Fine-tune model | train/valid FASTAs | LoRA weights (.pth) |
| `merge_and_save_model.py` | Merge LoRA weights | LoRA weights (.pth) | Complete HuggingFace model |
| `run_inference.sh` | Run inference | LoRA weights (.pth) + AA FASTAs | Predicted 3Di sequences |

---

## Dependencies

- Python 3.10+
- PyTorch
- Transformers
- BioPython
- ProstT5 model (downloaded automatically)

---

## Environment Setup

Activate your conda environment:
```bash
conda activate prostt5  # or your environment name
```

---

## Additional Resources

- `ProstT5/README.md` - Original ProstT5 documentation
