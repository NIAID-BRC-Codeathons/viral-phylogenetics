#!/bin/bash
# Run inference on A.fasta and B.fasta using fine-tuned ProstT5 model

set -e

# Default model (can be overridden)
MODEL_PATH="${1:-test_prostt5_lora_5k_1k/prostt5_finetuned_model_5k_1k.pth}"
MODEL_DIR="${MODEL_DIR:-model}"
DEVICE="${DEVICE:-cuda}"

echo "=========================================="
echo "ProstT5 Inference on A.fasta and B.fasta"
echo "=========================================="
echo ""
echo "Model: $MODEL_PATH"
echo "Model directory: $MODEL_DIR"
echo "Device: $DEVICE"
echo ""

# Check if model exists
if [ ! -f "$MODEL_PATH" ]; then
    echo "ERROR: Model not found at $MODEL_PATH"
    echo ""
    echo "Available models:"
    echo "  - test_prostt5_lora_5k_1k/prostt5_finetuned_model_5k_1k.pth"
    echo "  - test_prostt5_lora_9k_1k/prostt5_finetuned_model_9k_1k.pth"
    echo ""
    echo "Usage: ./run_inference.sh [model_path]"
    echo "  Example: ./run_inference.sh test_prostt5_lora_9k_1k/prostt5_finetuned_model_9k_1k.pth"
    exit 1
fi

# Check if input files exist
if [ ! -f "A.fasta" ]; then
    echo "ERROR: A.fasta not found"
    exit 1
fi

if [ ! -f "B.fasta" ]; then
    echo "ERROR: B.fasta not found"
    exit 1
fi

# Run inference on A.fasta
echo "Processing A.fasta..."
python inference_finetuned_prostt5.py \
    --input A.fasta \
    --output A_3di_predicted.fasta \
    --model "$MODEL_PATH" \
    --model_dir "$MODEL_DIR" \
    --device "$DEVICE"

echo ""
echo "✓ A.fasta processed -> A_3di_predicted.fasta"
echo ""

# Run inference on B.fasta
echo "Processing B.fasta..."
python inference_finetuned_prostt5.py \
    --input B.fasta \
    --output B_3di_predicted.fasta \
    --model "$MODEL_PATH" \
    --model_dir "$MODEL_DIR" \
    --device "$DEVICE"

echo ""
echo "✓ B.fasta processed -> B_3di_predicted.fasta"
echo ""
echo "=========================================="
echo "Inference complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - A_3di_predicted.fasta"
echo "  - B_3di_predicted.fasta"
echo ""
