#!/bin/bash
# Create dataset with 5000 training and 1000 validation samples
# Output files will be named: train_aa_5k_1k.faa, valid_aa_5k_1k.faa, etc.

# ./create_dataset_5k_1k.sh 2>&1 | tee create_dataset_5k_1k.log

set -e  # Exit on error

echo "=========================================="
echo "Creating 5k/1k Dataset"
echo "=========================================="
echo ""

# Default input files (adjust if needed)
TRAIN_AA="${TRAIN_AA:-train_aa.faa}"
TRAIN_3DI="${TRAIN_3DI:-train_3di.faa}"
VALID_AA="${VALID_AA:-valid_aa.faa}"
VALID_3DI="${VALID_3DI:-valid_3di.faa}"

# Sample sizes
TRAIN_SIZE=5000
VALID_SIZE=1000
SEED=42

echo "Configuration:"
echo "  Training samples: $TRAIN_SIZE"
echo "  Validation samples: $VALID_SIZE"
echo "  Random seed: $SEED"
echo ""
echo "Input files:"
echo "  Train AA: $TRAIN_AA"
echo "  Train 3Di: $TRAIN_3DI"
echo "  Valid AA: $VALID_AA"
echo "  Valid 3Di: $VALID_3DI"
echo ""

# Check input files exist
for file in "$TRAIN_AA" "$TRAIN_3DI" "$VALID_AA" "$VALID_3DI"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Input file not found: $file"
        exit 1
    fi
done

# Run the Python script
python create_small_dataset.py \
    --train_aa "$TRAIN_AA" \
    --train_3di "$TRAIN_3DI" \
    --valid_aa "$VALID_AA" \
    --valid_3di "$VALID_3DI" \
    --train_size "$TRAIN_SIZE" \
    --valid_size "$VALID_SIZE" \
    --seed "$SEED"

echo ""
echo "=========================================="
echo "Dataset creation complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  train_aa_5k_1k.faa"
echo "  train_3di_5k_1k.faa"
echo "  valid_aa_5k_1k.faa"
echo "  valid_3di_5k_1k.faa"
echo ""
echo "You can now use these with finetune_prostt5_small.sh"
echo "(update the script to use the new filenames)"
