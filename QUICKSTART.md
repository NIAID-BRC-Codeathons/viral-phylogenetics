# Quick Start Guide

Get up and running with polyprotein cleavage prediction in minutes.

## Prerequisites

✅ Conda environment installed (see [INSTALL.md](INSTALL.md))  
✅ Basic familiarity with Python and command line

## 5-Minute Quick Start

### 1. Setup Environment

```bash
# If not already done
conda env create -f environment.yml
conda activate polyprotein_cleavage
```

### 2. Test Installation

```bash
python simple_demo.py
```

This runs a quick demo to verify everything works.

### 3. Download Sample Data

```bash
# Download 50 viral polyproteins from RefSeq
python data_prep.py refseq \
    --output sample_data.json \
    --max-entries 50 \
    --email your.email@institution.edu
```

⚠️ **Important**: Use your real email address for NCBI requests.

### 4. Train a Model

```bash
python -c "
from esmretrain import PolyproteinCleavagePredictor

# Initialize predictor
predictor = PolyproteinCleavagePredictor(verbose=True)

# Train on sample data
predictor.train('sample_data.json', epochs=10, batch_size=8)

# Save model
predictor.save_model('my_first_model.pth')
print('✅ Training complete!')
"
```

### 5. Make Predictions

```bash
python -c "
from esmretrain import PolyproteinCleavagePredictor

# Load trained model
predictor = PolyproteinCleavagePredictor()
predictor.load_model('my_first_model.pth')

# Example sequence (SARS-CoV-2 polyprotein fragment)
sequence = 'MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVA'

# Predict cleavage sites
sites = predictor.get_cleavage_sites(sequence, threshold=0.5)
print(f'Predicted cleavage sites: {sites}')
"
```

🎉 **Congratulations!** You've successfully trained and used a cleavage prediction model.

## Common Workflows

### Workflow 1: Quick Evaluation

Test the framework with minimal data:

```bash
# 1. Create demo data
python simple_demo.py

# 2. Quick training (few epochs)
python -c "
from esmretrain import PolyproteinCleavagePredictor
predictor = PolyproteinCleavagePredictor()
predictor.train('demo_data.json', epochs=3, batch_size=2)
print('Quick training done!')
"
```

### Workflow 2: Serious Training

For actual research with proper data:

```bash
# 1. Download comprehensive dataset
python data_prep.py families \
    --families Coronaviridae Picornaviridae Flaviviridae \
    --output viral_polyproteins.json \
    --max-per-family 200 \
    --email your.email@institution.edu

# 2. Validate data quality
python data_prep.py validate --input viral_polyproteins.json

# 3. Create train/val/test splits
python data_prep.py split \
    --input viral_polyproteins.json \
    --train 0.7 --val 0.15 --test 0.15

# 4. Train with proper parameters
python -c "
from esmretrain import PolyproteinCleavagePredictor
predictor = PolyproteinCleavagePredictor('config.json', verbose=True)
predictor.train('viral_polyproteins_train.json', epochs=50)
predictor.save_model('production_model.pth')
"

# 5. Evaluate on test set
python -c "
from esmretrain import PolyproteinCleavagePredictor
predictor = PolyproteinCleavagePredictor()
predictor.load_model('production_model.pth')
metrics = predictor.evaluate('viral_polyproteins_test.json')
print('Test Results:', metrics)
"
```

### Workflow 3: Custom Data

If you have your own polyprotein data:

```bash
# 1. Prepare your data in the required JSON format
# See README_cleavage.md for format specification

# 2. Validate format
python data_prep.py validate --input your_data.json

# 3. Train model
python -c "
from esmretrain import PolyproteinCleavagePredictor
predictor = PolyproteinCleavagePredictor()
predictor.train('your_data.json')
predictor.save_model('custom_model.pth')
"
```

## Essential Commands

### Data Commands

```bash
# Download RefSeq data
python data_prep.py refseq --email your@email.com

# Download specific virus families
python data_prep.py families --families Coronaviridae

# Parse FASTA with annotations
python data_prep.py fasta --fasta seqs.fasta --annotations sites.csv

# Split data
python data_prep.py split --input data.json

# Validate data format
python data_prep.py validate --input data.json
```

### Training Commands

```python
from esmretrain import PolyproteinCleavagePredictor

# Basic training
predictor = PolyproteinCleavagePredictor()
predictor.train('data.json')

# Advanced training with custom config
predictor = PolyproteinCleavagePredictor('config.json')
predictor.train('data.json', epochs=100, batch_size=16)

# Training with validation
predictor.train('train.json', validation_file='val.json')
```

### Prediction Commands

```python
# Load model and predict
predictor = PolyproteinCleavagePredictor()
predictor.load_model('model.pth')

# Get cleavage probabilities
probs = predictor.predict_sequence(sequence)

# Get cleavage sites above threshold
sites = predictor.get_cleavage_sites(sequence, threshold=0.5)

# Evaluate on test set
metrics = predictor.evaluate('test.json')
```

## Configuration

### Quick Config Changes

Edit `config.json` for common adjustments:

```json
{
  "training": {
    "batch_size": 16,     // Reduce if memory issues
    "learning_rate": 0.001,
    "epochs": 50
  },
  "data": {
    "max_sequence_length": 2000  // Reduce for memory
  },
  "device": {
    "use_cuda": false     // Set to false for CPU-only
  }
}
```

### Performance Tuning

For different hardware setups:

```bash
# High-end GPU (RTX 3080+)
batch_size: 64, max_sequence_length: 5000

# Mid-range GPU (GTX 1660+)  
batch_size: 32, max_sequence_length: 3000

# Low-end GPU or CPU
batch_size: 8, max_sequence_length: 1000
```

## Troubleshooting

### Quick Fixes

**❌ CUDA out of memory**
```python
# Reduce batch size in config.json
"batch_size": 8
```

**❌ ESM model download fails**
```bash
# Clear cache and retry
rm -rf ~/.cache/huggingface/
python simple_demo.py
```

**❌ No cleavage sites predicted**
```python
# Lower threshold
sites = predictor.get_cleavage_sites(sequence, threshold=0.1)
```

**❌ Training is slow**
```bash
# Check if using GPU
python -c "import torch; print(torch.cuda.is_available())"

# Reduce sequence length
# Edit config.json: "max_sequence_length": 1000
```

### Data Issues

**❌ "No sequences found"**
- Check your email address for RefSeq downloads
- Try broader search terms
- Verify internet connection

**❌ "Validation failed"**
- Check data format in README_cleavage.md
- Ensure cleavage sites are within sequence bounds
- Verify JSON syntax

## Next Steps

### Learn More
- 📖 [README_cleavage.md](README_cleavage.md) - Detailed framework documentation
- 📖 [README_refseq.md](README_refseq.md) - Data download guide
- 📖 [INSTALL.md](INSTALL.md) - Installation troubleshooting

### Advanced Usage
- Experiment with different viral families
- Try different neural network architectures
- Implement custom evaluation metrics
- Add structure-based features

### Research Applications
- Analyze cleavage patterns across viral families
- Compare predictions with experimental data
- Study evolution of cleavage sites
- Design cleavage-resistant sequences

## Example Scripts

### Complete Training Script

Save as `train_model.py`:

```python
#!/usr/bin/env python3
"""Complete training script for polyprotein cleavage prediction"""

from esmretrain import PolyproteinCleavagePredictor
import json

def main():
    # Download data
    print("Downloading data...")
    import subprocess
    subprocess.run([
        "python", "data_prep.py", "families",
        "--families", "Coronaviridae", "Picornaviridae", 
        "--output", "training_data.json",
        "--max-per-family", "100",
        "--email", "your@email.com"  # CHANGE THIS
    ])
    
    # Split data
    print("Splitting data...")
    subprocess.run([
        "python", "data_prep.py", "split",
        "--input", "training_data.json"
    ])
    
    # Train model
    print("Training model...")
    predictor = PolyproteinCleavagePredictor(verbose=True)
    predictor.train("training_data_train.json", 
                   validation_file="training_data_val.json",
                   epochs=30)
    
    # Save model
    predictor.save_model("final_model.pth")
    
    # Evaluate
    print("Evaluating...")
    metrics = predictor.evaluate("training_data_test.json")
    
    # Save results
    with open("results.json", "w") as f:
        json.dump(metrics, f, indent=2)
    
    print("Training complete! Results saved to results.json")

if __name__ == "__main__":
    main()
```

Run with:
```bash
python train_model.py
```

### Batch Prediction Script

Save as `predict_batch.py`:

```python
#!/usr/bin/env python3
"""Batch prediction script"""

from esmretrain import PolyproteinCleavagePredictor
import json
import sys

def predict_file(model_path, input_file, output_file, threshold=0.5):
    # Load model
    predictor = PolyproteinCleavagePredictor()
    predictor.load_model(model_path)
    
    # Load sequences
    with open(input_file) as f:
        data = json.load(f)
    
    # Predict each sequence
    results = []
    for item in data:
        sequence = item["sequence"]
        sites = predictor.get_cleavage_sites(sequence, threshold=threshold)
        
        results.append({
            "protein_id": item.get("protein_id", "unknown"),
            "sequence": sequence,
            "predicted_sites": sites,
            "n_sites": len(sites)
        })
    
    # Save results
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"Predicted {len(results)} sequences")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python predict_batch.py model.pth input.json output.json")
        sys.exit(1)
    
    predict_file(sys.argv[1], sys.argv[2], sys.argv[3])
```

Run with:
```bash
python predict_batch.py final_model.pth sequences.json predictions.json
```

---

🚀 **You're ready to go!** Start with the 5-minute quick start and explore from there.