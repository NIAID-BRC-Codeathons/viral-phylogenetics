# Polyprotein Cleavage Site Prediction

This project uses ESM3 protein language model embeddings to predict cleavage sites in viral polyproteins. The approach treats cleavage prediction as a binary classification problem where each amino acid position is classified as either a cleavage site or not.

## Overview

Viral polyproteins are large proteins that are cleaved at specific sites to produce individual functional proteins. Accurate prediction of these cleavage sites is important for:

- Understanding viral protein processing
- Drug target identification
- Vaccine development
- Viral pathogenesis research

## Approach

1. **ESM3 Embeddings**: Extract contextual embeddings for each amino acid position using the ESM3 protein language model
2. **Binary Classification**: Train a neural network to predict cleavage probability for each position
3. **Sequence-level Prediction**: Process entire polyprotein sequences to identify all cleavage sites

## Model Architecture

- **Input**: ESM3 embeddings (1536-dimensional vectors per amino acid)
- **Network**: Feedforward neural network [1536 → 512 → 256 → 128 → 1]
- **Output**: Sigmoid activation producing cleavage probability (0-1)
- **Loss**: Binary cross-entropy with class weight balancing

## Data Format

The training data should be in JSON format with the following structure:

```json
[
  {
    "protein_id": "SARS_CoV_2_polyprotein",
    "sequence": "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQYGRSGETLGVLVPHVGEIPVAYRKVLLRKNGNKGAGGHSYGADLKSFDLGDELGTDPYEDFQENWNTKHSSGVTRELMRELNGGAYVTKSQDGSRQTQVAVGVKFSVKADNDSLFHPGGKHPHVNHLPNFNFTWNLKKIGASNYEFLKQSLKGLKFRSNLLQKGRPQGKFSFLKDTLEEELKDGELIFKSLLLQNHPRNP",
    "cleavage_sites": [15, 48, 89, 156, 234],
    "organism": "SARS-CoV-2"
  }
]
```

### Required Fields

- `protein_id`: Unique identifier for the polyprotein
- `sequence`: Full amino acid sequence (single letter code)
- `cleavage_sites`: List of 0-indexed positions where cleavage occurs
- `organism`: Source organism (optional but recommended)

## Usage

### Basic Training

```python
from esmretrain import PolyproteinCleavagePredictor

# Initialize predictor
predictor = PolyproteinCleavagePredictor()

# Train on your data
predictor.train("your_data.json", epochs=50, batch_size=32)

# Save the trained model
predictor.save_model("cleavage_model.pth")
```

### Making Predictions

```python
# Load trained model
predictor = PolyproteinCleavagePredictor()
predictor.load_model("cleavage_model.pth")

# Predict cleavage sites
sequence = "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEV..."
cleavage_probs = predictor.predict_sequence(sequence)

# Get high-confidence cleavage sites
threshold = 0.5
cleavage_sites = predictor.get_cleavage_sites(sequence, threshold=threshold)
print(f"Predicted cleavage sites: {cleavage_sites}")
```

## Performance Considerations

### Memory Usage

- ESM3 embeddings require significant memory (1536 floats per amino acid)
- For long sequences (>2000 AA), consider processing in chunks
- GPU memory requirements scale with batch size and sequence length

### Training Time

- ESM3 embedding extraction is the bottleneck (requires GPU for reasonable speed)
- Consider pre-computing embeddings for large datasets
- Training the classification network is relatively fast

### Model Size

- ESM3 model: ~2.8GB
- Classification network: <10MB
- Total memory requirement: ~3GB GPU memory for inference

## Evaluation Metrics

The model is evaluated using:

- **Precision**: TP / (TP + FP) - accuracy of predicted cleavage sites
- **Recall**: TP / (TP + FN) - fraction of actual cleavage sites found
- **F1-Score**: Harmonic mean of precision and recall
- **AUROC**: Area under ROC curve for threshold-independent evaluation
- **AUPRC**: Area under precision-recall curve (better for imbalanced data)

## Data Imbalance

Cleavage sites are rare (~1-5% of positions), leading to class imbalance:

- **Class Weights**: Automatically computed based on class frequencies
- **Focal Loss**: Optional alternative to handle hard examples
- **Threshold Tuning**: Optimize classification threshold on validation set

## File Structure

```
polyprotein_cleavage/
├── esmretrain.py           # Main framework
├── data_prep.py            # Data preparation utilities
├── environment.yml         # Full conda environment
├── environment-minimal.yml # Minimal environment
├── README_cleavage.md      # This file
├── config.json            # Model configuration
└── simple_demo.py         # Example usage
```

## Validation Strategy

1. **Cross-validation**: 5-fold CV to assess model stability
2. **Temporal splits**: If data has time stamps, use temporal validation
3. **Organism holdout**: Test on unseen viral species
4. **Family holdout**: Test on unseen viral families

## Troubleshooting

### Common Issues

1. **CUDA out of memory**: Reduce batch size or sequence length
2. **Slow embedding extraction**: Ensure GPU is available and CUDA is properly installed
3. **Poor performance**: Check data quality, class balance, and hyperparameters
4. **Model not learning**: Verify data format, check learning rate, ensure shuffling

### Debug Mode

Enable detailed logging:

```python
import logging
logging.basicConfig(level=logging.DEBUG)

predictor = PolyproteinCleavagePredictor(verbose=True)
```

## Future Improvements

1. **Architecture**: Experiment with transformer-based classifiers
2. **Multi-task Learning**: Predict cleavage type in addition to location
3. **Protein Structure**: Incorporate structural features if available
4. **Active Learning**: Iteratively improve with expert feedback
5. **Transfer Learning**: Pre-train on related protein processing tasks

## References

- Lin et al. "Evolutionary-scale prediction of atomic level protein structure with a language model." (ESM-Fold)
- Rives et al. "Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences." (ESM-1b)
- Hayes et al. "Simulating 500 million years of evolution with a language model." (ESM3)

## Citation

If you use this framework in your research, please cite:

```bibtex
@software{polyprotein_cleavage_2024,
  title={Polyprotein Cleavage Site Prediction with ESM3},
  author={Your Name},
  year={2024},
  url={https://github.com/your-repo/polyprotein-cleavage}
}
```