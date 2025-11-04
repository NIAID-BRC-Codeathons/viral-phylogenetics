# Installation Guide

This guide will help you set up the polyprotein cleavage prediction environment.

## Prerequisites

- Linux, macOS, or Windows with WSL2
- Miniconda or Anaconda installed
- 8GB+ RAM recommended
- GPU with 4GB+ VRAM (optional but recommended)
- 10GB+ free disk space

## Quick Installation

### 1. Clone or Download

If you have the project files, make sure you have:
```
polyprotein_cleavage/
├── esmretrain.py
├── data_prep.py
├── environment.yml
├── config.json
└── *.md documentation files
```

### 2. Create Conda Environment

For full installation with GPU support:
```bash
conda env create -f environment.yml
conda activate polyprotein_cleavage
```

For minimal CPU-only installation:
```bash
conda env create -f environment-minimal.yml
conda activate polyprotein_minimal
```

### 3. Verify Installation

Test the setup:
```bash
python -c "import torch; print(f'PyTorch: {torch.__version__}')"
python -c "import transformers; print(f'Transformers: {transformers.__version__}')"
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
```

### 4. Test ESM Installation

```bash
python -c "import esm; print('ESM package installed successfully')"
```

## Detailed Installation

### Option 1: Full Environment (Recommended)

This installs all dependencies including GPU support, visualization tools, and development utilities.

```bash
# Create environment
conda env create -f environment.yml

# Activate environment
conda activate polyprotein_cleavage

# Verify CUDA support (if you have a GPU)
python -c "import torch; print(torch.cuda.is_available())"
```

**Includes:**
- PyTorch with CUDA support
- ESM protein language models
- Full scientific Python stack
- Bioinformatics tools
- Visualization libraries
- Development tools

### Option 2: Minimal Environment

For testing or CPU-only systems:

```bash
# Create minimal environment
conda env create -f environment-minimal.yml

# Activate environment
conda activate polyprotein_minimal
```

**Includes:**
- PyTorch CPU-only
- Essential ML libraries
- ESM models
- Basic utilities

### Option 3: Manual Installation

If you prefer to install manually:

```bash
# Create new environment
conda create -n polyprotein_cleavage python=3.9

# Activate environment
conda activate polyprotein_cleavage

# Install PyTorch (check pytorch.org for latest versions)
# For GPU:
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia

# For CPU only:
# conda install pytorch torchvision torchaudio cpuonly -c pytorch

# Install other dependencies
conda install -c conda-forge transformers numpy pandas scikit-learn biopython tqdm
pip install esm fair-esm
```

## Troubleshooting

### Common Issues

**1. CUDA/GPU Issues**

```bash
# Check CUDA installation
nvidia-smi

# Check PyTorch CUDA support
python -c "import torch; print(torch.cuda.is_available())"

# If CUDA not available, install CPU version
conda install pytorch torchvision torchaudio cpuonly -c pytorch
```

**2. ESM Installation Issues**

```bash
# Try alternative ESM installation
pip uninstall esm fair-esm
pip install fair-esm

# Or install from source
git clone https://github.com/facebookresearch/esm.git
cd esm
pip install -e .
```

**3. Memory Issues**

For systems with limited memory:
```bash
# Edit config.json to reduce memory usage
{
  "training": {
    "batch_size": 8,  # Reduce from 32
  },
  "data": {
    "max_sequence_length": 2000,  # Reduce from 5000
  }
}
```

**4. Import Errors**

```bash
# Reinstall problematic packages
conda remove transformers
conda install -c huggingface transformers

# Or use pip
pip install --upgrade transformers
```

### Dependency Conflicts

If you encounter dependency conflicts:

```bash
# Remove environment and recreate
conda env remove -n polyprotein_cleavage
conda env create -f environment.yml

# Or update existing environment
conda env update -f environment.yml --prune
```

## Development Setup

For development and debugging:

```bash
# Install additional development tools
conda activate polyprotein_cleavage
conda install jupyter ipython pytest flake8 black

# Install pre-commit hooks (optional)
pip install pre-commit
pre-commit install
```

## Alternative Installations

### Docker (Advanced)

Create a Dockerfile:
```dockerfile
FROM nvidia/cuda:11.8-runtime-ubuntu20.04

RUN apt-get update && apt-get install -y \
    python3 python3-pip wget

RUN pip3 install torch torchvision torchaudio \
    transformers esm biopython numpy pandas scikit-learn

COPY . /app
WORKDIR /app
```

Build and run:
```bash
docker build -t polyprotein-cleavage .
docker run --gpus all -it polyprotein-cleavage
```

### Pip-only Installation

If you can't use conda:
```bash
# Create virtual environment
python -m venv polyprotein_env
source polyprotein_env/bin/activate  # Linux/Mac
# polyprotein_env\Scripts\activate  # Windows

# Install requirements
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install transformers esm biopython numpy pandas scikit-learn tqdm matplotlib
```

## Testing Installation

### Quick Test

```bash
python simple_demo.py
```

### Comprehensive Test

Create a test script:
```python
#!/usr/bin/env python3
import sys

def test_installation():
    try:
        import torch
        print(f"✓ PyTorch {torch.__version__}")
        
        import transformers
        print(f"✓ Transformers {transformers.__version__}")
        
        import esm
        print("✓ ESM package")
        
        import numpy as np
        print(f"✓ NumPy {np.__version__}")
        
        import pandas as pd
        print(f"✓ Pandas {pd.__version__}")
        
        from sklearn import __version__ as sklearn_version
        print(f"✓ Scikit-learn {sklearn_version}")
        
        if torch.cuda.is_available():
            print(f"✓ CUDA available: {torch.cuda.device_count()} devices")
        else:
            print("⚠ CUDA not available (CPU only)")
        
        print("\n🎉 Installation successful!")
        return True
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        return False

if __name__ == "__main__":
    success = test_installation()
    sys.exit(0 if success else 1)
```

Save as `test_installation.py` and run:
```bash
python test_installation.py
```

## Performance Optimization

### For Training

```bash
# Set environment variables for better performance
export OMP_NUM_THREADS=1
export CUDA_LAUNCH_BLOCKING=0

# For mixed precision (if supported)
export TORCH_CUDNN_V8_API_ENABLED=1
```

### Memory Optimization

```python
# In Python scripts
import torch
torch.backends.cudnn.benchmark = True  # For consistent input sizes
torch.backends.cudnn.deterministic = False  # For speed over reproducibility
```

## Uninstallation

To remove the environment:
```bash
conda env remove -n polyprotein_cleavage
```

To remove downloaded models:
```bash
rm -rf ~/.cache/torch/hub/checkpoints/
rm -rf ~/.cache/huggingface/
```

## Next Steps

After successful installation:

1. Read `README_cleavage.md` for usage instructions
2. Follow `QUICKSTART.md` for your first prediction
3. Check `README_refseq.md` for data download instructions
4. Run `simple_demo.py` to test the framework

## Getting Help

If you encounter issues:

1. Check this troubleshooting section
2. Verify your system meets the prerequisites  
3. Try the minimal installation first
4. Check the GitHub issues page
5. Create a new issue with:
   - Your operating system
   - Python version
   - Error messages
   - Installation method used