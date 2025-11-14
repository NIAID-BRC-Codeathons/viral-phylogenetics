#!/bin/bash --login

set -e

# conda create -n prostt5 python=3.10 pip --yes
# conda activate prostt5

pip install --upgrade pip
pip install torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 --index-url https://download.pytorch.org/whl/cu118

pip install "transformers>=4.44" peft accelerate datasets einops sentencepiece tokenizers evaluate
pip install protobuf
pip install biopython scikit-learn pandas numpy tqdm
pip install scikit-bio
pip install ete3 matplotlib

conda install anaconda::evaluate
conda install conda-forge::loguru
