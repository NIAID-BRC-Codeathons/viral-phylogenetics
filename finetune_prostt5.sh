#!/bin/bash
# Fine-tuning script for full dataset

# Set which GPU to use (0 for first GPU, 1 for second GPU, etc.)
export CUDA_VISIBLE_DEVICES=3

# Set port for distributed training (use different port if running multiple jobs)
export MASTER_PORT=9993

python finetune_prostt5_lora_script.py \
    --trainaafasta train_aa.faa \
    --trainthreedifasta train_3di.faa \
    --validaafasta valid_aa.faa \
    --validthreedifasta valid_3di.faa \
    -o test_prostt5_lora -b 1 \
    --finetune_name prostt5_finetuned_model -f 2>&1 | tee finetune_prostt5_log.txt
