#!/bin/bash
# Fine-tuning script for small dataset

# Set which GPU to use (0 for first GPU, 1 for second GPU, etc.)
export CUDA_VISIBLE_DEVICES=5

# Set port for distributed training (use different port if running multiple jobs)
export MASTER_PORT=9995

python finetune_prostt5_lora_script.py \
    --trainaafasta train_aa_9k_1k.faa \
    --trainthreedifasta train_3di_9k_1k.faa \
    --validaafasta valid_aa_9k_1k.faa \
    --validthreedifasta valid_3di_9k_1k.faa \
    -o test_prostt5_lora_9k_1k -b 1 \
    --finetune_name prostt5_finetuned_model_9k_1k -f 2>&1 | tee finetune_prostt5_9k_1k_log.txt
