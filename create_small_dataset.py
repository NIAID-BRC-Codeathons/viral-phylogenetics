#!/usr/bin/env python3
"""
Create smaller subsets of training and validation datasets for debugging.

This script takes the full train/valid FASTA files and creates smaller versions
by randomly sampling a specified number of records from each set.

Usage:
    python create_small_dataset.py \
        --train_aa /path_to/train_aa.faa \
        --train_3di /path_to/train_3di.faa \
        --valid_aa /path_to/valid_aa.faa \
        --valid_3di /path_to/valid_3di.faa \
        --train_size 5000 \
        --valid_size 1000 \
        --seed 42

    Or use a bash script:
    ./create_dataset_5k_1k.sh
"""
import random
import argparse
from pathlib import Path
from Bio import SeqIO


def create_small_dataset(
    input_aa_file,
    input_3di_file,
    output_aa_file,
    output_3di_file,
    num_records,
    seed=None):
    """Create a smaller subset of paired AA and 3Di FASTA files.
    
    Args:
        input_aa_file: Path to input amino acid FASTA file
        input_3di_file: Path to input 3Di FASTA file
        output_aa_file: Path to output amino acid FASTA file
        output_3di_file: Path to output 3Di FASTA file
        num_records: Number of records to sample
        seed: Random seed for reproducibility
    """
    if seed is not None:
        random.seed(seed)
    
    # Read input files
    print(f"Reading {input_aa_file}")
    aa_records = list(SeqIO.parse(input_aa_file, 'fasta'))
    print(f"Reading {input_3di_file}")
    di_records = list(SeqIO.parse(input_3di_file, 'fasta'))
    
    # Validate matching lengths and IDs
    if len(aa_records) != len(di_records):
        raise ValueError(f'FASTAs have different number of records: '
                        f'{len(aa_records)} vs {len(di_records)}')
    
    for aa, di in zip(aa_records, di_records):
        if aa.id != di.id:
            raise ValueError(f'ID mismatch: {aa.id} vs {di.id}')
    
    total_records = len(aa_records)
    num_records = min(num_records, total_records)
    
    if num_records == total_records:
        print(f"Warning: Requested {num_records} records, but only {total_records} available. Using all records.")
    else:
        print(f"Sampling {num_records} records from {total_records} total records")
    
    # Randomly sample indices
    indices = list(range(total_records))
    random.shuffle(indices)
    selected_indices = sorted(indices[:num_records])  # Sort to maintain order
    
    # Select records
    selected_aa = [aa_records[i] for i in selected_indices]
    selected_di = [di_records[i] for i in selected_indices]
    
    # Write output files
    print(f"Writing {output_aa_file}")
    with open(output_aa_file, 'w') as f:
        SeqIO.write(selected_aa, f, 'fasta')
    
    print(f"Writing {output_3di_file}")
    with open(output_3di_file, 'w') as f:
        SeqIO.write(selected_di, f, 'fasta')
    
    print(f"Created small dataset with {num_records} records")
    print(f"  - {output_aa_file}")
    print(f"  - {output_3di_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Create smaller subsets of train/val datasets for debugging'
    )
    parser.add_argument(
        '--train_aa', type=str,
        default='./train_aa.faa',
        help='Input training amino acid FASTA file'
    )
    parser.add_argument(
        '--train_3di', type=str,
        default='./train_3di.faa',
        help='Input training 3Di FASTA file'
    )
    parser.add_argument(
        '--valid_aa', type=str,
        default='./valid_aa.faa',
        help='Input validation amino acid FASTA file'
    )
    parser.add_argument(
        '--valid_3di', type=str,
        default='./valid_3di.faa',
        help='Input validation 3Di FASTA file'
    )
    parser.add_argument(
        '--train_size', type=int, default=2000,
        help='Number of training records to sample (default: 2000)'
    )
    parser.add_argument(
        '--valid_size', type=int, default=400,
        help='Number of validation records to sample (default: 400)'
    )
    parser.add_argument(
        '--seed', type=int, default=42,
        help='Random seed for reproducibility (default: 42)'
    )
    parser.add_argument(
        '--suffix', type=str, default='_small',
        help='Suffix to add to output filenames (default: _small)'
    )
    
    args = parser.parse_args()
    
    # Generate output filenames with sample sizes encoded
    train_aa_path = Path(args.train_aa)
    train_3di_path = Path(args.train_3di)
    valid_aa_path = Path(args.valid_aa)
    valid_3di_path = Path(args.valid_3di)
    
    # Format: train_aa_5k_1k.faa (5k train, 1k valid)
    train_size_str = f"{args.train_size // 1000}k" if args.train_size >= 1000 else str(args.train_size)
    valid_size_str = f"{args.valid_size // 1000}k" if args.valid_size >= 1000 else str(args.valid_size)
    size_suffix = f"_{train_size_str}_{valid_size_str}"
    
    train_aa_small = train_aa_path.parent / f"{train_aa_path.stem}{size_suffix}{train_aa_path.suffix}"
    train_3di_small = train_3di_path.parent / f"{train_3di_path.stem}{size_suffix}{train_3di_path.suffix}"
    valid_aa_small = valid_aa_path.parent / f"{valid_aa_path.stem}{size_suffix}{valid_aa_path.suffix}"
    valid_3di_small = valid_3di_path.parent / f"{valid_3di_path.stem}{size_suffix}{valid_3di_path.suffix}"
    
    # Create small training set
    print("=" * 60)
    print("Creating small TRAINING dataset")
    print("=" * 60)
    create_small_dataset(
        args.train_aa,
        args.train_3di,
        str(train_aa_small),
        str(train_3di_small),
        args.train_size,
        seed=args.seed
    )
    
    # Create small validation set
    print("\n" + "=" * 60)
    print("Creating small VALIDATION dataset")
    print("=" * 60)
    create_small_dataset(
        args.valid_aa,
        args.valid_3di,
        str(valid_aa_small),
        str(valid_3di_small),
        args.valid_size,
        seed=args.seed
    )
    
    print("\n" + "=" * 60)
    print("Small dataset creation complete!")
    print("=" * 60)
    print(f"\nOutput files:")
    print(f"  Training: {train_aa_small}, {train_3di_small}")
    print(f"  Validation: {valid_aa_small}, {valid_3di_small}")
    print(f"\nYou can now use these files with finetune_prostt5_small.sh")


if __name__ == '__main__':
    main()

