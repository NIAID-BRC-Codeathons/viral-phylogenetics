"""
Split FASTA files into train and validation sets.

Usage:
    python split_fastas.py \
        --aa_fasta_path /path_to/bfvd.fasta \
        --di_fasta_path /path_to/bfvd_ss.fasta \
        --train_ratio 0.8 \
        --output_path . \
        --random_seed 42
"""
import random
import argparse
from pathlib import Path
from Bio import SeqIO


def validate_records(aa_records, di_records):
    """Validate that AA and 3Di records match in length and IDs.
    
    Returns:
        True if validation passes
        
    Raises:
        ValueError if records don't match in length or IDs
    """
    if len(aa_records) != len(di_records):
        raise ValueError(f'FASTAs have different number of records: '
                        f'{len(aa_records)} vs {len(di_records)}')
    
    mismatches = []
    for aa, di in zip(aa_records, di_records):
        if aa.id != di.id:
            mismatches.append(f'{aa.id} vs {di.id}')
    
    if mismatches:
        raise ValueError(f'ID mismatches found: {", ".join(mismatches[:5])}'
                        f'{" (and more...)" if len(mismatches) > 5 else ""}')
    
    return True


def split_and_write(
    aa_records,
    di_records,
    train_ratio=0.8,
    output_dir=None,
    seed=None):
    """Split records into train/validation sets and write output files."""
    # Validate train_ratio
    if not 0 < train_ratio < 1:
        raise ValueError(f'train_ratio must be between 0 and 1, got {train_ratio}')

    if seed is not None:
        random.seed(seed)

    # Get indices and shuffle
    num_records = len(aa_records)
    if num_records == 0:
        raise ValueError('No records found in input files')
    if num_records == 1:
        raise ValueError('Cannot split single record into train/val sets.')

    indices = list(range(num_records))
    random.shuffle(indices)

    # Split according to ratio
    split_idx = int(train_ratio * num_records)
    train_indices = indices[:split_idx]
    valid_indices = indices[split_idx:]

    # Set output directory
    if output_dir is None:
        output_dir = Path('.')
    else:
        # Convert to Path if not already (handles both str and Path inputs)
        output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define output files
    output_files = {
        'train_aa': output_dir / 'train_aa.faa',
        'valid_aa': output_dir / 'valid_aa.faa',
        'train_3di': output_dir / 'train_3di.faa',
        'valid_3di': output_dir / 'valid_3di.faa',
    }

    # Write output files
    try:
        with open(output_files['train_aa'], 'w') as f:
            SeqIO.write([aa_records[i] for i in train_indices], f, 'fasta')

        with open(output_files['valid_aa'], 'w') as f:
            SeqIO.write([aa_records[i] for i in valid_indices], f, 'fasta')

        with open(output_files['train_3di'], 'w') as f:
            SeqIO.write([di_records[i] for i in train_indices], f, 'fasta')

        with open(output_files['valid_3di'], 'w') as f:
            SeqIO.write([di_records[i] for i in valid_indices], f, 'fasta')
    except IOError as e:
        raise IOError(f'Error writing output files: {e}')

    return train_indices, valid_indices, output_files


def main():
    parser = argparse.ArgumentParser(
        description='Split paired AA and 3Di FASTA files into train/validation sets'
    )
    parser.add_argument(
        '--aa_fasta_path', type=str,
        default='/nfs/ml_lab/projects/ml_lab/joverbeek/brc_codeathon/bfvd/bfvd_foldseekdb/bfvd.fasta',
        help='Path to amino acid FASTA file'
    )
    parser.add_argument(
        '--di_fasta_path', type=str,
        default='/nfs/ml_lab/projects/ml_lab/joverbeek/brc_codeathon/bfvd/bfvd_foldseekdb/bfvd_ss.fasta',
        help='Path to 3Di FASTA file'
    )
    parser.add_argument(
        '--train_ratio', type=float,
        default=0.8,
        help='Ratio of data for training (default: 0.8, must be between 0 and 1)'
    )
    parser.add_argument(
        '--output_path', type=str, default='.',
        help='Output directory for split files (default: current directory)'
    )
    parser.add_argument(
        '--random_seed', type=int, default=None,
        help='Random seed for reproducibility'
    )
    args = parser.parse_args()

    # Read input files
    try:
        aa_records = list(SeqIO.parse(args.aa_fasta_path, 'fasta'))
        di_records = list(SeqIO.parse(args.di_fasta_path, 'fasta'))
    except FileNotFoundError as e:
        print(f'Error: File not found - {e}')
        return 1
    except Exception as e:
        print(f'Error reading FASTA files: {e}')
        return 1

    # Validate records
    try:
        validate_records(aa_records, di_records)
    except ValueError as e:
        print(f'Validation error: {e}')
        return 1

    # Split and write
    try:
        train_indices, valid_indices, output_files = split_and_write(
            aa_records, di_records, 
            train_ratio=args.train_ratio,
            output_dir=Path(args.output_path),
            seed=args.random_seed
        )
    except Exception as e:
        print(f'Error during split/write: {e}')
        return 1

    # Print summary
    print(f'Train size: {len(train_indices)}, Valid size: {len(valid_indices)}')
    print(f'Files created in {Path(args.output_path).resolve()}:')
    for key, path in output_files.items():
        print(f'  - {path}')

    return 0


if __name__ == '__main__':
    exit(main())
