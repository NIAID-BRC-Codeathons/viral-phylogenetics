# RefSeq Data Download Guide

This document explains how to download viral polyprotein data from NCBI RefSeq using the automated tools provided in this project.

## Overview

The `data_prep.py` script provides automated tools to download viral polyprotein sequences with cleavage site annotations from NCBI's RefSeq database. This eliminates the need for manual data collection and ensures reproducible datasets.

## Quick Start

### Download General Viral Polyproteins

```bash
python data_prep.py refseq \
    --output polyproteins.json \
    --max-entries 500 \
    --email your.email@institution.edu
```

### Download Specific Viral Families

```bash
python data_prep.py families \
    --families Coronaviridae Picornaviridae Flaviviridae \
    --output family_polyproteins.json \
    --max-per-family 100 \
    --email your.email@institution.edu
```

## Available Viral Families

The script includes built-in support for major virus families known to have polyproteins:

- **Coronaviridae**: SARS-CoV-2, MERS-CoV, common cold coronaviruses
- **Picornaviridae**: Poliovirus, rhinovirus, enterovirus
- **Flaviviridae**: Dengue, Zika, hepatitis C, yellow fever
- **Caliciviridae**: Norovirus, sapovirus
- **Arteriviridae**: PRRSV (porcine reproductive and respiratory syndrome virus)
- **Astroviridae**: Astrovirus
- **Hepeviridae**: Hepatitis E virus

## How It Works

### 1. Search Strategy

The script uses NCBI E-utilities to search for sequences matching:

```
polyprotein[Title] AND virus[Organism]
polyprotein[Title] AND viral[Organism]
(polyprotein OR poly-protein)[Title] AND (virus OR viral)[All Fields]
mature_protein[Feature] AND polyprotein[Title]
```

### 2. Data Extraction

For each sequence, the script:

1. Downloads the full GenBank record
2. Extracts mature peptide annotations (mat_peptide features)
3. Calculates cleavage sites from peptide boundaries
4. Validates sequence and annotation quality

### 3. Cleavage Site Detection

Cleavage sites are inferred from mature peptide coordinates:

- **Adjacent peptides**: Cleavage between consecutive mature peptides
- **Gap regions**: Cleavage at boundaries where peptides are removed
- **Processing sites**: Sites marked explicitly in annotations

## Configuration

### Required Parameters

- `--email`: Your email address (required by NCBI for API access)

### Optional Parameters

- `--max-entries`: Maximum sequences to download (default: 1000)
- `--delay`: Delay between requests in seconds (default: 0.5)
- `--output`: Output JSON file name

### Rate Limiting

The script includes automatic rate limiting to be respectful to NCBI servers:

- Default delay: 0.5 seconds between requests
- Automatic retry with exponential backoff on errors
- Progress reporting every 10 sequences

## Output Format

Downloaded data is saved in JSON format compatible with the training framework:

```json
[
  {
    "protein_id": "Coronaviridae_YP_009724390.1",
    "sequence": "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEAR...",
    "cleavage_sites": [180, 818, 2763, 3263, 3569, 3859, 3942],
    "organism": "Severe acute respiratory syndrome coronavirus 2",
    "viral_family": "Coronaviridae",
    "mature_peptides": [
      {
        "start": 0,
        "end": 179,
        "product": "leader protein"
      },
      {
        "start": 180,
        "end": 817,
        "product": "nsp1"
      }
    ]
  }
]
```

## Quality Control

### Automatic Filtering

The script automatically filters out sequences that:

- Have no mature peptide annotations
- Have malformed or incomplete sequence data
- Lack clear cleavage site information
- Are too short (<100 amino acids) or too long (>10,000 amino acids)

### Manual Validation

Always validate downloaded data:

```bash
python data_prep.py validate --input polyproteins.json
```

This checks for:
- Required fields presence
- Sequence format validity
- Cleavage site coordinates within sequence bounds
- Data type consistency

## Data Splits

After downloading, create train/validation/test splits:

```bash
python data_prep.py split \
    --input polyproteins.json \
    --train 0.7 \
    --val 0.15 \
    --test 0.15 \
    --seed 42
```

This creates three files:
- `polyproteins_train.json`
- `polyproteins_val.json`
- `polyproteins_test.json`

## Advanced Usage

### Custom Search Terms

Modify the search terms in `search_refseq_polyproteins()` for specific needs:

```python
search_terms = [
    "polyprotein[Title] AND coronavirus[Organism]",
    "polyprotein[Title] AND (SARS OR MERS)[All Fields]",
    "structural protein[Title] AND flavivirus[Organism]"
]
```

### Batch Processing

For large downloads, consider running in batches:

```bash
# Download Coronaviridae
python data_prep.py families --families Coronaviridae --max-per-family 200

# Download Picornaviridae  
python data_prep.py families --families Picornaviridae --max-per-family 200

# Combine results
python -c "
import json
data1 = json.load(open('viral_family_polyproteins.json'))
# Load and combine other files...
"
```

## Troubleshooting

### Common Issues

1. **HTTP Errors**: NCBI servers are temporarily unavailable
   - **Solution**: Wait and retry, increase delay parameter

2. **Empty Results**: Search terms don't match expected data
   - **Solution**: Check search terms, try broader queries

3. **Rate Limiting**: Too many requests too quickly
   - **Solution**: Increase `--delay` parameter

4. **Parsing Errors**: GenBank format variations
   - **Solution**: Check specific sequences manually, update parser

### Debug Mode

Enable verbose output for troubleshooting:

```python
# In data_prep.py, add debug prints
print(f"Processing {accession}: {len(genbank_text)} characters")
print(f"Found {len(mature_peptides)} mature peptides")
```

### Manual Verification

For important datasets, manually verify a sample:

1. Check sequences on NCBI website
2. Verify cleavage sites against literature
3. Compare with existing databases (UniProt, etc.)

## Best Practices

### 1. Email Configuration

Use institutional email address and follow NCBI guidelines:

```bash
export NCBI_EMAIL="researcher@university.edu"
python data_prep.py refseq --email $NCBI_EMAIL
```

### 2. Incremental Downloads

For large datasets, download incrementally:

```bash
# Day 1: Coronaviridae
python data_prep.py families --families Coronaviridae

# Day 2: Picornaviridae  
python data_prep.py families --families Picornaviridae

# Continue...
```

### 3. Version Control

Keep track of download dates and parameters:

```bash
echo "$(date): Downloaded with --max-entries 500" >> download_log.txt
```

### 4. Backup Original Data

Always keep original downloaded data:

```bash
cp polyproteins.json polyproteins_original_$(date +%Y%m%d).json
```

## Data Statistics

After downloading, analyze your dataset:

```python
import json
data = json.load(open('polyproteins.json'))

print(f"Total sequences: {len(data)}")
print(f"Total amino acids: {sum(len(d['sequence']) for d in data)}")
print(f"Total cleavage sites: {sum(len(d['cleavage_sites']) for d in data)}")

# Organism distribution
organisms = {}
for d in data:
    org = d.get('organism', 'Unknown')
    organisms[org] = organisms.get(org, 0) + 1

print("Top organisms:")
for org, count in sorted(organisms.items(), key=lambda x: x[1], reverse=True)[:10]:
    print(f"  {org}: {count}")
```

## References

- [NCBI E-utilities Documentation](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [RefSeq Database](https://www.ncbi.nlm.nih.gov/refseq/)
- [GenBank Format Specification](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)