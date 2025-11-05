#!/usr/bin/env python3
"""
Data preparation utilities for polyprotein cleavage site prediction

This script helps convert various data formats into the format expected by
the cleavage prediction framework, including downloading from RefSeq.
"""

import json
import csv
import re
import time
import urllib.request
import urllib.parse
import urllib.error
import random
import os
import itertools
from collections import defaultdict, namedtuple
from typing import List, Dict, Tuple, Optional
from pathlib import Path
from dataclasses import dataclass

# Optional imports for NCBI Datasets API
try:
    import requests
    import pandas as pd
    DATASETS_AVAILABLE = True
except ImportError:
    DATASETS_AVAILABLE = False
    print("Warning: requests/pandas not available. NCBI Datasets API features disabled.")


def safe_api_request(url: str, max_retries: int = 3, 
                     base_delay: float = 1.0) -> Optional[str]:
    """
    Make a safe API request with rate limiting and retry logic
    
    Args:
        url: URL to request
        max_retries: Maximum number of retry attempts
        base_delay: Base delay between requests (seconds)
    
    Returns:
        Response text or None if failed
    """
    for attempt in range(max_retries + 1):
        try:
            # Add jitter to avoid thundering herd
            delay = base_delay + random.uniform(0, 0.5)
            if attempt > 0:
                # Exponential backoff for retries
                delay *= (2 ** attempt)
                print(f"  Retry {attempt}/{max_retries} after {delay:.1f}s...")
            
            time.sleep(delay)
            
            with urllib.request.urlopen(url, timeout=30) as response:
                return response.read().decode()
                
        except urllib.error.HTTPError as e:
            if e.code == 429:  # Too Many Requests
                if attempt < max_retries:
                    retry_delay = base_delay * (3 ** (attempt + 1))
                    print(f"  Rate limited (429). Waiting {retry_delay:.1f}s before retry...")
                    time.sleep(retry_delay)
                    continue
                else:
                    print(f"  Rate limit exceeded after {max_retries} retries")
                    return None
            elif e.code in [502, 503, 504]:  # Server errors
                if attempt < max_retries:
                    retry_delay = base_delay * (2 ** (attempt + 1))
                    print(f"  Server error ({e.code}). Retrying in {retry_delay:.1f}s...")
                    time.sleep(retry_delay)
                    continue
                else:
                    print(f"  Server error {e.code} after {max_retries} retries")
                    return None
            else:
                print(f"  HTTP error {e.code}: {e}")
                return None
                
        except urllib.error.URLError as e:
            if attempt < max_retries:
                retry_delay = base_delay * (2 ** (attempt + 1))
                print(f"  URL error: {e}. Retrying in {retry_delay:.1f}s...")
                time.sleep(retry_delay)
                continue
            else:
                print(f"  URL error after {max_retries} retries: {e}")
                return None
                
        except Exception as e:
            print(f"  Unexpected error: {e}")
            return None
    
    return None


def parse_fasta_with_cleavage_annotations(fasta_file: str,
                                          annotation_file: str,
                                          output_file: str):
    """
    Parse FASTA sequences with separate cleavage site annotations
    
    Args:
        fasta_file: Path to FASTA file with sequences
        annotation_file: CSV file with cleavage site annotations
        output_file: Output JSON file path
    
    Expected annotation format:
    protein_id,cleavage_sites,organism
    seq1,"10,25,67",SARS-CoV-2
    seq2,"15,89",Poliovirus
    """
    
    # Read FASTA sequences
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # Take first part of header
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    # Read cleavage annotations
    annotations = {}
    with open(annotation_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            protein_id = row['protein_id']
            cleavage_sites = [int(x.strip()) for x in row['cleavage_sites'].split(',') if x.strip()]
            organism = row.get('organism', '')
            annotations[protein_id] = {
                'cleavage_sites': cleavage_sites,
                'organism': organism
            }
    
    # Combine sequences and annotations
    output_data = []
    for protein_id, sequence in sequences.items():
        if protein_id in annotations:
            output_data.append({
                'protein_id': protein_id,
                'sequence': sequence,
                'cleavage_sites': annotations[protein_id]['cleavage_sites'],
                'organism': annotations[protein_id]['organism']
            })
        else:
            print(f"Warning: No annotations found for {protein_id}")
    
    # Save to JSON
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"Processed {len(output_data)} sequences")
    print(f"Output saved to {output_file}")


def parse_genbank_polyprotein(genbank_file: str,
                              cleavage_keywords: List[str] = None,
                              output_file: str = "polyproteins.json"):
    """
    Parse GenBank files to extract polyprotein sequences and cleavage sites
    
    This is a simplified parser - you may need to adapt it for your specific
    GenBank format.
    """
    
    if cleavage_keywords is None:
        cleavage_keywords = ['cleavage', 'mature_protein', 'mat_peptide']
    
    print("Note: This is a simplified GenBank parser.")
    print("You may need to adapt it for your specific data format.")
    print(f"Looking for keywords: {cleavage_keywords}")
    
    # Placeholder implementation
    # In practice, you'd use BioPython or similar to parse GenBank files
    output_data = []
    
    # Example structure (replace with actual parsing)
    example_data = [
        {
            "protein_id": "example_from_genbank",
            "sequence": "MKQHKAMIVALIVICITAVVAALVTRKDLCEVAKLRDGGKQKD",
            "cleavage_sites": [15, 25, 35],
            "organism": "Example virus from GenBank"
        }
    ]
    
    output_data.extend(example_data)
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"Parsed GenBank file: {genbank_file}")
    print(f"Output saved to {output_file}")
    print("Note: This is example data. Implement actual GenBank parsing for real use.")


@dataclass
class PolyproteinEntry:
    """Data structure for a polyprotein with cleavage information"""
    accession: str
    organism: str
    sequence: str
    mature_peptides: List[Tuple[int, int, str]]  # (start, end, product_name)
    cleavage_sites: List[int]  # Derived from mature peptides


# NCBI Datasets API configuration
DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
DATASETS_HEADERS = {"Accept": "application/json"}
if os.getenv("NCBI_API_KEY"):
    DATASETS_HEADERS["X-API-Key"] = os.getenv("NCBI_API_KEY")


def stream_annotation_report_by_taxon(taxon: str, params: Dict = None, 
                                      backoff: float = 2.0):
    """
    Yields Python dicts for each report from NCBI Datasets API
    /virus/taxon/{taxon}/annotation_report.
    
    Args:
        taxon: Viral family/taxon name or ID
        params: Additional query parameters
        backoff: Delay on JSON decode errors
    
    Yields:
        Dict records from the annotation report
    """
    if not DATASETS_AVAILABLE:
        raise RuntimeError("requests/pandas required for Datasets API")
    
    url = f"{DATASETS_API}/virus/taxon/{taxon}/annotation_report"
    
    # Add pagination parameters
    if params is None:
        params = {}
    
    # Start with first page
    page_token = None
    
    while True:
        current_params = params.copy()
        if page_token:
            current_params['page_token'] = page_token
        
        try:
            response = requests.get(url, headers=DATASETS_HEADERS, 
                                   params=current_params, timeout=120)
            
            if response.status_code != 200:
                raise RuntimeError(f"HTTP {response.status_code} for {url}: "
                                  f"{response.text[:500]}")
            
            data = response.json()
            
            # Yield each report in the reports array
            reports = data.get('reports', [])
            for report in reports:
                yield report
            
            # Check for next page
            page_token = data.get('next_page_token')
            if not page_token:
                break
                
            # Be nice to the API
            time.sleep(backoff)
            
        except requests.RequestException as e:
            raise RuntimeError(f"Request failed for {url}: {e}")
        except json.JSONDecodeError as e:
            raise RuntimeError(f"Invalid JSON response from {url}: {e}")


def collect_mat_peptides(record: Dict) -> List[Dict]:
    """
    From one NCBI Datasets API record, extract mature peptide information.
    
    Args:
        record: Virus record from Datasets API
    
    Returns:
        List of mature peptide feature dictionaries
    """
    out = []
    
    # Get record-level info
    accession = record.get("accession", "Unknown")
    genes = record.get("genes", [])
    
    for gene in genes:
        gene_name = gene.get("name", "Unknown")
        cds_list = gene.get("cds", [])
        
        for cds in cds_list:
            cds_name = cds.get("name", "Unknown")
            
            # Get nucleotide coordinates
            nucleotide = cds.get("nucleotide", {})
            ranges = nucleotide.get("range", [])
            
            # Get protein info
            protein = cds.get("protein", {})
            protein_name = protein.get("name", "Unknown")
            
            # Extract start/end from ranges (use first range if multiple)
            start, end = None, None
            if ranges and isinstance(ranges, list) and len(ranges) > 0:
                first_range = ranges[0]
                if isinstance(first_range, dict):
                    start = int(first_range.get("begin", 0))
                    end = int(first_range.get("end", 0))
            
            # EXPANDED INCLUSION CRITERIA - be more permissive
            # Include viral proteins that could be relevant for cleavage analysis
            include_keywords = [
                'polyprotein', 'protein', 'peptide', 'mature',
                'structural', 'non-structural', 'nonstructural', 
                'orf', 'open reading frame', 'glycoprotein',
                'capsid', 'envelope', 'replicase', 'spike',
                'nucleocapsid', 'membrane'
            ]
            
            # Also include if it's a significant protein (length > 50 aa)
            protein_length = 0
            if start and end:
                protein_length = abs(end - start) // 3  # Rough amino acid length
            
            should_include = (
                any(keyword in cds_name.lower() for keyword in include_keywords) or
                any(keyword in gene_name.lower() for keyword in include_keywords) or
                protein_length > 50  # Include larger proteins regardless of name
            )
            
            if should_include:
                out.append({
                    "nuc_accession": accession,
                    "gene_name": gene_name,
                    "cds_name": cds_name,
                    "protein_name": protein_name,
                    "start": start,
                    "end": end,
                    "length_aa": protein_length,
                    "organism": "Unknown",  # Not in this API response
                    "feature_type": "cds",
                    "product": cds_name,
                    "parent_name": (gene_name if "polyprotein" in gene_name.lower() 
                                   else cds_name),
                    "seq_id": accession,
                    "strand": None,  # Not provided
                    "protein_acc": protein.get("accession_version"),
                    "assembly_accession": None,
                    "tax_id": None
                })
    
    return out


def to_cleavage_junctions(mat_rows: List[Dict]) -> List[Dict]:
    """
    Group mat_peptides by (assembly_accession, parent_name, seq_id, strand),
    sort by genomic coordinate, and emit junctions between consecutive peptides.
    
    Args:
        mat_rows: List of mature peptide dictionaries
    
    Returns:
        List of cleavage junction dictionaries
    """
    Junction = namedtuple("Junction",
        "assembly_accession tax_id organism nuc_accession parent_name "
        "left_product right_product left_end right_start junction_genomic "
        "seq_id strand"
    )
    
    cleavages = []
    
    # Key that approximates a polyprotein context
    def keyfunc(d):
        return (d["assembly_accession"], d["parent_name"], 
                d["seq_id"], d["strand"])
    
    # Sort and group by polyprotein context
    sorted_rows = sorted(
        mat_rows, 
        key=lambda d: (d["assembly_accession"], d["parent_name"] or "", 
                      d["seq_id"] or "", d["start"] or 0)
    )
    
    for key, group in itertools.groupby(sorted_rows, key=keyfunc):
        group = list(group)
        
        # Sort by strand-aware coordinate
        if (group[0]["strand"] or "").lower() == "minus":
            group.sort(key=lambda d: (d["start"] or 0), reverse=True)
        else:
            group.sort(key=lambda d: (d["start"] or 0))
        
        # Create junctions between consecutive peptides
        for a, b in zip(group, group[1:]):
            j = a["end"] if (a["strand"] or "").lower() != "minus" else a["start"]
            
            cleavages.append(Junction(
                a["assembly_accession"], a["tax_id"], a["organism"], 
                a["nuc_accession"], a["parent_name"],
                a["product"], b["product"],
                a["end"], b["start"], j, a["seq_id"], a["strand"]
            )._asdict())
    
    return cleavages


def download_viral_families_datasets(
    families: List[str] = None,
    output_file: str = "datasets_polyproteins.json",
    max_records_per_family: int = 100
) -> List[Dict]:
    """
    Download viral polyprotein data using NCBI Datasets API
    
    Args:
        families: List of viral family names
        output_file: Output JSON file path
        max_records_per_family: Maximum records to process per family
    
    Returns:
        List of polyprotein dictionaries in standard format
    """
    if not DATASETS_AVAILABLE:
        raise RuntimeError("requests and pandas required for Datasets API")
    
    if families is None:
        families = ["Flaviviridae", "Picornaviridae", "Coronaviridae", 
                   "Togaviridae", "Caliciviridae"]
    
    print(f"🦠 Downloading from {len(families)} viral families using "
          f"NCBI Datasets API...")
    
    all_mats = []
    
    for family in families:
        print(f"  Processing {family}...")
        try:
            record_count = 0
            for record in stream_annotation_report_by_taxon(family):
                mat_peptides = collect_mat_peptides(record)
                if mat_peptides:  # Only process records with peptides
                    all_mats.extend(mat_peptides)
                    print(f"    Record {record_count+1}: Found {len(mat_peptides)} peptides")
                
                record_count += 1
                
                # Limit processing for testing/speed
                if record_count >= max_records_per_family:
                    print(f"    Reached limit of {max_records_per_family} records")
                    break
                    
        except Exception as e:
            print(f"  Error processing {family}: {e}")
            continue
    
    print(f"✓ Collected {len(all_mats)} mature peptide features")
    
    # Build cleavage junctions per polyprotein
    cleavages = to_cleavage_junctions(all_mats)
    print(f"✓ Found {len(cleavages)} cleavage sites")
    
    # Convert to standard format (group by polyprotein)
    polyprotein_data = defaultdict(lambda: {
        'cleavage_sites': [],
        'mature_peptides': [],
        'organism': None,
        'assembly_accession': None
    })
    
    # Group cleavages by parent polyprotein
    for cleavage in cleavages:
        key = (cleavage['nuc_accession'], cleavage['parent_name'])
        
        polyprotein_data[key]['cleavage_sites'].append(
            cleavage['junction_genomic']
        )
        polyprotein_data[key]['organism'] = cleavage['organism']
        polyprotein_data[key]['assembly_accession'] = cleavage.get('assembly_accession')
    
    # Convert to final format
    result_data = []
    for (nuc_acc, parent_name), data in polyprotein_data.items():
        if not data['cleavage_sites']:
            continue
            
        result_data.append({
            "protein_id": f"datasets_{nuc_acc}_{parent_name}",
            "sequence": "",  # Would need additional API call for sequence
            "cleavage_sites": sorted(set(data['cleavage_sites'])),
            "organism": data['organism'],
            "viral_family": "Unknown",  # Could be inferred from organism
            "assembly_accession": data['assembly_accession'],
            "parent_name": parent_name,
            "nucleotide_accession": nuc_acc,
            "source": "ncbi_datasets"
        })
    
    print(f"✓ Processed {len(result_data)} polyproteins")
    
    # Save to file
    with open(output_file, 'w') as f:
        json.dump(result_data, f, indent=2)
    
    print(f"✅ Saved data to {output_file}")
    return result_data


def download_refseq_viral_polyproteins(
    output_file: str = "refseq_polyproteins.json",
    max_entries: int = 1000,
    delay: float = 0.5,
    email: str = "dmoi@unil.ch"
) -> List[PolyproteinEntry]:
    """
    Download viral polyprotein data from RefSeq with cleavage annotations
    
    Args:
        output_file: Output JSON file path
        max_entries: Maximum number of entries to download
        delay: Delay between requests (seconds) to be respectful to NCBI
        email: Email for NCBI requests (required)
    
    Returns:
        List of PolyproteinEntry objects
    """
    
    print("Downloading viral polyprotein data from RefSeq...")
    
    # Step 1: Search for viral polyproteins
    search_results = search_refseq_polyproteins(max_entries, email)
    
    # Step 2: Download detailed records for each hit
    polyproteins = []
    
    for i, accession in enumerate(search_results):
        print(f"Processing {i+1}/{len(search_results)}: {accession}")
        
        try:
            entry = fetch_polyprotein_details(accession, email)
            if entry and len(entry.cleavage_sites) > 0:
                polyproteins.append(entry)
                print(f"  Found {len(entry.cleavage_sites)} cleavage sites")
            else:
                print("  No cleavage sites found, skipping")
        
        except Exception as e:
            print(f"  Error processing {accession}: {e}")
        
        # Be respectful to NCBI servers
        time.sleep(delay)
    
    print(f"\nDownloaded {len(polyproteins)} polyproteins with cleavage data")
    
    # Step 3: Convert to our standard format and save
    output_data = []
    for entry in polyproteins:
        output_data.append({
            "protein_id": entry.accession,
            "sequence": entry.sequence,
            "cleavage_sites": entry.cleavage_sites,
            "organism": entry.organism,
            "mature_peptides": [
                {
                    "start": start,
                    "end": end,
                    "product": product
                }
                for start, end, product in entry.mature_peptides
            ]
        })
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"Data saved to {output_file}")
    return polyproteins


def search_refseq_polyproteins(max_entries: int = 1000,
                               email: str = "dmoi@unil.ch") -> List[str]:
    """
    Search RefSeq for viral polyprotein sequences
    
    Returns list of accession numbers
    """
    
    # NCBI E-utilities search for viral polyproteins
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    # Search terms for viral polyproteins
    polyprotein_genes = ("( polyprotein[Protein Name] OR ORF1a[Gene] OR "
                         "ORF1ab[Gene] OR pp1a[Gene] OR pp1ab[Gene] )")
    
    search_terms = [
        "polyprotein[Title] AND virus[Organism]",
        "polyprotein[Title] AND viral[Organism]",
        "(polyprotein OR poly-protein)[Title] AND (virus OR viral)[All Fields]",
        f"{polyprotein_genes} AND virus[Organism]",
        f"{polyprotein_genes} AND viral[Organism]",
        "mature_protein[Feature] AND polyprotein[Title]",
        ("( replicase[Protein Name] OR nonstructural[Protein Name] ) AND "
         "polyprotein[Title] AND virus[Organism]")
    ]
    
    all_ids = set()
    
    for search_term in search_terms:
        print(f"Searching: {search_term}")
        
        # E-search to get IDs
        search_url = f"{base_url}esearch.fcgi"
        params = {
            'db': 'protein',
            'term': search_term,
            'retmax': max_entries // len(search_terms),
            'retmode': 'json',
            'sort': 'relevance',
            'email': email
        }
        
        try:
            url = f"{search_url}?{urllib.parse.urlencode(params)}"
            response_text = safe_api_request(url)
            
            if response_text:
                data = json.loads(response_text)
                
                if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                    ids = data['esearchresult']['idlist']
                    all_ids.update(ids)
                    print(f"  Found {len(ids)} entries")
                else:
                    print("  No results found")
            else:
                print("  Failed to get response")
        
        except Exception as e:
            print(f"  Error in search: {e}")
        
        # Additional delay between search terms
        time.sleep(1.0)  # Default search delay
    
    print(f"Total unique entries found: {len(all_ids)}")
    return list(all_ids)[:max_entries]


def fetch_polyprotein_details(protein_id: str,
                              email: str = "user@example.com") -> Optional[PolyproteinEntry]:
    """
    Fetch detailed information for a specific protein ID
    
    Args:
        protein_id: NCBI protein ID or accession
        email: Email for NCBI requests
    
    Returns:
        PolyproteinEntry object or None if not suitable
    """
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    # Fetch full record in GenBank format
    fetch_url = f"{base_url}efetch.fcgi"
    params = {
        'db': 'protein',
        'id': protein_id,
        'rettype': 'gb',
        'retmode': 'text',
        'email': email
    }
    
    try:
        url = f"{fetch_url}?{urllib.parse.urlencode(params)}"
        genbank_text = safe_api_request(url, max_retries=5, base_delay=1.5)
        
        if genbank_text:
            return parse_genbank_polyprotein_text(genbank_text)
        else:
            print(f"Failed to fetch data for {protein_id}")
            return None
    
    except Exception as e:
        print(f"Error fetching {protein_id}: {e}")
        return None


def parse_genbank_polyprotein_text(genbank_text: str) -> Optional[PolyproteinEntry]:
    """
    Parse GenBank format text to extract polyprotein information
    
    Args:
        genbank_text: Raw GenBank format text
    
    Returns:
        PolyproteinEntry object or None
    """
    
    lines = genbank_text.strip().split('\n')
    
    # Extract basic information
    accession = None
    organism = None
    sequence = []
    mature_peptides = []
    
    # Parse header information
    for line in lines:
        if line.startswith('ACCESSION'):
            accession = line.split()[1]
        elif line.startswith('  ORGANISM'):
            organism = line.replace('  ORGANISM  ', '').strip()
    
    if not accession:
        return None
    
    # Parse features for mature peptides
    in_features = False
    current_feature = None
    
    for line in lines:
        if line.startswith('FEATURES'):
            in_features = True
            continue
        elif line.startswith('ORIGIN'):
            break
        elif not in_features:
            continue
        
        # Look for mature_protein or mat_peptide features
        if re.match(r'\s+(mat_peptide|mature_protein)', line):
            # Extract coordinates
            coord_match = re.search(r'(\d+)\.\.(\d+)', line)
            if coord_match:
                start = int(coord_match.group(1)) - 1  # Convert to 0-indexed
                end = int(coord_match.group(2)) - 1
                current_feature = {'start': start, 'end': end, 'product': ''}
        
        elif current_feature and '/product=' in line:
            # Extract product name
            product_match = re.search(r'/product="([^"]+)"', line)
            if product_match:
                current_feature['product'] = product_match.group(1)
                mature_peptides.append((
                    current_feature['start'],
                    current_feature['end'],
                    current_feature['product']
                ))
                current_feature = None
    
    # Extract sequence
    in_origin = False
    for line in lines:
        if line.startswith('ORIGIN'):
            in_origin = True
            continue
        elif line.startswith('//'):
            break
        elif in_origin:
            # Remove line numbers and spaces
            seq_part = re.sub(r'[^a-zA-Z]', '', line)
            sequence.append(seq_part.upper())
    
    if not sequence or not mature_peptides:
        return None
    
    full_sequence = ''.join(sequence)
    
    # Calculate cleavage sites from mature peptides
    cleavage_sites = calculate_cleavage_sites(mature_peptides, len(full_sequence))
    
    return PolyproteinEntry(
        accession=accession,
        organism=organism or "Unknown virus",
        sequence=full_sequence,
        mature_peptides=mature_peptides,
        cleavage_sites=cleavage_sites
    )


def calculate_cleavage_sites(mature_peptides: List[Tuple[int, int, str]],
                             sequence_length: int) -> List[int]:
    """
    Calculate cleavage sites from mature peptide coordinates
    
    Cleavage occurs where there's a gap between consecutive mature peptides
    or at the boundaries where peptides are removed during processing.
    
    Args:
        mature_peptides: List of (start, end, product_name) tuples
        sequence_length: Length of the full polyprotein sequence
    
    Returns:
        List of 0-indexed cleavage site positions
    """
    
    if not mature_peptides:
        return []
    
    # Sort by start position
    sorted_peptides = sorted(mature_peptides, key=lambda x: x[0])
    
    cleavage_sites = []
    
    # Check for cleavage between consecutive peptides
    for i in range(len(sorted_peptides) - 1):
        current_end = sorted_peptides[i][1]
        next_start = sorted_peptides[i + 1][0]
        
        # If there's a gap, there might be cleavage sites
        if next_start > current_end + 1:
            # Cleavage at the end of current peptide
            cleavage_sites.append(current_end)
            # Cleavage at the start of next peptide
            cleavage_sites.append(next_start)
        elif next_start == current_end + 1:
            # Adjacent peptides - cleavage between them
            cleavage_sites.append(current_end)
    
    # Remove duplicates and sort
    cleavage_sites = sorted(list(set(cleavage_sites)))
    
    # Filter out invalid positions
    cleavage_sites = [site for site in cleavage_sites
                      if 0 <= site < sequence_length - 1]
    
    return cleavage_sites


def download_specific_viral_families(
    viral_families: List[str] = None,
    output_file: str = "viral_family_polyproteins.json",
    max_per_family: int = 100,
    email: str = "user@example.com",
    delay_config: Dict[str, float] = None
) -> None:
    """
    Download polyproteins from specific viral families
    
    Args:
        viral_families: List of viral family names
        output_file: Output file path
        max_per_family: Maximum entries per family
        email: Email for NCBI requests
        delay_config: Custom delay configuration for rate limiting
    """
    
    # Default delay configuration for NCBI rate limiting
    if delay_config is None:
        delay_config = {
            'base_delay': 2.0,      # Base delay between requests
            'search_delay': 1.0,    # Delay between search terms
            'batch_1_delay': 2.0,   # First 10 proteins
            'batch_2_delay': 3.0,   # Next 10 proteins  
            'batch_3_delay': 4.0,   # Remaining proteins
            'error_delay': 5.0      # Delay after errors
        }
    
    if viral_families is None:
        viral_families = [
            'Coronaviridae',    # SARS-CoV-2, MERS, etc.
            'Picornaviridae',   # Poliovirus, rhinovirus
            'Flaviviridae',     # Dengue, Zika, HCV
            'Caliciviridae',    # Norovirus
            'Arteriviridae',    # PRRSV
            'Astroviridae',     # Astrovirus
            'Hepeviridae'       # Hepatitis E
        ]
    
    all_data = []
    
    for family in viral_families:
        print(f"\nProcessing {family}...")
        
        # Try multiple search strategies to maximize results
        search_strategies = [
            f"polyprotein[Title] AND {family}[Organism]",
            f"ORF1ab[Gene] AND {family}[Organism]",
            f"ORF1a[Gene] AND {family}[Organism]",
            f"polyprotein[Protein Name] AND {family}[Organism]",
            f"replicase[Protein Name] AND {family}[Organism]"
        ]
        
        all_family_ids = set()
        for strategy in search_strategies:
            try:
                max_per_strategy = max_per_family // len(search_strategies) + 1
                ids = search_refseq_polyproteins_by_term(
                    strategy, max_per_strategy, email)
                all_family_ids.update(ids)
                if len(ids) > 0:
                    print(f"  Strategy '{strategy}' found {len(ids)} entries")
            except Exception as e:
                print(f"  Strategy '{strategy}' failed: {e}")
        
        family_ids = list(all_family_ids)[:max_per_family]
        print(f"  Total unique entries found: {len(family_ids)}")
        
        family_polyproteins = []
        print(family_ids[0:10] , '...')
        for i, protein_id in enumerate(family_ids):
            try:
                print(f"  Fetching protein {i+1}/{len(family_ids)}: "
                      f"{protein_id}")
                entry = fetch_polyprotein_details(protein_id, email)
                if entry and len(entry.cleavage_sites) > 0:
                    family_polyproteins.append(entry)
                    print(f"    ✓ Found {len(entry.cleavage_sites)} "
                          f"cleavage sites")
                else:
                    print("    ✗ No cleavage sites found")
                
                # Progressive delay using configuration
                if i < 10:
                    time.sleep(delay_config['batch_1_delay'])
                elif i < 20:
                    time.sleep(delay_config['batch_2_delay'])
                else:
                    time.sleep(delay_config['batch_3_delay'])
                    
            except Exception as e:
                print(f"Error processing {protein_id}: {e}")
                time.sleep(delay_config['error_delay'])
        
        count = len(family_polyproteins)
        print(f"Downloaded {count} polyproteins from {family}")
        
        # Convert to standard format
        for entry in family_polyproteins:
            all_data.append({
                "protein_id": f"{family}_{entry.accession}",
                "sequence": entry.sequence,
                "cleavage_sites": entry.cleavage_sites,
                "organism": entry.organism,
                "viral_family": family,
                "mature_peptides": [
                    {
                        "start": start,
                        "end": end,
                        "product": product
                    }
                    for start, end, product in entry.mature_peptides
                ]
            })
    
    # Save combined data
    with open(output_file, 'w') as f:
        json.dump(all_data, f, indent=2)
    
    print(f"\nTotal downloaded: {len(all_data)} polyproteins")
    print(f"Data saved to {output_file}")


def search_refseq_polyproteins_by_term(search_term: str,
                                       max_entries: int,
                                       email: str) -> List[str]:
    """Helper function to search by specific term"""
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    search_url = f"{base_url}esearch.fcgi"
    
    params = {
        'db': 'protein',
        'term': search_term,
        'retmax': max_entries,
        'retmode': 'json',
        'sort': 'relevance',
        'email': email
    }
    
    try:
        url = f"{search_url}?{urllib.parse.urlencode(params)}"
        response_text = safe_api_request(url, max_retries=3, base_delay=1.0)
        
        if response_text:
            data = json.loads(response_text)
            
            if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                return data['esearchresult']['idlist']
        else:
            print(f"Failed to get response for search: {search_term}")
    except Exception as e:
        print(f"Search error: {e}")
    
    return []


def create_train_val_test_splits(input_file: str,
                                 train_ratio: float = 0.7,
                                 val_ratio: float = 0.15,
                                 test_ratio: float = 0.15,
                                 random_seed: int = 42):
    """
    Split data into train/validation/test sets
    
    Args:
        input_file: JSON file with all data
        train_ratio: Fraction for training set
        val_ratio: Fraction for validation set
        test_ratio: Fraction for test set
        random_seed: Random seed for reproducibility
    """
    
    import random
    random.seed(random_seed)
    
    # Load data
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    # Shuffle data
    random.shuffle(data)
    
    # Calculate split indices
    n_total = len(data)
    n_train = int(n_total * train_ratio)
    n_val = int(n_total * val_ratio)
    
    # Split data
    train_data = data[:n_train]
    val_data = data[n_train:n_train + n_val]
    test_data = data[n_train + n_val:]
    
    # Save splits
    base_name = Path(input_file).stem
    
    train_file = f"{base_name}_train.json"
    val_file = f"{base_name}_val.json"
    test_file = f"{base_name}_test.json"
    
    with open(train_file, 'w') as f:
        json.dump(train_data, f, indent=2)
    
    with open(val_file, 'w') as f:
        json.dump(val_data, f, indent=2)
    
    with open(test_file, 'w') as f:
        json.dump(test_data, f, indent=2)
    
    print(f"Data split completed:")
    print(f"  Training: {len(train_data)} sequences -> {train_file}")
    print(f"  Validation: {len(val_data)} sequences -> {val_file}")
    print(f"  Test: {len(test_data)} sequences -> {test_file}")


def validate_data_format(json_file: str):
    """
    Validate that data is in the correct format for training
    """
    
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error loading JSON file: {e}")
        return False
    
    if not isinstance(data, list):
        print("Error: Data should be a list of dictionaries")
        return False
    
    required_fields = ['protein_id', 'sequence', 'cleavage_sites']
    
    for i, item in enumerate(data):
        if not isinstance(item, dict):
            print(f"Error: Item {i} is not a dictionary")
            return False
        
        for field in required_fields:
            if field not in item:
                print(f"Error: Item {i} missing required field '{field}'")
                return False
        
        # Validate sequence
        if not isinstance(item['sequence'], str):
            print(f"Error: Item {i} sequence is not a string")
            return False
        
        if not item['sequence']:
            print(f"Error: Item {i} has empty sequence")
            return False
        
        # Validate cleavage sites
        if not isinstance(item['cleavage_sites'], list):
            print(f"Error: Item {i} cleavage_sites is not a list")
            return False
        
        for site in item['cleavage_sites']:
            if not isinstance(site, int):
                print(f"Error: Item {i} has non-integer cleavage site: {site}")
                return False
            
            if site < 0 or site >= len(item['sequence']):
                print(f"Error: Item {i} cleavage site {site} out of range for sequence length {len(item['sequence'])}")
                return False
    
    print(f"Data validation passed! {len(data)} sequences are correctly formatted.")
    return True


def download_viral_polyproteins_unified(
    method: str = "refseq",
    families: List[str] = None,
    output_file: str = "unified_polyproteins.json",
    max_per_family: int = 100,
    email: str = "user@example.com",
    **kwargs
) -> List[Dict]:
    """
    Unified interface for downloading viral polyprotein data using either
    RefSeq E-utilities or NCBI Datasets API.
    
    Args:
        method: "refseq" or "datasets" 
        families: List of viral family names
        output_file: Output JSON file path
        max_per_family: Maximum entries per family
        email: Email for NCBI requests
        **kwargs: Additional method-specific parameters
    
    Returns:
        List of polyprotein dictionaries
    """
    
    if method.lower() == "datasets":
        if not DATASETS_AVAILABLE:
            print("❌ Datasets API requires 'requests' and 'pandas' packages")
            print("   Install with: pip install requests pandas")
            return []
        
        print("🔬 Using NCBI Datasets API (faster, structured data)...")
        return download_viral_families_datasets(families, output_file)
        
    elif method.lower() == "refseq":
        print("🧬 Using RefSeq E-utilities (comprehensive, slower)...")
        
        # Use existing RefSeq download function
        if families is None:
            families = [
                'Coronaviridae', 'Flaviviridae', 'Picornaviridae',
                'Caliciviridae', 'Togaviridae'
            ]
        
        delay_config = kwargs.get('delay_config', {
            'base_delay': 2.0,
            'search_delay': 1.0,
            'batch_1_delay': 2.0,
            'batch_2_delay': 3.0,
            'batch_3_delay': 4.0,
            'error_delay': 5.0
        })
        
        download_specific_viral_families(
            families, output_file, max_per_family, email, delay_config
        )
        
        # Load and return the data
        try:
            with open(output_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            print(f"❌ Output file {output_file} not found")
            return []
            
    else:
        raise ValueError(f"Unknown method: {method}. Use 'refseq' or 'datasets'")


def compare_download_methods(families: List[str] = None, 
                           sample_size: int = 5) -> Dict:
    """
    Compare the two download methods on a small sample.
    
    Args:
        families: Viral families to test
        sample_size: Maximum entries per family for testing
    
    Returns:
        Comparison statistics
    """
    if families is None:
        families = ['Coronaviridae']  # Small test
    
    print("🔬 Comparing download methods...")
    
    results = {
        'refseq': {'data': [], 'time': 0, 'error': None},
        'datasets': {'data': [], 'time': 0, 'error': None}
    }
    
    # Test RefSeq method
    print("\n📡 Testing RefSeq E-utilities...")
    start_time = time.time()
    try:
        results['refseq']['data'] = download_viral_polyproteins_unified(
            method="refseq",
            families=families,
            output_file="test_refseq.json",
            max_per_family=sample_size,
            email="test@example.com"
        )
        results['refseq']['time'] = time.time() - start_time
        print(f"✓ RefSeq: {len(results['refseq']['data'])} polyproteins "
              f"in {results['refseq']['time']:.1f}s")
    except Exception as e:
        results['refseq']['error'] = str(e)
        print(f"❌ RefSeq failed: {e}")
    
    # Test Datasets method  
    print("\n🗄️  Testing NCBI Datasets API...")
    start_time = time.time()
    try:
        results['datasets']['data'] = download_viral_polyproteins_unified(
            method="datasets", 
            families=families,
            output_file="test_datasets.json"
        )
        results['datasets']['time'] = time.time() - start_time
        print(f"✓ Datasets: {len(results['datasets']['data'])} polyproteins "
              f"in {results['datasets']['time']:.1f}s")
    except Exception as e:
        results['datasets']['error'] = str(e)
        print(f"❌ Datasets failed: {e}")
    
    # Summary
    print(f"\n📊 COMPARISON SUMMARY:")
    for method, result in results.items():
        if result['error']:
            print(f"   {method}: FAILED - {result['error']}")
        else:
            count = len(result['data'])
            time_taken = result['time']
            print(f"   {method}: {count} polyproteins in {time_taken:.1f}s")
    
    # Cleanup
    for file in ['test_refseq.json', 'test_datasets.json']:
        if os.path.exists(file):
            os.remove(file)
    
    return results


def main():
    """Command line interface for data preparation"""
    
    import argparse
    
    parser = argparse.ArgumentParser(description="Data preparation for cleavage prediction")
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # FASTA + annotations parser
    fasta_parser = subparsers.add_parser('fasta', help='Parse FASTA with annotation file')
    fasta_parser.add_argument('--fasta', required=True, help='FASTA file with sequences')
    fasta_parser.add_argument('--annotations', required=True, help='CSV file with cleavage annotations')
    fasta_parser.add_argument('--output', default='polyproteins.json', help='Output JSON file')
    
    # GenBank parser
    genbank_parser = subparsers.add_parser('genbank', help='Parse GenBank file')
    genbank_parser.add_argument('--input', required=True, help='GenBank file')
    genbank_parser.add_argument('--output', default='polyproteins.json', help='Output JSON file')
    
    # RefSeq downloader
    refseq_parser = subparsers.add_parser('refseq', help='Download from RefSeq')
    refseq_parser.add_argument('--output', default='refseq_polyproteins.json', help='Output JSON file')
    refseq_parser.add_argument('--max-entries', type=int, default=1000, help='Maximum entries to download')
    refseq_parser.add_argument('--delay', type=float, default=0.5, help='Delay between requests (seconds)')
    refseq_parser.add_argument('--email', default='user@example.com', help='Email for NCBI requests')
    
    # Viral families downloader
    families_parser = subparsers.add_parser('families', help='Download specific viral families')
    families_parser.add_argument('--families', nargs='+', help='Viral family names')
    families_parser.add_argument('--output', default='viral_family_polyproteins.json', help='Output JSON file')
    families_parser.add_argument('--max-per-family', type=int, default=100, help='Max entries per family')
    families_parser.add_argument('--email', default='user@example.com', help='Email for NCBI requests')
    
    # Unified downloader (NEW)
    unified_parser = subparsers.add_parser('download', help='Unified download interface')
    unified_parser.add_argument('--method', choices=['refseq', 'datasets'], default='refseq',
                                help='Download method to use')
    unified_parser.add_argument('--families', nargs='+', help='Viral family names')
    unified_parser.add_argument('--output', default='unified_polyproteins.json', help='Output JSON file')
    unified_parser.add_argument('--max-per-family', type=int, default=100, help='Max entries per family')
    unified_parser.add_argument('--email', default='user@example.com', help='Email for NCBI requests')
    
    # Method comparison (NEW)
    compare_parser = subparsers.add_parser('compare', help='Compare download methods')
    compare_parser.add_argument('--families', nargs='+', help='Viral family names to test')
    compare_parser.add_argument('--sample-size', type=int, default=5, help='Max entries per family for testing')
    
    # Data splitter
    split_parser = subparsers.add_parser('split', help='Split data into train/val/test')
    split_parser.add_argument('--input', required=True, help='Input JSON file')
    split_parser.add_argument('--train', type=float, default=0.7, help='Training fraction')
    split_parser.add_argument('--val', type=float, default=0.15, help='Validation fraction')
    split_parser.add_argument('--test', type=float, default=0.15, help='Test fraction')
    split_parser.add_argument('--seed', type=int, default=42, help='Random seed')
    
    # Data validator
    validate_parser = subparsers.add_parser('validate', help='Validate data format')
    validate_parser.add_argument('--input', required=True, help='JSON file to validate')
    
    args = parser.parse_args()
    
    if args.command == 'fasta':
        parse_fasta_with_cleavage_annotations(
            args.fasta,
            args.annotations,
            args.output
        )
    elif args.command == 'genbank':
        parse_genbank_polyprotein(args.input, output_file=args.output)
    elif args.command == 'refseq':
        download_refseq_viral_polyproteins(
            args.output,
            args.max_entries,
            args.delay,
            args.email
        )
    elif args.command == 'families':
        download_specific_viral_families(
            args.families,
            args.output,
            args.max_per_family,
            args.email
        )
    elif args.command == 'download':
        download_viral_polyproteins_unified(
            method=args.method,
            families=args.families,
            output_file=args.output,
            max_per_family=args.max_per_family,
            email=args.email
        )
    elif args.command == 'compare':
        compare_download_methods(args.families, args.sample_size)
    elif args.command == 'split':
        if abs(args.train + args.val + args.test - 1.0) > 1e-6:
            print("Error: Train, validation, and test ratios must sum to 1.0")
            return
        create_train_val_test_splits(
            args.input,
            args.train,
            args.val,
            args.test,
            args.seed
        )
    elif args.command == 'validate':
        validate_data_format(args.input)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()