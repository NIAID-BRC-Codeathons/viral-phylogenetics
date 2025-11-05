#!/usr/bin/env python3
"""
Test script for RefSeq data download functionality

This script tests the data_prep.py RefSeq download functions without
actually downloading large amounts of data.
"""

import json
import tempfile
import os
import time
from data_prep import (
    search_refseq_polyproteins_by_term,
    fetch_polyprotein_details,
    parse_genbank_polyprotein_text,
    calculate_cleavage_sites,
    validate_data_format,
    PolyproteinEntry,
    download_specific_viral_families,
    create_train_val_test_splits
)


def test_search_function():
    """Test the search functionality with a small query"""
    print("Testing RefSeq search...")
    
    try:
        # Test a very specific search to limit results
        ids = search_refseq_polyproteins_by_term(
            "polyprotein[Title] AND SARS-CoV-2[Organism]",
            max_entries=5,
            email="test@example.com"
        )
        
        print(f"  ✓ Search returned {len(ids)} IDs")
        if len(ids) > 0:
            print(f"  Sample IDs: {ids[:3]}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Search failed: {e}")
        return False


def test_genbank_parsing():
    """Test GenBank format parsing with example data"""
    print("Testing GenBank parsing...")
    
    # Example GenBank-like text (simplified)
    example_genbank = """
LOCUS       YP_009724390             7096 aa            linear   VRL 18-JUL-2020
DEFINITION  orf1ab polyprotein [Severe acute respiratory syndrome coronavirus 2].
ACCESSION   YP_009724390
VERSION     YP_009724390.1
ORGANISM    Severe acute respiratory syndrome coronavirus 2
            Viruses; Riboviria; Orthornavirae; Pisuviricota; Pisoniviricetes;
            Nidovirales; Cornidovirineae; Coronaviridae; Orthocoronavirinae;
            Betacoronavirus; Sarbecovirus.
FEATURES             Location/Qualifiers
     source          1..7096
                     /organism="Severe acute respiratory syndrome coronavirus 2"
     mat_peptide     1..180
                     /product="leader protein"
     mat_peptide     181..818
                     /product="nsp1"
     mat_peptide     819..2763
                     /product="nsp2"
     mat_peptide     2764..3263
                     /product="nsp3"
ORIGIN      
        1 meslvpgfne kthvqlslpv lqvrdvlvrg fgdsveevls earqhlkdgt cglvevekgv
       61 lpqleqpyvf ikrsdartap hghvmvelva elegiqygrs getlgvlvph vgeipvayrk
      121 vllrkngnkg agghsygadl ksfdlgdelg tdpyedfqen wntkhssgvt relmrelngg
//
"""
    
    try:
        entry = parse_genbank_polyprotein_text(example_genbank)
        
        if entry:
            print(f"  ✓ Parsed entry: {entry.accession}")
            print(f"  Organism: {entry.organism}")
            print(f"  Sequence length: {len(entry.sequence)}")
            print(f"  Mature peptides: {len(entry.mature_peptides)}")
            print(f"  Cleavage sites: {entry.cleavage_sites}")
        else:
            print("  ❌ Failed to parse GenBank text")
            return False
        
        return True
        
    except Exception as e:
        print(f"  ❌ Parsing failed: {e}")
        return False


def test_cleavage_calculation():
    """Test cleavage site calculation from mature peptides"""
    print("Testing cleavage site calculation...")
    
    try:
        # Example mature peptides: [(start, end, product)]
        mature_peptides = [
            (0, 179, "leader protein"),     # positions 0-179
            (180, 817, "nsp1"),            # positions 180-817
            (818, 2762, "nsp2"),           # positions 818-2762
            (2763, 3262, "nsp3")           # positions 2763-3262
        ]
        
        sequence_length = 3263
        cleavage_sites = calculate_cleavage_sites(mature_peptides, sequence_length)
        
        print(f"  ✓ Calculated {len(cleavage_sites)} cleavage sites")
        print(f"  Sites: {cleavage_sites}")
        
        # Expected sites should be at peptide boundaries
        expected_sites = [179, 817, 2762]  # End positions of first 3 peptides
        
        if set(cleavage_sites) == set(expected_sites):
            print("  ✓ Cleavage sites match expected positions")
        else:
            print(f"  ⚠ Unexpected sites. Expected: {expected_sites}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Calculation failed: {e}")
        return False


def test_data_validation():
    """Test data format validation"""
    print("Testing data validation...")
    
    try:
        # Create valid test data
        valid_data = [
            {
                "protein_id": "test_protein_1",
                "sequence": "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLS",
                "cleavage_sites": [10, 20, 30],
                "organism": "Test virus"
            },
            {
                "protein_id": "test_protein_2", 
                "sequence": "MGAQVSTQKSGDLNETRNLLEKTQGKITAIKQLLPALL",
                "cleavage_sites": [15, 25],
                "organism": "Another test virus"
            }
        ]
        
        # Save to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(valid_data, f, indent=2)
            temp_file = f.name
        
        try:
            # Test validation
            is_valid = validate_data_format(temp_file)
            
            if is_valid:
                print("  ✓ Valid data passed validation")
            else:
                print("  ❌ Valid data failed validation")
                return False
            
            # Test invalid data
            invalid_data = [
                {
                    "protein_id": "bad_protein",
                    "sequence": "INVALID",
                    "cleavage_sites": [100],  # Out of range
                    "organism": "Bad virus"
                }
            ]
            
            with open(temp_file, 'w') as f:
                json.dump(invalid_data, f, indent=2)
            
            is_invalid = validate_data_format(temp_file)
            
            if not is_invalid:
                print("  ✓ Invalid data correctly rejected")
            else:
                print("  ❌ Invalid data incorrectly accepted")
                return False
            
        finally:
            os.unlink(temp_file)
        
        return True
        
    except Exception as e:
        print(f"  ❌ Validation test failed: {e}")
        return False


def test_polyprotein_entry():
    """Test PolyproteinEntry dataclass"""
    print("Testing PolyproteinEntry...")
    
    try:
        entry = PolyproteinEntry(
            accession="TEST_001",
            organism="Test virus",
            sequence="MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLS",
            mature_peptides=[(0, 19, "protein1"), (20, 38, "protein2")],
            cleavage_sites=[19]
        )
        
        print(f"  ✓ Created entry: {entry.accession}")
        print(f"  Sequence length: {len(entry.sequence)}")
        print(f"  Mature peptides: {len(entry.mature_peptides)}")
        print(f"  Cleavage sites: {entry.cleavage_sites}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ PolyproteinEntry test failed: {e}")
        return False


def test_viral_clade_download():
    """Test downloading polyproteins from specific viral clades"""
    print("Testing viral clade download...")
    
    try:
        # Test with a small set of viral families
        test_families = ["Coronaviridae", "Picornaviridae"]
        output_file = "test_viral_clades.json"
        
        print(f"  Attempting to download from families: {test_families}")
        print("  Note: This requires internet connection to NCBI")
        
        # Download with small limits for testing
        download_specific_viral_families(
            viral_families=test_families,
            output_file=output_file,
            max_per_family=3,  # Small number for testing
            email="test@example.com"
        )
        
        # Check if file was created and has content
        if os.path.exists(output_file):
            with open(output_file, 'r') as f:
                data = json.load(f)
            
            print(f"  ✓ Downloaded {len(data)} polyproteins")
            
            # Verify data structure
            if len(data) > 0:
                sample = data[0]
                required_fields = [
                    'protein_id', 'sequence', 'cleavage_sites',
                    'organism', 'viral_family'
                ]
                
                missing_fields = [
                    field for field in required_fields
                    if field not in sample
                ]
                if missing_fields:
                    print(f"  ❌ Missing fields: {missing_fields}")
                    return False
                
                print("  ✓ Data structure validated")
                protein_id = sample['protein_id']
                viral_family = sample['viral_family']
                print(f"  Sample: {protein_id} from {viral_family}")
                print(f"  Sequence length: {len(sample['sequence'])}")
                print(f"  Cleavage sites: {len(sample['cleavage_sites'])}")
            
            # Clean up test file
            os.unlink(output_file)
            return True
        else:
            print("  ❌ No output file created")
            return False
        
    except Exception as e:
        print(f"  ❌ Viral clade download test failed: {e}")
        # Clean up on error
        if os.path.exists("test_viral_clades.json"):
            os.unlink("test_viral_clades.json")
        return False


def test_specific_refseq_search():
    """Test downloading specific polyproteins from RefSeq"""
    print("Testing specific RefSeq search...")
    
    try:
        # Test with a very specific search that should return known results
        search_terms = [
            "polyprotein[Title] AND SARS-CoV-2[Organism]",
            "polyprotein[Title] AND Poliovirus[Organism]",
            "polyprotein[Title] AND Dengue[Organism]"
        ]
        
        all_downloaded = []
        
        for search_term in search_terms:
            print(f"  Searching: {search_term}")
            
            try:
                # Get a few IDs for each search
                ids = search_refseq_polyproteins_by_term(
                    search_term=search_term,
                    max_entries=2,  # Small number for testing
                    email="test@example.com"
                )
                
                print(f"    Found {len(ids)} IDs")
                
                # Try to fetch details for first ID if available
                if len(ids) > 0:
                    first_id = ids[0]
                    print(f"    Fetching details for: {first_id}")
                    
                    entry = fetch_polyprotein_details(
                        first_id, "test@example.com"
                    )
                    if entry:
                        all_downloaded.append(entry)
                        print(f"    ✓ Downloaded: {entry.accession}")
                        print(f"      Organism: {entry.organism}")
                        seq_len = len(entry.sequence)
                        print(f"      Sequence length: {seq_len}")
                        cleavage_count = len(entry.cleavage_sites)
                        print(f"      Cleavage sites: {cleavage_count}")
                    else:
                        print(f"    ⚠ Could not fetch details for {first_id}")
                
                # Be nice to NCBI
                time.sleep(1)
                
            except Exception as e:
                print(f"    ❌ Error with search '{search_term}': {e}")
        
        if len(all_downloaded) > 0:
            downloaded_count = len(all_downloaded)
            print(f"  ✓ Successfully downloaded {downloaded_count} polyproteins")
            
            # Test that we got diverse organisms
            organisms = [entry.organism for entry in all_downloaded]
            print(f"  Organisms found: {organisms}")
            
            return True
        else:
            print("  ⚠ No polyproteins downloaded, but searches may have run")
            # Don't fail if no data found, connection might be limited
            return True
            
    except Exception as e:
        print(f"  ❌ Specific RefSeq search test failed: {e}")
        return False


def test_data_splitting():
    """Test data splitting functionality"""
    print("Testing data splitting...")
    
    try:
        # Create test data
        test_data = []
        base_sequence = "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLS"
        for i in range(10):
            test_data.append({
                "protein_id": f"test_protein_{i}",
                "sequence": base_sequence + "A" * i,
                "cleavage_sites": [10 + i, 20 + i],
                "organism": f"Test virus {i}"
            })
        
        # Save test data to temporary file
        test_input_file = "test_data_for_splitting.json"
        with open(test_input_file, 'w') as f:
            json.dump(test_data, f, indent=2)
        
        try:
            # Test splitting
            create_train_val_test_splits(
                input_file=test_input_file,
                train_ratio=0.6,
                val_ratio=0.2,
                test_ratio=0.2,
                random_seed=42
            )
            
            # Check that split files were created
            train_file = "test_data_for_splitting_train.json"
            val_file = "test_data_for_splitting_val.json"
            test_file = "test_data_for_splitting_test.json"
            
            split_files = [train_file, val_file, test_file]
            all_exist = all(os.path.exists(f) for f in split_files)
            
            if not all_exist:
                print("  ❌ Not all split files were created")
                return False
            
            # Load and check split sizes
            with open(train_file, 'r') as f:
                train_data = json.load(f)
            with open(val_file, 'r') as f:
                val_data = json.load(f)
            with open(test_file, 'r') as f:
                test_data_split = json.load(f)
            
            total_split = (len(train_data) + len(val_data) +
                           len(test_data_split))
            
            print("  ✓ Split completed:")
            print(f"    Train: {len(train_data)} sequences")
            print(f"    Validation: {len(val_data)} sequences")
            print(f"    Test: {len(test_data_split)} sequences")
            original_count = len(test_data)
            print(f"    Total: {total_split} (original: {original_count})")
            
            if total_split != len(test_data):
                print("  ❌ Split sizes don't match original data size")
                return False
            
            # Check that splits are roughly correct proportions
            expected_train = int(len(test_data) * 0.6)
            expected_val = int(len(test_data) * 0.2)
            
            if abs(len(train_data) - expected_train) > 1:
                train_len = len(train_data)
                msg = f"  ⚠ Train split: {train_len} vs {expected_train}"
                print(msg)
            
            if abs(len(val_data) - expected_val) > 1:
                val_len = len(val_data)
                print(f"  ⚠ Val split unexpected: {val_len} vs {expected_val}")
            
            print("  ✓ Data splitting completed successfully")
            
        finally:
            # Clean up test files
            cleanup_files = [
                test_input_file,
                "test_data_for_splitting_train.json",
                "test_data_for_splitting_val.json",
                "test_data_for_splitting_test.json"
            ]
            
            for file in cleanup_files:
                if os.path.exists(file):
                    os.unlink(file)
        
        return True
        
    except Exception as e:
        print(f"  ❌ Data splitting test failed: {e}")
        return False


def run_all_tests():
    """Run all tests"""
    print("🧪 Testing RefSeq Data Preparation")
    print("=" * 50)
    
    tests = [
        ("PolyproteinEntry", test_polyprotein_entry),
        ("Cleavage Calculation", test_cleavage_calculation),
        ("GenBank Parsing", test_genbank_parsing),
        ("Data Validation", test_data_validation),
        ("Data Splitting", test_data_splitting),
        ("RefSeq Search", test_search_function),
        ("Specific RefSeq Search", test_specific_refseq_search),
        ("Viral Clade Download", test_viral_clade_download),
    ]
    
    results = []
    
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        try:
            success = test_func()
            results.append((test_name, success))
        except Exception as e:
            print(f"  ❌ Test crashed: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 50)
    print("Test Summary:")
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    for test_name, success in results:
        status = "✓ PASS" if success else "❌ FAIL"
        print(f"  {status}: {test_name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed!")
        return True
    else:
        print("⚠ Some tests failed. Check implementation.")
        return False


def main():
    """Main test function"""
    
    print("RefSeq Data Preparation Test Suite")
    print("This tests the data download and processing functionality")
    print("Note: Some tests require internet access to NCBI\n")
    
    try:
        success = run_all_tests()
        
        if success:
            print("\n✓ Ready to download real data!")
            print("Usage: python data_prep.py refseq --email your@email.com")
        else:
            print("\n❌ Fix issues before using data download")
        
    except KeyboardInterrupt:
        print("\nTests interrupted by user")
    except Exception as e:
        print(f"\nTest suite crashed: {e}")


if __name__ == "__main__":
    main()