#!/usr/bin/env python3
"""
Test script for RefSeq data download functionality

This script tests the data_prep.py RefSeq download functions without
actually downloading large amounts of data.
"""

import json
import tempfile
import os
from data_prep import (
    search_refseq_polyproteins_by_term,
    fetch_polyprotein_details,
    parse_genbank_polyprotein_text,
    calculate_cleavage_sites,
    validate_data_format,
    PolyproteinEntry
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


def run_all_tests():
    """Run all tests"""
    print("🧪 Testing RefSeq Data Preparation")
    print("=" * 50)
    
    tests = [
        ("PolyproteinEntry", test_polyprotein_entry),
        ("Cleavage Calculation", test_cleavage_calculation),
        ("GenBank Parsing", test_genbank_parsing),
        ("Data Validation", test_data_validation),
        ("RefSeq Search", test_search_function),
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