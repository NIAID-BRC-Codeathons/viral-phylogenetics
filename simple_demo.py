#!/usr/bin/env python3
"""
Simple demonstration of polyprotein cleavage site prediction

This script shows basic usage of the framework with example data.
"""

import json
import torch
from esmretrain import PolyproteinCleavagePredictor

def create_example_data():
    """Create some example polyprotein data for demonstration"""
    
    # Example sequences with known cleavage sites
    # These are simplified examples - real polyproteins are much longer
    example_data = [
        {
            "protein_id": "example_coronavirus",
            "sequence": "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVA",
            "cleavage_sites": [15, 35, 55, 75],
            "organism": "Example coronavirus"
        },
        {
            "protein_id": "example_picornavirus", 
            "sequence": "MGAQVSTQKSGDLNETRNLLEKTQGKITAIKQLLPALLQRGLTSKDPFGKATLRTLVLAKQ",
            "cleavage_sites": [20, 40],
            "organism": "Example picornavirus"
        },
        {
            "protein_id": "example_flavivirus",
            "sequence": "MRCVGVGNRDFVEGLSGATWVDVVLEHGGCVTMAQGKPSLDIKETDVFCVKVLAPYMPGVVQKQTCEQCRRDLKFEEGFDVFRQCLHGGMTIDGR",
            "cleavage_sites": [25, 50, 80],
            "organism": "Example flavivirus"
        }
    ]
    
    return example_data

def run_demo():
    """Run a complete demonstration of the framework"""
    
    print("🧬 Polyprotein Cleavage Prediction Demo")
    print("=" * 50)
    
    # Create example data
    print("1. Creating example dataset...")
    data = create_example_data()
    
    # Save to file
    with open("demo_data.json", "w") as f:
        json.dump(data, f, indent=2)
    
    print(f"   Created {len(data)} example sequences")
    for item in data:
        seq_len = len(item["sequence"])
        n_sites = len(item["cleavage_sites"])
        print(f"   - {item['protein_id']}: {seq_len} AA, {n_sites} cleavage sites")
    
    # Initialize predictor
    print("\n2. Initializing predictor...")
    predictor = PolyproteinCleavagePredictor(
        config_file="config.json",
        verbose=True
    )
    
    print("   ✓ Predictor initialized")
    
    # Check if we have enough data and computational resources
    print("\n3. Checking system requirements...")
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"   Device: {device}")
    
    if device.type == "cuda":
        gpu_memory = torch.cuda.get_device_properties(0).total_memory
        print(f"   GPU Memory: {gpu_memory / 1e9:.1f} GB")
    
    # For this demo, we'll use CPU-only mode to avoid dependency issues
    if len(data) < 10:
        print("   Note: Demo dataset is very small - not suitable for real training")
        print("   This demo will only show the prediction interface")
        demo_prediction_only()
        return
    
    print("\n4. Training model...")
    print("   (This is a demo with minimal data - real training needs more sequences)")
    
    try:
        # Train the model
        predictor.train(
            data_file="demo_data.json",
            epochs=3,  # Very short training for demo
            batch_size=2,
            save_path="demo_model.pth"
        )
        print("   ✓ Training completed")
        
        # Make predictions
        print("\n5. Making predictions...")
        demo_predictions(predictor, data)
        
    except Exception as e:
        print(f"   Error during training: {e}")
        print("   Falling back to prediction demo...")
        demo_prediction_only()

def demo_prediction_only():
    """Demonstrate the prediction interface without training"""
    
    print("\n4. Demonstrating prediction interface...")
    print("   (Using randomly initialized model - predictions will be meaningless)")
    
    # Initialize predictor
    predictor = PolyproteinCleavagePredictor(verbose=False)
    
    # Example sequence
    test_sequence = "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVA"
    
    print(f"   Test sequence: {test_sequence[:50]}...")
    print(f"   Length: {len(test_sequence)} amino acids")
    
    try:
        # Get cleavage probabilities (this will fail without ESM3, but we'll catch it)
        print("\n5. Computing predictions...")
        print("   Note: This requires ESM3 model download (~2.8GB)")
        print("   Skipping actual prediction in this demo")
        
        # Simulate predictions
        import random
        random.seed(42)
        simulated_probs = [random.random() * 0.3 for _ in range(len(test_sequence))]
        
        # Find high-probability sites
        threshold = 0.15
        predicted_sites = [i for i, prob in enumerate(simulated_probs) if prob > threshold]
        
        print(f"   Threshold: {threshold}")
        print(f"   Predicted cleavage sites: {predicted_sites}")
        
        # Show context around predictions
        for site in predicted_sites[:3]:  # Show first 3
            start = max(0, site - 5)
            end = min(len(test_sequence), site + 6)
            context = test_sequence[start:end]
            pos_in_context = site - start
            print(f"   Site {site}: {context[:pos_in_context]}|{context[pos_in_context+1:]}")
        
    except Exception as e:
        print(f"   Could not compute predictions: {e}")
    
    print("\n6. Demo completed!")
    print("   To run with real data:")
    print("   1. Install dependencies: conda env create -f environment.yml")
    print("   2. Download data: python data_prep.py refseq --email your@email.com")
    print("   3. Train model: python simple_demo.py")

def demo_predictions(predictor, data):
    """Demonstrate predictions on trained model"""
    
    for item in data:
        print(f"\n   Testing: {item['protein_id']}")
        sequence = item["sequence"]
        true_sites = item["cleavage_sites"]
        
        try:
            # Predict cleavage sites
            predicted_sites = predictor.get_cleavage_sites(sequence, threshold=0.5)
            
            print(f"   True sites: {true_sites}")
            print(f"   Predicted: {predicted_sites}")
            
            # Calculate simple metrics
            true_set = set(true_sites)
            pred_set = set(predicted_sites)
            
            if len(true_set) > 0:
                overlap = len(true_set & pred_set)
                precision = overlap / len(pred_set) if len(pred_set) > 0 else 0
                recall = overlap / len(true_set)
                
                print(f"   Precision: {precision:.2f}, Recall: {recall:.2f}")
            
        except Exception as e:
            print(f"   Prediction error: {e}")

def show_requirements():
    """Show system requirements"""
    
    print("\n📋 System Requirements")
    print("-" * 30)
    print("Hardware:")
    print("  - GPU with 4GB+ memory (recommended)")
    print("  - 8GB+ RAM")
    print("  - 10GB+ disk space (for ESM3 model)")
    
    print("\nSoftware:")
    print("  - Python 3.9+")
    print("  - PyTorch with CUDA support")
    print("  - HuggingFace Transformers")
    print("  - ESM package")
    
    print("\nData:")
    print("  - Polyprotein sequences with cleavage annotations")
    print("  - Minimum 100+ sequences for training")
    print("  - Recommended 1000+ sequences")

def main():
    """Main demo function"""
    
    print("Starting polyprotein cleavage prediction demo...")
    
    # Show requirements
    show_requirements()
    
    # Run the demo
    try:
        run_demo()
    except KeyboardInterrupt:
        print("\nDemo interrupted by user")
    except Exception as e:
        print(f"\nDemo failed with error: {e}")
        print("This is expected if dependencies are not fully installed")
        print("See README_cleavage.md for installation instructions")
    
    print("\nDemo finished!")

if __name__ == "__main__":
    main()