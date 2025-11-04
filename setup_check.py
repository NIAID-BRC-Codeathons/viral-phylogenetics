#!/usr/bin/env python3
"""
Setup verification script

This script checks if all components are properly installed and configured
for polyprotein cleavage site prediction.
"""

import sys
import importlib
import subprocess
import json
import os
from pathlib import Path


def check_python_version():
    """Check Python version"""
    print("🐍 Python Version:")
    version = sys.version_info
    print(f"   Version: {version.major}.{version.minor}.{version.micro}")
    
    if version.major < 3 or (version.major == 3 and version.minor < 8):
        print("   ❌ Python 3.8+ required")
        return False
    else:
        print("   ✓ Python version OK")
        return True


def check_package(package_name, import_name=None, version_attr='__version__'):
    """Check if a package is installed and get version"""
    if import_name is None:
        import_name = package_name
    
    try:
        module = importlib.import_module(import_name)
        
        # Try to get version
        version = "unknown"
        if hasattr(module, version_attr):
            version = getattr(module, version_attr)
        elif hasattr(module, 'version'):
            version = module.version
        elif hasattr(module, '__version__'):
            version = module.__version__
        
        print(f"   ✓ {package_name}: {version}")
        return True, version
        
    except ImportError as e:
        print(f"   ❌ {package_name}: Not installed ({e})")
        return False, None


def check_dependencies():
    """Check all required dependencies"""
    print("\n📦 Dependencies:")
    
    required_packages = [
        ("torch", "torch"),
        ("transformers", "transformers"),
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("scikit-learn", "sklearn"),
        ("biopython", "Bio"),
        ("tqdm", "tqdm"),
        ("requests", "requests"),
    ]
    
    optional_packages = [
        ("esm", "esm"),
        ("fair-esm", "esm"),  # Alternative ESM package
        ("matplotlib", "matplotlib"),
        ("seaborn", "seaborn"),
        ("jupyter", "jupyter"),
    ]
    
    all_good = True
    missing_required = []
    missing_optional = []
    
    print("   Required packages:")
    for package_name, import_name in required_packages:
        success, version = check_package(package_name, import_name)
        if not success:
            all_good = False
            missing_required.append(package_name)
    
    print("   Optional packages:")
    for package_name, import_name in optional_packages:
        success, version = check_package(package_name, import_name)
        if not success:
            missing_optional.append(package_name)
    
    if missing_required:
        print(f"\n   ❌ Missing required packages: {', '.join(missing_required)}")
    
    if missing_optional:
        print(f"   ⚠ Missing optional packages: {', '.join(missing_optional)}")
    
    return all_good, missing_required, missing_optional


def check_pytorch_cuda():
    """Check PyTorch CUDA support"""
    print("\n🚀 PyTorch CUDA:")
    
    try:
        import torch
        
        print(f"   PyTorch version: {torch.__version__}")
        
        cuda_available = torch.cuda.is_available()
        print(f"   CUDA available: {cuda_available}")
        
        if cuda_available:
            device_count = torch.cuda.device_count()
            print(f"   CUDA devices: {device_count}")
            
            for i in range(device_count):
                props = torch.cuda.get_device_properties(i)
                memory_gb = props.total_memory / 1e9
                print(f"   Device {i}: {props.name} ({memory_gb:.1f} GB)")
            
            return True
        else:
            print("   ⚠ No CUDA devices available (will use CPU)")
            return False
            
    except ImportError:
        print("   ❌ PyTorch not installed")
        return False


def check_esm_models():
    """Check ESM model availability"""
    print("\n🧬 ESM Models:")
    
    try:
        import esm
        print("   ✓ ESM package installed")
        
        # Try to load a model (this will download if not cached)
        print("   Testing model loading...")
        try:
            model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
            print("   ✓ ESM2 model loaded successfully")
            del model  # Free memory
            return True
        except Exception as e:
            print(f"   ❌ ESM model loading failed: {e}")
            return False
            
    except ImportError:
        print("   ❌ ESM package not installed")
        return False


def check_project_files():
    """Check if project files exist"""
    print("\n📁 Project Files:")
    
    required_files = [
        "esmretrain.py",
        "data_prep.py",
        "environment.yml",
        "config.json"
    ]
    
    optional_files = [
        "environment-minimal.yml",
        "simple_demo.py",
        "README_cleavage.md",
        "README_refseq.md",
        "INSTALL.md",
        "QUICKSTART.md"
    ]
    
    all_present = True
    
    print("   Required files:")
    for filename in required_files:
        if Path(filename).exists():
            print(f"   ✓ {filename}")
        else:
            print(f"   ❌ {filename}")
            all_present = False
    
    print("   Optional files:")
    for filename in optional_files:
        if Path(filename).exists():
            print(f"   ✓ {filename}")
        else:
            print(f"   ⚠ {filename}")
    
    return all_present


def check_config_file():
    """Check config.json validity"""
    print("\n⚙️ Configuration:")
    
    try:
        with open("config.json", 'r') as f:
            config = json.load(f)
        
        required_sections = ['model', 'training', 'data', 'evaluation']
        
        for section in required_sections:
            if section in config:
                print(f"   ✓ {section} section present")
            else:
                print(f"   ❌ {section} section missing")
                return False
        
        # Check some key parameters
        batch_size = config.get('training', {}).get('batch_size', 32)
        print(f"   Batch size: {batch_size}")
        
        use_cuda = config.get('device', {}).get('use_cuda', True)
        print(f"   CUDA enabled: {use_cuda}")
        
        return True
        
    except FileNotFoundError:
        print("   ❌ config.json not found")
        return False
    except json.JSONDecodeError as e:
        print(f"   ❌ Invalid JSON in config.json: {e}")
        return False


def test_basic_functionality():
    """Test basic framework functionality"""
    print("\n🧪 Basic Functionality:")
    
    try:
        # Test imports
        print("   Testing imports...")
        from esmretrain import PolyproteinCleavagePredictor
        print("   ✓ Core framework imports OK")
        
        from data_prep import validate_data_format
        print("   ✓ Data preparation imports OK")
        
        # Test predictor initialization
        print("   Testing predictor initialization...")
        predictor = PolyproteinCleavagePredictor(verbose=False)
        print("   ✓ Predictor initialization OK")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Functionality test failed: {e}")
        return False


def test_data_validation():
    """Test data validation with sample data"""
    print("\n🔍 Data Validation:")
    
    try:
        # Create test data
        test_data = [
            {
                "protein_id": "test_seq",
                "sequence": "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLS",
                "cleavage_sites": [10, 20, 30],
                "organism": "Test virus"
            }
        ]
        
        # Save to temp file
        test_file = "test_validation.json"
        with open(test_file, 'w') as f:
            json.dump(test_data, f)
        
        # Test validation
        from data_prep import validate_data_format
        is_valid = validate_data_format(test_file)
        
        # Clean up
        os.unlink(test_file)
        
        if is_valid:
            print("   ✓ Data validation working")
            return True
        else:
            print("   ❌ Data validation failed")
            return False
            
    except Exception as e:
        print(f"   ❌ Data validation test error: {e}")
        return False


def check_disk_space():
    """Check available disk space"""
    print("\n💾 Disk Space:")
    
    try:
        import shutil
        
        # Check current directory
        total, used, free = shutil.disk_usage(".")
        free_gb = free / (1024**3)
        
        print(f"   Available space: {free_gb:.1f} GB")
        
        if free_gb < 5:
            print("   ⚠ Low disk space (recommend 10GB+ for models)")
            return False
        elif free_gb < 10:
            print("   ⚠ Limited disk space (10GB+ recommended)")
            return True
        else:
            print("   ✓ Sufficient disk space")
            return True
            
    except Exception as e:
        print(f"   ❌ Could not check disk space: {e}")
        return False


def generate_report(results):
    """Generate a setup report"""
    print("\n" + "=" * 60)
    print("📋 SETUP REPORT")
    print("=" * 60)
    
    total_checks = len(results)
    passed_checks = sum(1 for result in results.values() if result)
    
    print(f"Overall Status: {passed_checks}/{total_checks} checks passed")
    print()
    
    for check_name, status in results.items():
        icon = "✓" if status else "❌"
        print(f"{icon} {check_name}")
    
    print()
    
    if passed_checks == total_checks:
        print("🎉 All checks passed! System is ready for training.")
        
        print("\nNext steps:")
        print("1. Download data: python data_prep.py refseq --email your@email.com")
        print("2. Run demo: python simple_demo.py")
        print("3. Start training!")
        
    elif passed_checks >= total_checks * 0.8:
        print("⚠ Most checks passed. Address any ❌ issues before training.")
        
        print("\nRecommended actions:")
        print("1. Install missing packages")
        print("2. Check CUDA setup if needed")
        print("3. Verify file integrity")
        
    else:
        print("❌ Multiple issues detected. Setup needs attention.")
        
        print("\nRequired actions:")
        print("1. Check installation guide: INSTALL.md")
        print("2. Install missing dependencies")
        print("3. Re-run this check")


def main():
    """Main setup check function"""
    
    print("🔧 Polyprotein Cleavage Prediction Setup Check")
    print("This script verifies your installation and configuration")
    print("=" * 60)
    
    # Run all checks
    results = {}
    
    results["Python Version"] = check_python_version()
    
    deps_ok, missing_req, missing_opt = check_dependencies()
    results["Dependencies"] = deps_ok
    
    results["PyTorch CUDA"] = check_pytorch_cuda()
    results["ESM Models"] = check_esm_models()
    results["Project Files"] = check_project_files()
    results["Configuration"] = check_config_file()
    results["Basic Functionality"] = test_basic_functionality()
    results["Data Validation"] = test_data_validation()
    results["Disk Space"] = check_disk_space()
    
    # Generate report
    generate_report(results)
    
    # Return success status
    return all(results.values())


if __name__ == "__main__":
    try:
        success = main()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\nSetup check interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nSetup check failed: {e}")
        sys.exit(1)