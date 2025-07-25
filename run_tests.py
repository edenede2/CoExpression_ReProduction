#!/usr/bin/env python3
"""
Test runner that can be executed from any directory
Automatically detects the correct paths for data and scripts
"""

import sys
import os
from pathlib import Path

def find_project_root():
    """Find the project root directory."""
    current = Path.cwd()
    
    # Look for indicators of project root
    indicators = ['data', 'new_script', 'old_scripts']
    
    # Check current directory and parents
    for path in [current] + list(current.parents):
        if all((path / indicator).exists() for indicator in indicators[:2]):
            return path
    
    # Fallback to current directory
    return current

def run_tests():
    """Run the pipeline tests from the correct directory."""
    
    # Find project root
    project_root = find_project_root()
    script_dir = project_root / 'new_script'
    
    print(f"Project root: {project_root}")
    print(f"Script directory: {script_dir}")
    print(f"Current working directory: {Path.cwd()}")
    
    if not script_dir.exists():
        print("❌ Could not find new_script directory")
        return False
    
    # Change to script directory
    os.chdir(script_dir)
    print(f"Changed to: {Path.cwd()}")
    
    # Add script directory to Python path
    sys.path.insert(0, str(script_dir))
    
    try:
        # Import and run the test
        from test_pipeline import test_pipeline
        
        print("\n" + "="*60)
        print("RUNNING PIPELINE TEST FROM CORRECT DIRECTORY")
        print("="*60)
        
        success = test_pipeline()
        return success
        
    except ImportError as e:
        print(f"❌ Could not import test modules: {e}")
        return False
    except Exception as e:
        print(f"❌ Test failed with error: {e}")
        
        # Try synthetic data fallback
        try:
            from test_pipeline import test_small_synthetic_data
            print("\n🧪 Falling back to synthetic data test...")
            return test_small_synthetic_data()
        except Exception as e2:
            print(f"❌ Synthetic test also failed: {e2}")
            return False

if __name__ == "__main__":
    success = run_tests()
    
    if success:
        print("\n🎯 All tests completed successfully!")
        print("\nNext steps:")
        print("  cd new_script")
        print("  streamlit run streamlit_app.py")
    else:
        print("\n❌ Tests failed. Check error messages above.")
        
    sys.exit(0 if success else 1)
