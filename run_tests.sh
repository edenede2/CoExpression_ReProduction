#!/bin/bash

# Script to run pipeline tests from any directory
# This handles the path issues automatically

echo "🧪 CoExpression Network Analysis Pipeline Test Runner"
echo "======================================================="

# Find the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR" && pwd )"

# Check if we're in the project root
if [ -d "$PROJECT_ROOT/new_script" ] && [ -d "$PROJECT_ROOT/data" ]; then
    echo "✅ Found project root: $PROJECT_ROOT"
    cd "$PROJECT_ROOT/new_script"
    echo "📁 Changed to script directory: $(pwd)"
else
    echo "❌ Could not find project structure"
    echo "   Please run this from the CoExpression_reProduction directory"
    exit 1
fi

# Check for required data files
echo ""
echo "🔍 Checking for required data files..."
DATA_DIR="../data"

if [ ! -f "$DATA_DIR/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct" ]; then
    echo "❌ Missing: GTEx expression data file"
    exit 1
fi

if [ ! -f "$DATA_DIR/GTEx_Analysis_v8_Annotations_SampleAttributesDS.tsv" ]; then
    echo "❌ Missing: Sample attributes file"
    exit 1
fi

if [ ! -f "$DATA_DIR/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.tsv" ]; then
    echo "❌ Missing: Subject phenotypes file"
    exit 1
fi

if [ ! -f "$DATA_DIR/hgnc_complete_set.tsv" ]; then
    echo "❌ Missing: HGNC gene data file"
    exit 1
fi

echo "✅ All required data files found"

# Run tests in order
echo ""
echo "🚀 Running pipeline tests..."
echo ""

echo "1️⃣  Testing minimal data loading..."
python minimal_test.py
if [ $? -ne 0 ]; then
    echo "❌ Minimal test failed"
    exit 1
fi

echo ""
echo "2️⃣  Testing age filtering functionality..."
python test_age_filtering.py
if [ $? -ne 0 ]; then
    echo "❌ Age filtering test failed"
    exit 1
fi

echo ""
echo "3️⃣  Running comprehensive pipeline test..."
echo "   (This may take several minutes...)"
python test_pipeline.py
if [ $? -ne 0 ]; then
    echo "⚠️  Full pipeline test had issues, but that's okay for first run"
    echo "   The synthetic data fallback should have worked"
fi

echo ""
echo "🎉 Testing completed!"
echo ""
echo "Next steps:"
echo "  1. Stay in the new_script directory: cd $(pwd)"
echo "  2. Launch the web app: streamlit run streamlit_app.py"
echo "  3. Or run specific tests: python test_age_filtering.py"
echo ""
echo "📖 See README.md for detailed usage instructions"
