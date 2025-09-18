#!/bin/bash

set -euxo pipefail

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Define paths relative to the script location
CODON_PROGRAM="${SCRIPT_DIR}/code/main.py"        # Python file in same directory as script
PYTHON_PROGRAM="${SCRIPT_DIR}/genome-assembly-copy/main.py"  # Codon file in same directory as script
DATA_BASE_DIR="${SCRIPT_DIR}/data"                 # Data folder is sibling to code folder

# Check if Python program exists
if [ ! -f "$PYTHON_PROGRAM" ]; then
    echo "Error: Python program not found at $PYTHON_PROGRAM"
    exit 1
fi

# Check if data base directory exists
if [ ! -d "$DATA_BASE_DIR" ]; then
    echo "Error: Data directory not found at $DATA_BASE_DIR"
    exit 1
fi

echo "Dataset   Language    Runtime   N50"
echo "----------------------------------------"

# Process each data folder (data1 to data4)
for i in {1..4}; do

    if [ "$i" -eq 4 ]; then
      ulimit -s 8192000
    fi

    data_folder="${DATA_BASE_DIR}/data${i}"

    # Check if data folder exists
    if [ ! -d "$data_folder" ]; then
        echo "Warning: Data folder $data_folder not found, skipping..."
        continue
    fi

    c_output=$(codon run -release "$CODON_PROGRAM" "$data_folder")
    echo "data${i}    codon   $c_output"

    p_output=$(python "$PYTHON_PROGRAM" "$data_folder")
    echo "data${i}    python   $p_output"
done