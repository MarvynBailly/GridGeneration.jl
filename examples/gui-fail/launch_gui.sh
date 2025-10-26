#!/bin/bash
# GridGeneration GUI Launcher for Unix/Mac
# Run: ./launch_gui.sh

echo "========================================"
echo "  GridGeneration.jl GUI Launcher"
echo "========================================"
echo ""

# Check if Julia is available
if ! command -v julia &> /dev/null
then
    echo "ERROR: Julia not found in PATH"
    echo ""
    echo "Please install Julia or add it to your PATH:"
    echo "https://julialang.org/downloads/"
    echo ""
    exit 1
fi

echo "Julia found!"
echo ""

# Navigate to script directory
cd "$(dirname "$0")"

echo "Starting GUI launcher..."
echo ""

julia launcher.jl

echo ""
echo "GUI closed."
