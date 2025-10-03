#!/bin/bash

# Exporta as variaveis de ambiente
. /app/env.sh

if [ ! -d "$PIPELINES_DIR" ]; then
    echo "Error: PIPELINES_DIR not defined."
    exit 1
fi

INSTALL_PIPE="$PIPELINES_DIR/predict_occultation/install.sh"

if [ ! -f "$INSTALL_PIPE" ]; then
    echo "Error: Installation script not found."
    exit 1
fi

# Installing pipeline
echo "Installing pipeline..."
. "$INSTALL_PIPE"

set -xe

echo "Environment: Local"


# Run the Python code with the given argument
# predict-run $ARGS || { echo "Failed to predict-run"; exit 1; }
python3 "$PIPE_BASE/local.py" || { echo "Failed to run local.py"; exit 1; }

echo "Done."