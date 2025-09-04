#!/bin/bash

# Check if the argument was given
if [ $# -eq 0 ]; then
    echo "Error: No arguments provided."
    exit 1
fi

ARGS=$@
shift $#

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

echo "Environment: "$EXECUTION_ENV

if [[ "$EXECUTION_ENV" = "linea" ]]
then
    echo "Setup remote env at LIneA (Slurm)"

    export EUPS_USERDATA=/tmp/`whoami`/eups
    . /opt/eups/bin/setups.sh

    echo "Eups Setup: gcc 11.1.0+0"
    setup gcc 11.1.0+0

    echo "Eups Setup: geradata 20240101+0"
    setup geradata 20240101+0

    echo "Eups Setup: elimina 20240101+0"
    setup elimina 20240101+0

    echo "Eups Setup: praia_occ_star_search_12 20240101+0"
    setup praia_occ_star_search_12 20240101+0    

    ulimit -s 100000
    ulimit -u 100000

    echo "Setup environment at LineA done."
fi

# Run the Python code with the given argument
predict-run $ARGS || { echo "Failed to predict-run"; exit 1; }

echo "Done."