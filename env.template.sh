#!/bin/bash

export PIPELINES_DIR="<PIPELINES_DIR>"
export PYTHONPATH=$PIPELINES_DIR/predict_occultation/packages
export PATH=$PATH:$PIPELINES_DIR/predict_occultation/scripts

# Execution Environment local or linea (linea use eups for setup fortrans)
export EXECUTION_ENV="local"

# Database Connections URI
export DB_ADMIN_URI="postgresql+psycopg2://postgres:postgres@localhost:5432/postgres"
export DB_CATALOG_URI="postgresql+psycopg2://postgres:postgres@localhost:5432/prod_gavo"

# External Inputs Path

# de440.bsp complete filepath
PLANETARY_EPHEMERIS_FILEPATH="/app/de440.bsp"
# naif0012.tls complete filepath
LEAP_SECONDS_FILEPATH="/app/naif0012.tls"