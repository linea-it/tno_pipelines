# Pipeline template
Pipeline template for the Orchestration.

## Install

The only requirement is to have miniconda or anaconda previously installed:

```bash
git clone https://github.com/linea-it/pipeline_template && cd pipeline_template
cp env.template.sh env.sh
```

Edit `env.sh` setting the path to this repo in `$PIPELINES_DIR` and execute:

```bash
source env.sh
```

## Run a pipeline

Create config:

```bash
pip install pydantic
cd $PIPELINES_DIR/predict_occultation
python config.py
```

To execute, simply:

```bash
# execute Hello World
cd $PIPELINES_DIR/predict_occultation
mkdir process001
# ./run.sh config.yaml process001
./run.sh config.yaml /app/predict_occultation/process001
```
