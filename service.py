# import os
# original_path = os.getcwd()
# os.environ["EXECUTION_PATH"] = original_path

from pathlib import Path
import uuid
from library import get_configs, write_job_file
from dao import OrbitTraceJobDao, OrbitTraceJobResultDao, PredictOccultationJobDao
# from orbit_trace import orbit_trace_job_to_run, orbit_trace_has_job_running, orbit_trace_make_job_json_file
from datetime import timezone
import shutil
import subprocess
from datetime import datetime, timezone
from predict_occultation import predict_job_queue, run_job as run_predict_job, rerun_job, check_tasks, ingest_job_results
# from orbit_trace import run_job, ingest_job_results, main

# from external_inputs.jpl import get_bsp_from_jpl
# get_bsp_from_jpl("Chiron", "2023-01-01", "2023-12-31", "/lustre/t1/tmp/tno/tmp", "Chiron.bsp")

run_predict_job(74)

# main('/lustre/t1/tmp/tno/orbit_trace/17-cdcb626a')
# ingest_job_results('/lustre/t1/tmp/tno/orbit_trace/16-6930783e', 16)

# Como iniciar o Celery
# celery -A tno_celery worker -l INFO
# celery -A tno_celery beat -l INFO

# celery -A tno_celery worker -Q single -c 1 --detach -l INFO --pidfile="/lustre/t1/tmp/tno/pipelines/tmp/%n.pid" --logfile="/lustre/t1/tmp/tno/pipelines/tmp/%n%I.log"
# celery -A tno_celery worker -Q default --detach -l INFO --pidfile="/lustre/t1/tmp/tno/pipelines/tmp/%n.pid" --logfile="/lustre/t1/tmp/tno/pipelines/tmp/%n%I.log"
# celery -A tno_celery beat --detach -l INFO --pidfile="/lustre/t1/tmp/tno/pipelines/tmp/celeryd.pid" --logfile="/lustre/t1/tmp/tno/pipelines/tmp/celeryd.log"



# Listar todos processos do meu usuario 
# ps -f -U 15161

# Comando para listar os processos do celery worker
# ps aux|grep 'celery worker'
# ps aux|grep 'celery beat'

# Comando para matar todos os processos do celery worker
# ps auxww | grep 'celery worker' | awk '{print $2}' | xargs kill -9
# ps auxww | grep 'celery beat' | awk '{print $2}' | xargs kill -9

# run_predict_job(57)
# rerun_job(57)
# check_tasks(56)
# predict_job_queue()

# from pathlib import Path

# rerun_job(64)
# from time import sleep
# flag = False
# while flag != True:
#     flag = check_tasks(64)
#     sleep(30)

# fp = Path('/lustre/t1/tmp/tno/predict_occultation/57')
# count_results_ingested = ingest_job_results(fp, 57)

# check_tasks(61)

# from predict_occultation import get_job_running

# running = get_job_running()
# print(running)
