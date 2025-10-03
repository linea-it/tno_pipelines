#!/usr/bin/env python3

import os
import pathlib
import subprocess
import sys
import time

import yaml
from config import Config
from packages.dao.task import PredictionState, TaskDao
from packages.dao.worker import WorkerDao
from sqlalchemy.orm import Session



def run_command_and_stream_output(command):
    """
    Executes a shell command and streams its output line by line.
    """
    try:
        # Popen starts the process and allows interaction with its I/O
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,  # Decodes stdout/stderr as text
            shell=True,  # Executes the command through the shell
        )

        # Stream stdout line by line
        for line in process.stdout:
            print(line.strip())
            sys.stdout.flush()  # Ensure immediate output

        # After stdout is exhausted, check for any stderr output
        stderr_output = process.stderr.read()
        if stderr_output:
            print(stderr_output.strip())

        # Wait for the process to terminate and get the return code
        process.wait()
        if process.returncode != 0:
            print(f"Command '{command}' exited with error code {process.returncode}")

    except FileNotFoundError:
        print(f"Error: Command '{command}' not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def main():
    # Check for tasks with status QUEUED
    # print("Checking for QUEUED tasks...")
    task_dao = TaskDao()
    db_session = Session(task_dao.get_db_engine())
    command = None
    try:
        task = task_dao.get_next_task(db_session, PredictionState.QUEUED)
        if not task:
            # print("No QUEUED tasks found.")
            time.sleep(interval)
            return
        print("=" * 60)
        print(f"Found task id: [{task.id}]")  # type: ignore

        if not task.workdir:
            raise ValueError("Task workdir is not set.")

        cfg = Config()
        cfg.task_id = task.id  # type: ignore
        cfg.asteroid_name = task.asteroid_id  # type: ignore
        cfg.asteroid_path = task.workdir  # type: ignore

        cfg_filepath = pathlib.Path(task.workdir).joinpath("config.yaml")  # type: ignore

        with open(cfg_filepath, "w") as outfile:
            data_json = cfg.model_dump()
            print(f"Writing config to {cfg_filepath}")
            print(data_json)

            yaml.dump(data_json, outfile)

        # Execute the run.sh script
        print("Executing run.sh")
        ret = (
            pathlib.Path(os.getenv("PIPELINES_DIR"))
            .joinpath("predict_occultation")
            .joinpath("run.sh")
            .resolve()
            .as_posix()
        )
        command = f"{ret} {str(cfg_filepath)}"
        print(f"Running command: {command}")
        print("-" * 60)

    except Exception as e:
        print(f"Error: {e}")
        db_session.rollback()
    finally:
        db_session.close()

    if command:
        run_command_and_stream_output(command)


if __name__ == "__main__":

    interval = int(os.getenv("INTERVAL", 2))

    print("=" * 60)
    print("Local Runner")
    worker_dao = WorkerDao()
    worker_name = os.getenv("WORKER_NAME", "worker_runner_1")
    worker_dao.initialize_heartbeat(worker_name)

    while True:
        try:
            main()
        except Exception as e:
            print(f"Unexpected error in main loop: {e}")
        finally:
            worker_dao.send_heartbeat(worker_name)
