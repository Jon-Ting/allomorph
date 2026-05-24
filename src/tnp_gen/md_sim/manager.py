import os
import re
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

def check_log_done(log_file: Path, stage: int) -> bool:
    """Checks if a LAMMPS log file indicates completion."""
    if not log_file.exists():
        return False
    
    try:
        with open(log_file, 'rb') as f:
            # Read the last few bytes to check for "DONE!"
            f.seek(0, os.SEEK_END)
            size = f.tell()
            chunk_size = 1024
            if size < chunk_size:
                chunk_size = size
            f.seek(-chunk_size, os.SEEK_END)
            content = f.read().decode('utf-8', errors='ignore')
            
            if stage == 2:
                return "ALL DONE!" in content or "DONE!" in content
            else:
                return "DONE!" in content
    except Exception as e:
        logger.error(f"Error reading log file {log_file}: {e}")
        return False

def update_config(config_file: Path, key: str, value: str = "true"):
    """Updates the config.yml file with a key-value pair."""
    content = ""
    if config_file.exists():
        with open(config_file, 'r') as f:
            content = f.read()
    
    if f"{key}: {value}" in content:
        return

    with open(config_file, 'a') as f:
        f.write(f"{key}: {value}\n")
    logger.info(f"Updated {config_file}: {key}: {value}")

def get_config_value(config_file: Path, key: str) -> str:
    """Gets a value from config.yml."""
    if not config_file.exists():
        return ""
    
    with open(config_file, 'r') as f:
        for line in f:
            if line.startswith(f"{key}:"):
                return line.split(":")[1].strip()
    return ""

def generate_job_list(
    stage: int,
    sim_data_dir: str,
    job_list_file: str,
    queue_list_file: str = "queueList",
):
    """
    Generates a list of jobs to be run based on the current state of simulations.
    Replaces jobList.sh.
    """
    sim_data_path = Path(sim_data_dir)
    job_list_path = Path(job_list_file)
    queue_list_path = sim_data_path / queue_list_file
    
    # Read currently queued jobs
    queued_jobs = set()
    if queue_list_path.exists():
        with open(queue_list_path, 'r') as f:
            queued_jobs = {line.strip() for line in f if line.strip()}

    # Read currently listed jobs to avoid duplicates
    listed_jobs = set()
    if job_list_path.exists():
        with open(job_list_path, 'r') as f:
            listed_jobs = {line.strip() for line in f if line.strip()}

    new_jobs = []

    # Find all Stage X input files
    # Structure is assumed to be SIM_DATA_DIR/Type/NP_Dir/NP_DirSX.in
    for in_file in sim_data_path.glob("*/*/*S%d.in" % stage):
        job_path = in_file.parent / in_file.stem
        job_path_str = str(job_path)
        dir_path = in_file.parent
        config_file = dir_path / "config.yml"
        run_lock = dir_path / "run.lock"
        log_file = dir_path / (in_file.stem + ".log")

        # Check prerequisites
        if stage == 1:
            if get_config_value(config_file, "S0eq") != "true":
                logger.debug(f"{job_path_str} unequilibrated, skipping...")
                continue
        elif stage == 2:
            if get_config_value(config_file, "S1ok") != "true":
                logger.debug(f"{job_path_str} unmelted, skipping...")
                continue

        # Check if already listed or queued
        if job_path_str in listed_jobs:
            logger.debug(f"{job_path_str} already on job list, skipping...")
            continue
        if job_path_str in queued_jobs:
            logger.debug(f"{job_path_str} already in queue, skipping...")
            continue
        
        # Check if running
        if run_lock.exists():
            logger.debug(f"{dir_path} is running a job, skipping...")
            continue

        # Check if already done
        if log_file.exists():
            if check_log_done(log_file, stage):
                key = f"S{stage}ok" if stage > 0 else "S0eq"
                update_config(config_file, key, "true")
                logger.info(f"{job_path_str} already done, skipping...")
                continue
            else:
                logger.info(f"{job_path_str} log exists but not DONE!, resubmitting...")
        else:
            logger.info(f"{job_path_str} ready, adding to job list...")

        new_jobs.append(job_path_str)

    # Append new jobs to the job list file
    if new_jobs:
        with open(job_list_path, 'a') as f:
            for job in new_jobs:
                f.write(job + "\n")
        logger.info(f"Added {len(new_jobs)} jobs to {job_list_file}")
    else:
        logger.info("No new jobs to add.")
