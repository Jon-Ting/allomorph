import os
import subprocess
import logging
import re
from pathlib import Path

logger = logging.getLogger(__name__)

def get_num_atoms(lmp_file: Path) -> int:
    """Reads the number of atoms from a LAMMPS data file."""
    try:
        with open(lmp_file, 'r') as f:
            for line in f:
                if 'atoms' in line:
                    return int(line.split()[0])
    except Exception as e:
        logger.error(f"Error reading {lmp_file}: {e}")
    return 0

def calculate_resources(num_atoms: int, stage: int):
    """Calculates ncpus, memory, and walltime based on number of atoms and stage."""
    ncpus = ((num_atoms - 1) // 64000 + 1) * 4
    mem = (num_atoms // 360000 + 1) * ncpus // 2  # GB
    wall_time = (36 * num_atoms) // ncpus  # s
    
    if stage == 1:
        mem = int(mem * 0.6) + 1
        wall_time = wall_time * 3
    elif stage == 2:
        mem = int(mem * 0.8) + 1
        wall_time = wall_time * 4
        
    if wall_time > 172800:
        wall_time = 172800
        
    hours = wall_time // 3600
    minutes = (wall_time % 3600) // 60
    seconds = wall_time % 60
    
    wall_time_str = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
    
    return ncpus, mem, wall_time_str

def generate_pbs_script(
    job_name: str,
    job_dir: Path,
    ncpus: int,
    wall_time: str,
    mem: int,
    project: str,
    email: str = None,
):
    """Generates a PBS script for a single LAMMPS job."""
    pbs_content = f"""#!/bin/bash
#PBS -P {project}
#PBS -q normal
#PBS -l ncpus={ncpus},walltime={wall_time},mem={mem}GB,jobfs=2GB
#PBS -l storage=scratch/{project}+gdata/{project}
#PBS -l wd
#PBS -N {job_name}
"""
    if email:
        pbs_content += f"#PBS -M {email}\n#PBS -m a\n"
    
    pbs_content += f"""
module load lammps/29Sep2021

cd {job_dir}
mpirun -np {ncpus} lmp_openmpi -in {job_name}.in > {job_name}.log
"""
    return pbs_content

def submit_jobs(
    job_list_file: str,
    sim_data_dir: str,
    project: str,
    max_queue_num: int = 15,
    email: str = None,
):
    """
    Submits jobs from a list to the PBS queue.
    Replaces the logic in subAnneal.sh.
    """
    job_list_path = Path(job_list_file)
    sim_data_path = Path(sim_data_dir)
    
    if not job_list_path.exists():
        logger.error(f"Job list {job_list_file} not found!")
        return

    # Check current queue
    try:
        qselect_out = subprocess.check_output(["qselect", "-u", os.environ.get("USER")], text=True)
        num_in_queue = len(qselect_out.splitlines())
    except Exception:
        num_in_queue = 0
        logger.warning("Could not check queue status, assuming 0.")

    num_to_sub = max(0, max_queue_num - num_in_queue)
    logger.info(f"maxQueueNum: {max_queue_num}, numInQueue: {num_in_queue}, numToSub: {num_to_sub}")

    with open(job_list_path, 'r') as f:
        jobs = f.read().splitlines()

    if not jobs:
        logger.info("No more jobs in job list.")
        return

    jobs_to_run = jobs[:num_to_sub]
    remaining_jobs = jobs[num_to_sub:]

    for job_path_str in jobs_to_run:
        # job_path_str might look like: /scratch/q27/jt5911/CL10S/AuPdPt30_.../AuPdPt30_...S0
        # We need to extract the stage and the directory
        
        full_job_path = Path(job_path_str)
        job_dir = full_job_path.parent
        job_name = full_job_path.name # e.g., AuPdPt30_...S0
        
        # Extract stage from job name
        match = re.search(r'S(\d+)$', job_name)
        if not match:
            logger.warning(f"Could not determine stage from {job_name}, skipping...")
            continue
        stage = int(match.group(1))
        
        # Get number of atoms from the .lmp file (which should have been copied in setup)
        # The .lmp file name is the part of job_name before S{stage}
        lmp_file_name = job_name[:match.start()] + ".lmp"
        lmp_file = job_dir / lmp_file_name
        
        if not lmp_file.exists():
            logger.warning(f"LMP file {lmp_file} not found, skipping...")
            continue
            
        num_atoms = get_num_atoms(lmp_file)
        if num_atoms == 0:
            logger.warning(f"Could not get number of atoms for {job_name}, skipping...")
            continue
            
        ncpus, mem, wall_time = calculate_resources(num_atoms, stage)
        
        pbs_script = generate_pbs_script(job_name, job_dir, ncpus, wall_time, mem, project, email)
        
        pbs_file = job_dir / f"{job_name}.pbs"
        with open(pbs_file, 'w') as f:
            f.write(pbs_script)
            
        logger.info(f"Generated PBS script for {job_name}: ncpus={ncpus}, mem={mem}GB, walltime={wall_time}")
        
        # Submit the job
        try:
            subprocess.run(["qsub", str(pbs_file)], check=True)
            logger.info(f"Submitted {job_name}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to submit {job_name}: {e}")

    # Update job list
    with open(job_list_path, 'w') as f:
        f.write("\n".join(remaining_jobs) + ("\n" if remaining_jobs else ""))
