"""Setup and run NCPac feature extraction."""

import logging
import os
import re
import shutil
import subprocess
from multiprocessing import Pool
from pathlib import Path
from zipfile import ZipFile

from np_gen.constants import TNP_DISTRIB_LIST

logger = logging.getLogger(__name__)

# Default configuration (override via function arguments)
DEFAULT_NUM_FRAME_PER_NP = 11
DEFAULT_ZFILL_NUM = 5
DEFAULT_DONE_FILE = "DONE.txt"
DEFAULT_OUT_MD_FILE = "MDout.csv"
DEFAULT_NCPAC_EXE_NAME = "NCPac.exe"
DEFAULT_NCPAC_INP_NAME = "NCPac.inp"
DEFAULT_HEADER_LINE = "CSIRO Nanostructure Databank - {} Nanoparticle Data Set"


def setup_ncpac(
    sim_data_dir: str,
    target_dir: str,
    path2_ncpac_exe: str,
    path2_ncpac_inp: str,
    ele_comb: str = None,
    source_dirs: list = None,
    num_frame_per_np=DEFAULT_NUM_FRAME_PER_NP,
    zfill_num=DEFAULT_ZFILL_NUM,
    done_file=DEFAULT_DONE_FILE,
    out_md_file=DEFAULT_OUT_MD_FILE,
    ncpac_exe_name=DEFAULT_NCPAC_EXE_NAME,
    ncpac_inp_name=DEFAULT_NCPAC_INP_NAME,
    header_line=None,
):
    """Copy xyz files to individual directories and relabel numerically.

    Args:
        sim_data_dir: Path to the simulation data directory.
        target_dir: Path where NCPac directories will be created.
        path2_ncpac_exe: Path to the NCPac executable.
        path2_ncpac_inp: Path to the NCPac input template.
        ele_comb: Element combination string (e.g. "AuPtPd"). If None, it is
            inferred from the subdirectories of *sim_data_dir*.
        source_dirs: List of source subdirectories to process. If None, all
            subdirectories of *sim_data_dir* are used.
    """
    if source_dirs is None:
        sim_data_path = Path(sim_data_dir)
        if sim_data_path.exists():
            source_dirs = sorted([d.name for d in sim_data_path.iterdir() if d.is_dir()])
        else:
            source_dirs = TNP_DISTRIB_LIST
    if ele_comb is None:
        # Best-effort inference: use the parent directory name if it looks like
        # an element combination, otherwise leave it for the caller to specify.
        sim_data_path = Path(sim_data_dir)
        ele_comb = sim_data_path.name if sim_data_path.name else "Nanoparticle"
    header_line = header_line or DEFAULT_HEADER_LINE.format(ele_comb)

    sim_data_path = Path(sim_data_dir)
    target_path = Path(target_dir)
    path2_ncpac_exe = Path(path2_ncpac_exe)
    path2_ncpac_inp = Path(path2_ncpac_inp)

    source_paths = [sim_data_path / d for d in source_dirs]

    logger.info("Copying xyz files to individual directories and relabelling numerically...")
    np_cnt = 0
    working_list = []

    target_path.mkdir(parents=True, exist_ok=True)

    out_md_path = target_path / out_md_file
    if not out_md_path.exists():
        with open(out_md_path, "w") as f:
            f.write("confID,T,P,PE,KE,TE\n")

    for source_path in source_paths:
        if not source_path.exists():
            logger.warning(f"Source path {source_path} does not exist, skipping...")
            continue

        for np_dir_name in os.listdir(source_path):
            np_dir_path = source_path / np_dir_name
            if not np_dir_path.is_dir():
                continue

            logger.info(f"  Nanoparticle: {np_dir_name}")

            # If not done for the nanoparticle yet, reextract Stage 2 files
            if not (np_dir_path / done_file).exists():
                zip_file = np_dir_path / f"{np_dir_name}S2.zip"
                if zip_file.exists():
                    with ZipFile(zip_file, "r") as zf:
                        zf.extractall(np_dir_path)
                    logger.info("    Extracted Stage 2 files...")
                else:
                    logger.warning(f"    Zip file {zip_file} not found, skipping...")
                    continue
            else:
                np_cnt += 1
                continue

            conf_cnt = np_cnt * num_frame_per_np
            s2_subdir = np_dir_path / f"{np_dir_name}S2"
            if not s2_subdir.exists():
                logger.warning(f"    Subdirectory {s2_subdir} not found, skipping...")
                continue

            all_s2_nps = [
                np for np in os.listdir(s2_subdir)
                if "min" in np
            ]

            for np_conf in sorted(
                all_s2_nps,
                key=lambda key: [int(i) for i in re.findall(r"min\.([0-9]+)", key)],
            ):
                ori_file_path = s2_subdir / np_conf
                conf_id = str(conf_cnt).zfill(zfill_num)
                logger.info(f"  Conformation ID: {conf_id}")
                conf_dir = target_path / conf_id
                conf_dir.mkdir(parents=True, exist_ok=True)

                logger.info("    Copying files...")
                with open(ori_file_path, "r") as f1:
                    with open(conf_dir / f"{conf_id}.xyz", "w") as f2:
                        f2.write(f1.readline())
                        f1.readline()  # Replace second line with CSIRO header
                        f2.write(f"{header_line} - {np_dir_name}\n")
                        f2.write("".join([line for line in f1.readlines()]))

                logger.info("    Replaced header...")
                shutil.copy(path2_ncpac_exe, conf_dir / ncpac_exe_name)
                shutil.copy(path2_ncpac_inp, conf_dir / ncpac_inp_name)

                # Update NCPac.inp with the correct xyz file name
                with open(conf_dir / ncpac_inp_name, "r") as f:
                    inp_lines = f.readlines()

                # Find the line that needs modification (usually the first one)
                if inp_lines:
                    inp_lines[0] = f"{conf_id}.xyz       - name of xyz input file                                              [in_filexyz]\n"

                with open(conf_dir / ncpac_inp_name, "w") as f:
                    f.writelines(inp_lines)

                working_list.append((str(conf_dir), conf_id))
                conf_cnt += 1

            # Extract outputs from MD for each configuration
            log_file = np_dir_path / f"{np_dir_name}S2.log"
            if log_file.exists():
                with open(log_file, "r") as f1:
                    with open(out_md_path, "a") as f2:
                        found_min_line, prev_line, current_conf_cnt = False, None, conf_cnt - num_frame_per_np
                        for line in f1:
                            if "- MINIMISATION -" in line and not found_min_line:
                                found_min_line = True
                            elif "Loop time of" in line and found_min_line:
                                if prev_line:
                                    conf_id = str(current_conf_cnt).zfill(zfill_num)
                                    parts = prev_line.split()
                                    if len(parts) >= 7:
                                        temp, pres, pot_e, kin_e, tot_e = parts[2:7]
                                        f2.write(f"{conf_id},{temp},{pres},{pot_e},{kin_e},{tot_e}\n")
                                        current_conf_cnt += 1
                                found_min_line = False
                            prev_line = line

            np_cnt += 1

            # Clean up directory and mark as done
            shutil.rmtree(s2_subdir)
            (np_dir_path / done_file).touch()
            logger.info(f"   {np_dir_name} Done!")

    return working_list


def run_ncpac(work, final_data_path: str, ncpac_exe_name=DEFAULT_NCPAC_EXE_NAME, verbose=False):
    """Execute NCPac.exe for a single work item."""
    conf_dir, conf_id = work
    final_data_path = Path(final_data_path)

    if verbose:
        logger.info(f"    Running NCPac for {conf_id}...")

    try:
        # Run NCPac using subprocess
        result = subprocess.run(
            [f"./{ncpac_exe_name}"],
            cwd=conf_dir,
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            logger.error(f"NCPac failed for {conf_id}: {result.stderr}")

        feature_file = Path(conf_dir) / "od_FEATURESET.csv"
        if not feature_file.exists():
            # If execution unsuccessful
            logger.warning(f"Feature file not found for {conf_id}, creating empty file.")
            (final_data_path / "Features" / f"{conf_id}.csv").touch()
        else:
            (final_data_path / "Structures").mkdir(parents=True, exist_ok=True)
            (final_data_path / "Features").mkdir(parents=True, exist_ok=True)
            shutil.copy(Path(conf_dir) / f"{conf_id}.xyz", final_data_path / "Structures")
            shutil.copy(feature_file, final_data_path / "Features" / f"{conf_id}.csv")

        # Remove unnecessary files
        conf_path = Path(conf_dir)
        for pattern in ["*.mod", "fort.*", "ov_*"]:
            for f in conf_path.glob(pattern):
                f.unlink()
        for f in conf_path.glob("od_*"):
            if f.name != "od_FEATURESET.csv":
                f.unlink()

        (conf_path / DEFAULT_DONE_FILE).touch()
        logger.info(f"   {conf_id} Done!")

    except Exception as e:
        logger.error(f"Error running NCPac for {conf_id}: {e}")


def run_ncpac_parallel(remaining_work, final_data_path, ncpac_exe_name=None):
    """Run NCPac in parallel over remaining work items."""
    final_data_path = Path(final_data_path)
    final_data_path.mkdir(parents=True, exist_ok=True)
    (final_data_path / "Structures").mkdir(parents=True, exist_ok=True)
    (final_data_path / "Features").mkdir(parents=True, exist_ok=True)

    with Pool() as p:
        p.starmap(
            run_ncpac,
            [(w, str(final_data_path), ncpac_exe_name) for w in remaining_work],
        )
