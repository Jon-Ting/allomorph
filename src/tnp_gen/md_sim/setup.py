import os
import shutil
import glob
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

def setup_md_sim(init_struct_dir: str, target_dir: str):
    """
    Sets up MD simulation directories by copying initial structure files.

    Args:
        init_struct_dir: Path to the directory containing initial structures (e.g., InitStruct/TNP).
        target_dir: Path where simulation directories should be created.
    """
    init_struct_path = Path(init_struct_dir)
    target_path = Path(target_dir)

    if not target_path.exists():
        target_path.mkdir(parents=True)
        logger.info(f"Created target directory: {target_path}")

    # Find all directories in InitStruct/TNP
    tnp_dirs = [d for d in init_struct_path.iterdir() if d.is_dir()]
    
    for tnp_dir in tnp_dirs:
        dir_name = tnp_dir.name
        tnp_target_dir = target_path / dir_name
        
        if not tnp_target_dir.exists():
            tnp_target_dir.mkdir(parents=True)
            logger.info(f"Created TNP target directory: {tnp_target_dir}")

        # Find all .lmp files
        lmp_files = list(tnp_dir.glob("*.lmp"))
        
        for lmp_file in lmp_files:
            file_name = lmp_file.name
            file_stem = lmp_file.stem
            
            sim_dir = tnp_target_dir / file_stem
            if not sim_dir.exists():
                sim_dir.mkdir(parents=True)
                logger.info(f"Created simulation directory: {sim_dir}")
            
            dest_file = sim_dir / file_name
            if not dest_file.exists():
                shutil.copy(lmp_file, dest_file)
                # Create an empty config.yml as in the original script
                (sim_dir / "config.yml").touch()
                logger.info(f"Copied {file_name} to {sim_dir}")
            else:
                logger.info(f"{file_name} already exists in {sim_dir}, skipped!")
