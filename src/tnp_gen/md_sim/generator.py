import os
import re
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

def generate_lammps_input(
    stage: int,
    sim_data_dir: str,
    init_struct_base_dir: str,
    eam_dir: str,
    template_dir: str,
    tnp_types: list = None,
    init_temp: int = 300,
    heat_temp: int = 2100,
    heat_rate: int = 10, # K/ps
    total_dumps: int = 10,
    s0_period: int = 20000,
    s0_ther_int: int = 100,
    s2_period: int = 20000,
    s2_ther_int: int = 100,
    s1_ther_int: int = 500,
):
    """
    Generates LAMMPS input files for each nanoparticle.

    Args:
        stage: Simulation stage (0, 1, or 2).
        sim_data_dir: Directory where simulation data is stored.
        init_struct_base_dir: Base directory for initial structure files.
        eam_dir: Directory containing EAM potential files.
        template_dir: Directory containing LAMMPS input templates.
        tnp_types: List of TNP types to process.
    """
    if tnp_types is None:
        tnp_types = ['CL10S/', 'CRALS/', 'CRSR/', 'CS/', 'CSL10/', 'CSRAL/', 'L10R/', 'LL10/', 'RRAL/']

    sim_data_path = Path(sim_data_dir)
    init_struct_base_path = Path(init_struct_base_dir)
    eam_path = Path(eam_dir)
    template_path = Path(template_dir)

    template_file = template_path / f"annealS{stage}.in"
    if not template_file.exists():
        logger.error(f"Template file {template_file} not found!")
        return

    with open(template_file, 'r') as f:
        template_content = f.read()

    for tnp_type in tnp_types:
        tnp_type_dir = sim_data_path / tnp_type.strip('/')
        if not tnp_type_dir.exists():
            logger.warning(f"Directory {tnp_type_dir} not found, skipping...")
            continue

        logger.info(f"Processing {tnp_type} directory...")
        
        # Iterate over each nanoparticle directory in the simulation data directory
        for tnp_dir in tnp_type_dir.iterdir():
            if not tnp_dir.is_dir():
                continue
            
            inp_file_name = tnp_dir.name
            logger.info(f"  Nanoparticle: {inp_file_name}")

            target_in_file = tnp_dir / f"{inp_file_name}S{stage}.in"
            if target_in_file.exists():
                logger.info(f"    {target_in_file} already exists, skipping...")
                continue

            # Identify elements from name (e.g., AuPdPt30_...)
            elements = re.findall(r'[A-Z][a-z]?', inp_file_name)
            if len(elements) < 3:
                logger.warning(f"    Could not identify 3 elements from {inp_file_name}, skipping...")
                continue
            
            element1, element2, element3 = elements[:3]
            logger.info(f"    Elements: {element1}, {element2}, {element3}")

            pot_file = eam_path / "setfl_files" / f"{element1}{element2}{element3}.set"
            init_struct_dir = init_struct_base_path / tnp_type.strip('/')
            
            # Variables for substitution
            subs = {
                "{INP_FILE_NAME}": inp_file_name,
                "{ELEMENT1}": element1,
                "{ELEMENT2}": element2,
                "{ELEMENT3}": element3,
                "{POT_FILE}": str(pot_file),
            }

            if stage == 0:
                s0_dump_int = s0_period // total_dumps
                subs.update({
                    "{INP_DIR_NAME}/": f"{init_struct_dir}/",
                    "{INIT_TEMP}": str(init_temp),
                    "{S0_PERIOD}": str(s0_period),
                    "{S0_THER_INT}": str(s0_ther_int),
                    "{S0_DUMP_INT}": str(s0_dump_int),
                })
            elif stage == 1:
                s1_period = int((heat_temp - init_temp) / heat_rate * 1000)
                s1_dump_int = s1_period // total_dumps
                subs.update({
                    "{INIT_TEMP}": str(init_temp),
                    "{HEAT_TEMP}": str(heat_temp),
                    "{S1_PERIOD}": str(s1_period),
                    "{S1_THER_INT}": str(s1_ther_int),
                    "{S1_DUMP_INT}": str(s1_dump_int),
                })
            elif stage == 2:
                s1_period = int((heat_temp - init_temp) / heat_rate * 1000)
                s1_dump_int = s1_period // total_dumps
                s2_dump_int = s2_period // total_dumps
                subs.update({
                    "{S1_DUMP_INT}": str(s1_dump_int),
                    "{HEAT_TEMP}": str(heat_temp),
                    "{S2_PERIOD}": str(s2_period),
                    "{S2_THER_INT}": str(s2_ther_int),
                    "{S2_DUMP_INT}": str(s2_dump_int),
                })

            # Perform substitution
            content = template_content
            for key, value in subs.items():
                content = content.replace(key, value)

            with open(target_in_file, 'w') as f:
                f.write(content)
            
            logger.info(f"    Generated {target_in_file}")
