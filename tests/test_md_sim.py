from np_gen.md_sim.generator import generate_lammps_input
from np_gen.md_sim.manager import generate_job_list, get_config_value, update_config
from np_gen.md_sim.setup import setup_md_sim


def test_setup_md_sim(tmp_path):
    # Create mock InitStruct directory
    init_dir = tmp_path / "InitStruct" / "TNP" / "Type1"
    init_dir.mkdir(parents=True)
    (init_dir / "NP1.lmp").touch()

    target_dir = tmp_path / "MDsim"

    setup_md_sim(str(tmp_path / "InitStruct" / "TNP"), str(target_dir))

    assert (target_dir / "Type1" / "NP1" / "NP1.lmp").exists()
    assert (target_dir / "Type1" / "NP1" / "config.yml").exists()

def test_generate_lammps_input(tmp_path):
    # Create mock directories and files
    sim_dir = tmp_path / "MDsim"
    np_dir = sim_dir / "Type1" / "AuPdPt30_test"
    np_dir.mkdir(parents=True)
    (np_dir / "AuPdPt30_test.lmp").write_text("atoms 100")

    init_dir = tmp_path / "InitStruct"
    eam_dir = tmp_path / "EAM" / "setfl_files"
    eam_dir.mkdir(parents=True)

    template_dir = tmp_path / "templates"
    template_dir.mkdir()
    (template_dir / "annealS0.in").write_text("name {INP_FILE_NAME} pot {POT_FILE} elements {ELEMENT1} {ELEMENT2} {ELEMENT3}")

    generate_lammps_input(
        stage=0,
        sim_data_dir=str(sim_dir),
        init_struct_base_dir=str(init_dir),
        eam_dir=str(tmp_path / "EAM"),
        template_dir=str(template_dir),
        tnp_types=["Type1/"]
    )

    in_file = np_dir / "AuPdPt30_testS0.in"
    assert in_file.exists()
    content = in_file.read_text()
    assert "AuPdPt30_test" in content
    assert "Au" in content
    assert "Pd" in content
    assert "Pt" in content


def test_generate_lammps_input_bimetallic(tmp_path):
    """Test that LAMMPS input generation works for bimetallic names too."""
    sim_dir = tmp_path / "MDsim"
    np_dir = sim_dir / "BNP" / "AuPd30_test"
    np_dir.mkdir(parents=True)
    (np_dir / "AuPd30_test.lmp").write_text("atoms 80")

    init_dir = tmp_path / "InitStruct"
    template_dir = tmp_path / "templates"
    template_dir.mkdir()
    (template_dir / "annealS0.in").write_text("name {INP_FILE_NAME} elements {ELEMENT1} {ELEMENT2}")

    generate_lammps_input(
        stage=0,
        sim_data_dir=str(sim_dir),
        init_struct_base_dir=str(init_dir),
        eam_dir=str(tmp_path / "EAM"),
        template_dir=str(template_dir),
        tnp_types=["BNP/"]
    )

    in_file = np_dir / "AuPd30_testS0.in"
    assert in_file.exists()
    content = in_file.read_text()
    assert "AuPd30_test" in content
    assert "Au" in content
    assert "Pd" in content

def test_config_management(tmp_path):
    config_file = tmp_path / "config.yml"
    update_config(config_file, "S0eq", "true")
    assert get_config_value(config_file, "S0eq") == "true"

    update_config(config_file, "S1ok", "false")
    assert get_config_value(config_file, "S1ok") == "false"

def test_generate_lammps_input_cu_ag(tmp_path):
    """Test LAMMPS input generation for CuAg nanoparticle."""
    sim_dir = tmp_path / "MDsim"
    np_dir = sim_dir / "BNP" / "CuAg20_test"
    np_dir.mkdir(parents=True)
    (np_dir / "CuAg20_test.lmp").write_text("atoms 50")

    init_dir = tmp_path / "InitStruct"
    template_dir = tmp_path / "templates"
    template_dir.mkdir()
    (template_dir / "annealS0.in").write_text("mapping {MAPPING}")

    generate_lammps_input(
        stage=0,
        sim_data_dir=str(sim_dir),
        init_struct_base_dir=str(init_dir),
        eam_dir=str(tmp_path / "EAM"),
        template_dir=str(template_dir),
        tnp_types=["BNP/"]
    )

    in_file = np_dir / "CuAg20_testS0.in"
    assert in_file.exists()
    content = in_file.read_text()
    assert "Cu Ag" in content
