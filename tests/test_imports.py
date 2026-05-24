"""Smoke tests that the package modules are importable."""


def test_import_constants():
    from tnp_gen.constants import ELE_DICT, DIAMETER_LIST, SHAPE_LIST

    assert "Au" in ELE_DICT
    assert "Pd" in ELE_DICT
    assert "Pt" in ELE_DICT


def test_import_eam():
    from tnp_gen.eam import AtType, Database, create_eam

    assert "Au" in Database
    assert "Pd" in Database
    assert "Pt" in Database


def test_import_init_struct():
    from tnp_gen.init_struct import (
        gen_mnp,
        write_mnp,
        gen_bnp,
        write_bnp,
        gen_tnp,
        write_tnp,
    )

    assert callable(gen_mnp)
    assert callable(write_mnp)
    assert callable(gen_bnp)
    assert callable(write_bnp)
    assert callable(gen_tnp)
    assert callable(write_tnp)


def test_import_feat_ext_eng():
    from tnp_gen.feat_ext_eng import (
        setup_ncpac,
        run_ncpac,
        merge_reformat_data,
        concat_np_feats,
    )

    assert callable(setup_ncpac)
    assert callable(run_ncpac)
    assert callable(merge_reformat_data)
    assert callable(concat_np_feats)


def test_cli_import():
    from tnp_gen.cli import main

    assert callable(main)
