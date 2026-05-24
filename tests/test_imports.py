"""Smoke tests that the package modules are importable."""


def test_import_constants():
    from allomorph.constants import (
        ELE_DICT,
        SHAPE_LIST,
    )

    assert "Au" in ELE_DICT
    assert "Pd" in ELE_DICT
    assert "Pt" in ELE_DICT
    assert "TH" in SHAPE_LIST
    assert "RD" in SHAPE_LIST
    assert "TO" in SHAPE_LIST
    assert "CO" in SHAPE_LIST


def test_import_eam():
    from allomorph.eam import Database

    assert "Au" in Database
    assert "Pd" in Database
    assert "Pt" in Database


def test_import_init_struct():
    from allomorph.init_struct import (
        gen_bnp,
        gen_mnp,
        gen_tnp,
        write_bnp,
        write_mnp,
        write_tnp,
    )

    assert callable(gen_mnp)
    assert callable(write_mnp)
    assert callable(gen_bnp)
    assert callable(write_bnp)
    assert callable(gen_tnp)
    assert callable(write_tnp)


def test_import_feat_ext_eng():
    from allomorph.feat_ext_eng import (
        concat_np_feats,
        merge_reformat_data,
        run_ncpac,
        setup_ncpac,
    )

    assert callable(setup_ncpac)
    assert callable(run_ncpac)
    assert callable(merge_reformat_data)
    assert callable(concat_np_feats)


def test_cli_import():
    from allomorph.cli import main

    assert callable(main)
