"""Command-line interface for np-gen."""

import argparse
import sys

from np_gen.constants import load_config, parse_ele_comb, update_constants


def main(argv=None):
    """Entry point for the np-gen CLI."""
    parser = argparse.ArgumentParser(
        prog="np-gen",
        description="Toolkit for generating monometallic to trimetallic nanoparticle structural datasets.",
    )
    parser.add_argument(
        "--config",
        help="Path to a configuration file (JSON, YAML, or TOML) to override default constants.",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # EAM subcommand
    eam_parser = subparsers.add_parser(
        "eam",
        help="Create EAM alloy potential files.",
    )
    eam_parser.add_argument(
        "-n", "--names",
        dest="name",
        nargs="+",
        required=True,
        help="Element names (e.g., Au Pt Pd)",
    )
    eam_parser.add_argument(
        "-nr",
        dest="nr",
        type=int,
        default=2000,
        help="Number of points in r space (default: 2000)",
    )
    eam_parser.add_argument(
        "-nrho",
        dest="nrho",
        type=int,
        default=2000,
        help="Number of points in rho space (default: 2000)",
    )
    eam_parser.add_argument(
        "-o", "--output-dir",
        dest="output_dir",
        default="setfl_files",
        help="Output directory for .set files (default: setfl_files)",
    )
    eam_parser.set_defaults(func=_eam_cmd)

    # init-struct subcommand
    init_parser = subparsers.add_parser(
        "init-struct",
        help="Generate initial nanoparticle structures.",
    )
    init_parser.add_argument(
        "--stage",
        choices=["mnp", "bnp", "tnp", "cs", "all"],
        default="all",
        help="Which structure generation stage to run (default: all)",
    )
    init_parser.add_argument(
        "--replace",
        action="store_true",
        help="Overwrite existing files.",
    )
    init_parser.add_argument(
        "--vis",
        action="store_true",
        help="Visualise generated structures (opens ASE GUI).",
    )
    init_parser.set_defaults(func=_init_struct_cmd)

    # md-sim subcommand
    md_parser = subparsers.add_parser(
        "md-sim",
        help="Manage MD simulations.",
    )
    md_subparsers = md_parser.add_subparsers(dest="action", help="MD actions")

    # md-sim setup
    md_setup_parser = md_subparsers.add_parser("setup", help="Set up simulation directories.")
    md_setup_parser.add_argument("--init-dir", required=True, help="Path to initial structures.")
    md_setup_parser.add_argument("--target-dir", required=True, help="Target simulation directory.")

    # md-sim gen-input
    md_gen_parser = md_subparsers.add_parser("gen-input", help="Generate LAMMPS input files.")
    md_gen_parser.add_argument("--stage", type=int, required=True, help="Simulation stage (0, 1, 2).")
    md_gen_parser.add_argument("--sim-dir", required=True, help="Simulation data directory.")
    md_gen_parser.add_argument("--init-dir", required=True, help="Initial structures base directory.")
    md_gen_parser.add_argument("--eam-dir", required=True, help="EAM potential directory.")
    md_gen_parser.add_argument("--template-dir", required=True, help="LAMMPS templates directory.")

    # md-sim gen-joblist
    md_joblist_parser = md_subparsers.add_parser("gen-joblist", help="Generate list of jobs to run.")
    md_joblist_parser.add_argument("--stage", type=int, required=True, help="Simulation stage.")
    md_joblist_parser.add_argument("--sim-dir", required=True, help="Simulation data directory.")
    md_joblist_parser.add_argument("--joblist", default="jobList", help="Output job list file.")

    # md-sim submit
    md_submit_parser = md_subparsers.add_parser("submit", help="Submit jobs to PBS queue.")
    md_submit_parser.add_argument("--joblist", default="jobList", help="Input job list file.")
    md_submit_parser.add_argument("--sim-dir", required=True, help="Simulation data directory.")
    md_submit_parser.add_argument("--project", required=True, help="PBS project code.")
    md_submit_parser.add_argument("--max-queue", type=int, default=15, help="Max jobs in queue.")
    md_submit_parser.add_argument("--email", help="Notification email.")

    md_parser.set_defaults(func=_md_sim_cmd)

    # feat-ext subcommand
    feat_parser = subparsers.add_parser(
        "feat-ext",
        help="Manage feature extraction.",
    )
    feat_subparsers = feat_parser.add_subparsers(dest="action", help="Feature extraction actions")

    # feat-ext setup
    feat_setup_parser = feat_subparsers.add_parser("setup", help="Set up NCPac directories.")
    feat_setup_parser.add_argument("--sim-dir", required=True, help="Simulation data directory.")
    feat_setup_parser.add_argument("--target-dir", required=True, help="Target directory for NCPac.")
    feat_setup_parser.add_argument("--exe", required=True, help="Path to NCPac executable.")
    feat_setup_parser.add_argument("--inp", required=True, help="Path to NCPac template input file.")
    feat_setup_parser.add_argument("--ele-comb", help="Element combination (e.g. AuPtPd). If omitted, inferred from sim-dir.")

    # feat-ext run
    feat_run_parser = feat_subparsers.add_parser("run", help="Run NCPac in parallel.")
    feat_run_parser.add_argument("--target-dir", required=True, help="Directory where NCPac is set up.")
    feat_run_parser.add_argument("--final-dir", required=True, help="Directory to store final features.")

    # feat-ext merge
    feat_merge_parser = feat_subparsers.add_parser("merge", help="Merge MD outputs with NCPac features.")
    feat_merge_parser.add_argument("--md-out", required=True, help="Path to MDout.csv.")
    feat_merge_parser.add_argument("--feat-dir", required=True, help="Directory with raw NCPac features.")
    feat_merge_parser.add_argument("--output-dir", required=True, help="Directory to store merged features.")
    feat_merge_parser.add_argument("--ele-comb", help="Element combination (e.g. AuPtPd).")

    # feat-ext concat
    feat_concat_parser = feat_subparsers.add_parser("concat", help="Concatenate all merged features into one CSV.")
    feat_concat_parser.add_argument("--input-dir", required=True, help="Directory with merged features.")
    feat_concat_parser.add_argument("--output-file", required=True, help="Output CSV file path.")

    feat_parser.set_defaults(func=_feat_ext_cmd)

    args = parser.parse_args(argv)

    if args.config:
        config = load_config(args.config)
        update_constants(config)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    args.func(args)


def _eam_cmd(args):
    """Run the EAM creation command."""
    from np_gen.eam.create_eam import create_eam as _create_eam
    argv = []
    argv.extend(["-n", *args.name])
    argv.extend(["-nr", str(args.nr)])
    argv.extend(["-nrho", str(args.nrho)])
    argv.extend(["-o", args.output_dir])
    _create_eam(argv)


def _init_struct_cmd(args):
    """Run the initial structure generation command."""
    from np_gen.init_struct.gen_bnp_al import main as gen_bnp_main
    from np_gen.init_struct.gen_bnp_cs import write_hard_core_shell as gen_bnp_cs_main
    from np_gen.init_struct.gen_mnp import main as gen_mnp_main
    from np_gen.init_struct.gen_tnp_al import main as gen_tnp_main

    if args.stage in ("mnp", "all"):
        print("=== Generating monometallic nanoparticles (MNP) ===")
        gen_mnp_main(replace=args.replace, vis=args.vis)
    if args.stage in ("bnp", "all"):
        print("=== Generating bimetallic nanoparticles (BNP) ===")
        gen_bnp_main(replace=args.replace, vis=args.vis)
    if args.stage in ("cs", "all"):
        print("=== Generating hard core-shell nanoparticles (CS) ===")
        gen_bnp_cs_main(replace=args.replace, vis=args.vis)
    if args.stage in ("tnp", "all"):
        print("=== Generating trimetallic nanoparticles (TNP) ===")
        gen_tnp_main(replace=args.replace, vis=args.vis)


def _md_sim_cmd(args):
    """Run the MD simulation management commands."""
    from np_gen.md_sim.generator import generate_lammps_input
    from np_gen.md_sim.manager import generate_job_list
    from np_gen.md_sim.setup import setup_md_sim
    from np_gen.md_sim.submission import submit_jobs

    if args.action == "setup":
        setup_md_sim(args.init_dir, args.target_dir)
    elif args.action == "gen-input":
        generate_lammps_input(
            args.stage, args.sim_dir, args.init_dir, args.eam_dir, args.template_dir
        )
    elif args.action == "gen-joblist":
        generate_job_list(args.stage, args.sim_dir, args.joblist)
    elif args.action == "submit":
        submit_jobs(args.joblist, args.sim_dir, args.project, args.max_queue, args.email)
    else:
        print("Unknown md-sim action.")


def _feat_ext_cmd(args):
    """Run the feature extraction management commands."""
    from np_gen.feat_ext_eng.gen_csvs import run_ncpac_parallel, setup_ncpac
    from np_gen.feat_ext_eng.merge_features import (
        concat_np_feats,
        generate_headers,
        run_merge_reformat_parallel,
    )

    if args.action == "setup":
        setup_ncpac(args.sim_dir, args.target_dir, args.exe, args.inp, ele_comb=args.ele_comb)
    elif args.action == "run":
        from pathlib import Path
        target_path = Path(args.target_dir)
        working_list = []
        if target_path.exists():
            for d in target_path.iterdir():
                if d.is_dir() and (d / "NCPac.exe").exists() and not (d / "DONE.txt").exists():
                    working_list.append((str(d), d.name))
        run_ncpac_parallel(working_list, args.final_dir)
    elif args.action == "merge":
        if args.ele_comb is None:
            # Try to infer from md-out's directory name
            args.ele_comb = Path(args.md_out).parent.name
            if not args.ele_comb or len(args.ele_comb) < 2:
                args.ele_comb = "AuPtPd"
        elements = parse_ele_comb(args.ele_comb)
        headers_list = generate_headers(elements)
        run_merge_reformat_parallel(
            args.md_out, args.feat_dir, args.output_dir, args.ele_comb, headers_list
        )
    elif args.action == "concat":
        concat_np_feats(args.feat_dir, args.output_file)
    else:
        print("Unknown feat-ext action.")


if __name__ == "__main__":
    main()
