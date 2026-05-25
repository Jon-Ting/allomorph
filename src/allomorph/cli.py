"""Command-line interface for AlloMorph."""

import argparse
import sys

from allomorph.constants import load_config, update_constants


def main(argv=None):
    """Entry point for the AlloMorph CLI."""
    parser = argparse.ArgumentParser(
        prog="allomorph",
        description="Toolkit for generating monometallic to trimetallic nanoparticle structural datasets.",
    )
    parser.add_argument(
        "--config",
        help="Path to a configuration file (JSON, YAML, or TOML) to override default constants.",
    )
    
    # Subcommands
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # init-struct subcommand (Main functionality)
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

    args = parser.parse_args(argv)

    if args.config:
        config = load_config(args.config)
        update_constants(config)

    if args.vis:
        print("Warning: --vis flag enabled. Visualization will run serially to prevent system hang.")
        print("Many windows may be opened sequentially. Close one to see the next.")

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    args.func(args)


def _init_struct_cmd(args):
    """Run the initial structure generation command."""
    from allomorph.init_struct.gen_bnp_al import main as gen_bnp_main
    from allomorph.init_struct.gen_bnp_cs import write_hard_core_shell as gen_bnp_cs_main
    from allomorph.init_struct.gen_mnp import main as gen_mnp_main
    from allomorph.init_struct.gen_tnp_al import main as gen_tnp_main

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


if __name__ == "__main__":
    main()
