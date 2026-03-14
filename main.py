"""CLI entry point for the local ancestry inference benchmarking pipeline."""

import argparse
import sys


def cmd_prep(args):
    from src.data_prep import run_prep

    run_prep(args.chrom)


def cmd_simulate(args):
    from src.simulation import run_simulate

    generations = (
        [int(g) for g in args.generations.split(",")] if args.generations else None
    )
    run_simulate(
        chrom=args.chrom,
        scenario=args.scenario,
        generations=generations,
        seed=args.seed,
    )


def cmd_run(args):
    from src.runners import run_tools

    generations = (
        [int(g) for g in args.generations.split(",")] if args.generations else None
    )
    run_tools(
        chrom=args.chrom,
        tool=args.tool,
        mode=args.mode,
        scenario=args.scenario,
        generations=generations,
    )


def cmd_evaluate(args):
    from src.evaluation import run_evaluate

    run_evaluate(args.chrom)


def cmd_visualize(args):
    from src.visualization import run_visualize

    run_visualize(args.chrom)


def cmd_all(args):
    """Run the full pipeline end-to-end."""
    print("=== Step 1/5: Data Preparation ===")
    cmd_prep(args)

    print("\n=== Step 2/5: Simulation ===")
    cmd_simulate(args)

    print("\n=== Step 3/5: Running LAI Tools ===")
    cmd_run(args)

    print("\n=== Step 4/5: Evaluation ===")
    cmd_evaluate(args)

    print("\n=== Step 5/5: Visualization ===")
    cmd_visualize(args)

    print("\nPipeline complete.")


def main():
    parser = argparse.ArgumentParser(
        description="Local Ancestry Inference Benchmarking Pipeline",
    )
    subparsers = parser.add_subparsers(dest="command", help="Pipeline step to run")

    # Common arguments
    def add_common(p):
        p.add_argument(
            "--chrom", type=int, default=22, help="Chromosome number (default: 22)"
        )

    # prep
    p_prep = subparsers.add_parser(
        "prep",
        help="Prepare data: reformat genetic maps, create sample maps, split VCFs",
    )
    add_common(p_prep)
    p_prep.set_defaults(func=cmd_prep)

    # simulate
    p_sim = subparsers.add_parser(
        "simulate", help="Simulate admixed genomes with ground-truth ancestry"
    )
    add_common(p_sim)
    p_sim.add_argument(
        "--scenario",
        choices=["2way", "3way"],
        default=None,
        help="Admixture scenario (default: all)",
    )
    p_sim.add_argument(
        "--generations",
        type=str,
        default=None,
        help="Comma-separated generation values (default: 5,10,20,50)",
    )
    p_sim.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    p_sim.set_defaults(func=cmd_simulate)

    # run
    p_run = subparsers.add_parser(
        "run", help="Run LAI tools (RFMix/Gnomix) on query data"
    )
    add_common(p_run)
    p_run.add_argument(
        "--tool",
        choices=["rfmix", "gnomix"],
        default=None,
        help="Tool to run (default: both)",
    )
    p_run.add_argument(
        "--mode",
        choices=["sim", "real"],
        default=None,
        help="Data mode (default: both)",
    )
    p_run.add_argument(
        "--scenario",
        choices=["2way", "3way"],
        default=None,
        help="Admixture scenario (default: all)",
    )
    p_run.add_argument(
        "--generations",
        type=str,
        default=None,
        help="Comma-separated generation values (default: 5,10,20,50)",
    )
    p_run.set_defaults(func=cmd_run)

    # evaluate
    p_eval = subparsers.add_parser(
        "evaluate", help="Compute accuracy, concordance, and global ancestry metrics"
    )
    add_common(p_eval)
    p_eval.set_defaults(func=cmd_evaluate)

    # visualize
    p_viz = subparsers.add_parser("visualize", help="Generate all benchmark plots")
    add_common(p_viz)
    p_viz.set_defaults(func=cmd_visualize)

    # all
    p_all = subparsers.add_parser("all", help="Run the full pipeline end-to-end")
    add_common(p_all)
    p_all.add_argument(
        "--scenario",
        choices=["2way", "3way"],
        default=None,
        help="Admixture scenario (default: all)",
    )
    p_all.add_argument(
        "--generations",
        type=str,
        default=None,
        help="Comma-separated generation values (default: 5,10,20,50)",
    )
    p_all.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    p_all.add_argument(
        "--tool",
        choices=["rfmix", "gnomix"],
        default=None,
        help="Tool to run (default: both)",
    )
    p_all.add_argument(
        "--mode",
        choices=["sim", "real"],
        default=None,
        help="Data mode (default: both)",
    )
    p_all.set_defaults(func=cmd_all)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
