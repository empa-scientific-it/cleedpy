from pathlib import Path

import numpy as np
import typer

from .. import search
from ..config import OLD_FORMAT_TEMPLATE, load_parameters


def link_parameters_to_x_init(parameters) -> np.ndarray:
    """Link the parameters to the x_init array for the optimization."""

    x_init = []
    for atom in parameters.overlayers:
        x_init.extend([atom.position.z])

    return np.array(x_init)


def run_search(
    parameters_file: str = typer.Option(  # noqa: B008
        "leed.yaml",
        "--input",
        "-i",
        help="Input file with parameters. Accepts ONLY .yaml",
    ),
    phase_path: str = typer.Option(  # noqa: B008
        "PHASE", "--phase", "-p", help="Phase path"
    ),
    output_file: str = typer.Option(  # noqa: B008
        "leed.out", "--output", "-o", help="Output file"
    ),  # noqa: B008
    optimize_overlayer_atoms: str = typer.Option(  # noqa: B008
        "xyz",
        "--optimize-overlayer-atoms",
        "-a",
        help="Which overlayer atom coordinates to optimize. Options: x, y, z, or any combination of them, e.g. 'xy', 'yz', 'xyz'",
    ),
    optimize_shift: bool = typer.Option(  # noqa: B008
        True,
        "--optimize-shift/--no-optimize-shift",
        "-s/-S",
        help="Whether to optimize the overall shift of the IV curves",
    ),
):
    """Command line interface for the search tool"""
    parameters_file = Path(parameters_file)
    if parameters_file.suffix in [".yml", ".yaml"]:
        config = load_parameters(parameters_file)
    else:
        message = f"Input file {parameters_file} has unsupported format, only .yml/.yaml are supported."
        raise ValueError(message)

    old_format = OLD_FORMAT_TEMPLATE.render(**config.model_dump())
    with open("initial_parameters.inp", "w") as f:
        f.write(old_format)

    searcher = search.CleedSearchCoordinator(config=config, phase_path=phase_path)
    searcher.set_search_parameters(
        overlayer_atoms=optimize_overlayer_atoms, optimize_shift=optimize_shift
    )
    searcher.start_optimization(method="Nelder-Mead")

    print("Optimization result:", searcher.result)
    print("Final parameters:", searcher.x)
    print("Final R-factor:", searcher.function_to_minimize(searcher.result.x))


def cli():
    """Search CLI."""
    typer.run(run_search)


if __name__ == "__main__":
    cli()
