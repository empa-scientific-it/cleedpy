from pathlib import Path

import typer

from .. import search
from ..config import OLD_FORMAT_TEMPLATE, load_parameters


def parse_atoms_to_optimize(value: str | None) -> set[str] | None:
    """Parse the atoms to optimize from a string."""
    if value is None:
        return None

    indices = []
    parts = value.split(",")
    for part in parts:
        if "-" in part:
            start, end = part.split("-")
            indices.extend(range(int(start) - 1, int(end)))  # Convert to 0-based index
        else:
            indices.append(int(part) - 1)  # Convert to 0-based index

    return set(indices)


def run_search(
    parameters_file: Path = typer.Option(  # noqa: B008
        "leed.yaml",
        "--input",
        "-i",
        file_okay=True,
        dir_okay=False,
        readable=True,
        exists=True,
        help="Input file with parameters. Accepts ONLY .yaml",
    ),
    phase_path: Path = typer.Option(  # noqa: B008
        "PHASE",
        "--phase",
        "-p",
        file_okay=False,
        dir_okay=True,
        readable=True,
        exists=True,
        help="Directory with phase shifts.",
    ),
    experimental_iv: Path = typer.Option(  # noqa: B008
        "experimental.txt",
        "--experimental-iv",
        "-e",
        file_okay=True,
        dir_okay=False,
        readable=True,
        exists=True,
        help="File with experimental IV curves",
    ),
    optimization_axes: str = typer.Option(  # noqa: B008
        "xyz",
        "--optimization-axes",
        "-a",
        help="Which coordinate axes to optimize. Options: x, y, z, or combinations like 'xy', 'yz', 'xyz'.",
    ),
    optimize_shift: bool = typer.Option(  # noqa: B008
        True,
        "--optimize-shift/--no-optimize-shift",
        help="Whether to optimize the overall shift of the IV curves",
    ),
    atoms_to_optimize: str = typer.Option(  # noqa: B008
        None,
        "--atoms-to-optimize",
        "-A",
        help="Which overlayer (!) atoms to optimize, e.g. '1-20,24-30'. The indexing starts from 1.",
    ),
) -> None:
    """Command line interface for the search tool"""

    if parameters_file.suffix in [".yml", ".yaml"]:
        config = load_parameters(parameters_file)
    else:
        message = f"Input file {parameters_file} has unsupported format, only .yml/.yaml are supported."
        raise ValueError(message)

    old_format = OLD_FORMAT_TEMPLATE.render(**config.model_dump())
    with open("initial_parameters.inp", "w") as f:
        f.write(old_format)

    searcher = search.CleedSearchCoordinator(
        config=config,
        phase_path=str(phase_path),
        experimental_iv_file=str(experimental_iv),
    )

    # Prepare what to optimize.
    atoms_to_optimize = parse_atoms_to_optimize(atoms_to_optimize)
    if atoms_to_optimize is None:  # If None, optimize all atoms
        atoms_to_optimize = set(range(len(config.overlayers)))
    atoms_and_axes = [
        optimization_axes if i in atoms_to_optimize else ""
        for i in range(len(config.overlayers))
    ]

    searcher.set_search_parameters(
        overlayer_atoms=atoms_and_axes, optimize_shift=optimize_shift
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
