from pathlib import Path

import numpy as np
import typer

from .. import search
from ..config import OLD_FORMAT_TEMPLATE, load_parameters
from ..interface.cleed import call_cleed
from ..rfactor import compute_rfactor


def link_parameters_to_x_init(parameters) -> np.ndarray:
    """Link the parameters to the x_init array for the optimization."""

    x_init = []
    for atom in parameters.overlayers:
        x_init.extend(atom.position)

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
):
    """Command line interface for the search tool"""
    parameters_file = Path(parameters_file)
    if parameters_file.suffix in [".yml", ".yaml"]:
        config = load_parameters(parameters_file)
        old_format = OLD_FORMAT_TEMPLATE.render(**config.model_dump())
        parameters_file = parameters_file.with_suffix(".inp")
        with open(parameters_file, "w") as f:
            f.write(old_format)
    else:
        message = f"Input file {parameters_file} has unsupported format, only .yml/.yaml are supported."
        raise ValueError(message)
    iteration = 0

    def f(x):
        nonlocal iteration
        result = call_cleed(str(parameters_file), str(parameters_file), phase_path)
        iteration += 1
        theor_iv = []
        with open(output_file, "w") as f:
            n_beams = result.n_beams
            n_energies = result.n_energies
            for i in range(n_energies):
                for j in range(n_beams):
                    theor_iv.append(
                        [
                            result.beam_index1[j],
                            result.beam_index2[j],
                            0,
                            result.energies[i] * 27.2114,
                            result.iv_curves[i * n_beams + j],
                        ]
                    )
        exp_iv = np.loadtxt("experimental.txt")
        r = compute_rfactor(
            theoretical_iv=np.array(theor_iv),
            experimental_iv=exp_iv,
            rfactor_type="pendry",
        )

        with open("search.log", "a") as f:
            f.write(f"Iteration {iteration}: Rfactor={r}\n")
        return r

    x_init = link_parameters_to_x_init(config)
    result = search.simplex(f, x_init)
    print("Optimization result:", result)


def cli():
    """Search CLI."""
    typer.run(run_search)


if __name__ == "__main__":
    cli()
