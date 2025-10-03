from pathlib import Path

import numpy as np
import typer
from scipy import optimize

from .. import rfactor
from ..config import OLD_FORMAT_TEMPLATE, load_parameters
from ..interface.cleed import call_cleed
from ..physics import constants


def link_parameters_to_x_init(parameters) -> np.ndarray:
    """Link the parameters to the x_init array for the optimization."""

    x_init = []
    for atom in parameters.overlayers:
        x_init.extend([atom.position.z])

    return np.array(x_init)


def cleed_result_to_iv(result) -> np.ndarray:
    """Convert CLEED result to IV array."""
    iv = []
    n_beams = result.n_beams
    n_energies = result.n_energies
    for i in range(n_energies):
        for j in range(n_beams):
            iv.append(
                [
                    result.beam_index1[j],
                    result.beam_index2[j],
                    0,
                    result.energies[i] * constants.HART,
                    result.iv_curves[i * n_beams + j],
                ]
            )
    return np.array(iv)


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
        initial_parameters_file = parameters_file.with_stem(
            parameters_file.stem + "_init"
        ).with_suffix(".inp")
        current_parameters_file = parameters_file.with_suffix(".inp")
        with open(initial_parameters_file, "w") as f:
            f.write(old_format)
    else:
        message = f"Input file {parameters_file} has unsupported format, only .yml/.yaml are supported."
        raise ValueError(message)

    # Some useful variables.
    iteration = 0
    largest_rfactor = 1.0

    def f(x):
        nonlocal iteration, largest_rfactor

        # Get x and update parameters
        for i, atom in enumerate(config.overlayers):
            atom.position.z = x[i]

        geometrical_r = rfactor.compute_geometrical_rfactor(config)

        if geometrical_r > 1.0:
            return geometrical_r + largest_rfactor

        old_format = OLD_FORMAT_TEMPLATE.render(**config.model_dump())
        with open(current_parameters_file, "w") as f:
            f.write(old_format)

        result = call_cleed(
            str(current_parameters_file), str(current_parameters_file), phase_path
        )
        iteration += 1
        theor_iv = cleed_result_to_iv(result)
        np.savetxt(
            output_file,
            theor_iv,
            header="1st beam index, 2nd beam index, beam set, energy, intensity",
            fmt=["%4d", "%4d", "%4d", "%8.2f", "%.10e"],
        )

        exp_iv = np.loadtxt("experimental.txt")

        optimize_shift_iteration = 0

        def for_shift_optimization(shift):
            nonlocal optimize_shift_iteration
            optimize_shift_iteration += 1
            r = rfactor.compute_rfactor(
                theoretical_iv=np.array(theor_iv),
                experimental_iv=exp_iv,
                shift=shift,
                rfactor_type="pendry",
            )
            with open("search.log", "a") as f:
                f.write(
                    f"Shift optimization iteration {optimize_shift_iteration}: Shift={shift} Rfactor={r}\n"
                )

            return r

        res = optimize.minimize_scalar(
            for_shift_optimization, bounds=(-10, 10), method="bounded"
        )
        with open("search.log", "a") as f:
            f.write(f"Optimal shift found: {res.x} with Rfactor={res.fun}\n")
        iv_r = rfactor.compute_rfactor(
            theoretical_iv=np.array(theor_iv),
            experimental_iv=exp_iv,
            rfactor_type="pendry",
            shift=res.x,
        )

        largest_rfactor = max(iv_r, largest_rfactor)

        with open("search.log", "a") as f:
            f.write(f"Iteration {iteration}: Rfactor={iv_r}\n")

        return iv_r + geometrical_r

    x_init = link_parameters_to_x_init(config)

    result = optimize.minimize(f, x_init, method="Nelder-Mead", tol=5e-4)
    print("Optimization result:", result.x)


def cli():
    """Search CLI."""
    typer.run(run_search)


if __name__ == "__main__":
    cli()
