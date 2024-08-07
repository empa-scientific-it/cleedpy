import datetime
from pathlib import Path

import typer

from ..config import load_parameters
from ..interface.cleed import call_cleed
from ..physics.constants import BOHR_TO_ANGSTROM, HART
from ..preprocess import extract_bulk_parameters, transform_cells
from ..version import __version__


def leed_bak(inptut_file: Path):
    """Leed CLI."""
    config = load_parameters(inptut_file)
    transformation_matrix = transform_cells(config)

    bulk_parameters = extract_bulk_parameters(config, transformation_matrix)
    print(
        f"Optical potential: {bulk_parameters.vr} eV (Vr), {bulk_parameters.vi} eV (Vi)"
    )
    print(f"Temperature: {bulk_parameters.temp} K")
    print(
        f"""Bulk 2-dim. unit cell: {bulk_parameters.a[0]}
    {bulk_parameters.a[1]*BOHR_TO_ANGSTROM} {bulk_parameters.a[3]*BOHR_TO_ANGSTROM}
    {bulk_parameters.a[2]*BOHR_TO_ANGSTROM} {bulk_parameters.a[4]*BOHR_TO_ANGSTROM}"""
    )


CLEED_OUT_HEADER = f"""# ####################################### #
#            output from CLEED            #
# ####################################### #
# vn {__version__}
# ts {datetime.datetime.now().ctime()}
#
"""


def print_cleed_results(result, output_file):
    """Function that prints the results of cleed."""
    with open(output_file, "w") as f:
        f.write(CLEED_OUT_HEADER)
        energy_step = (result.energies[result.n_energies - 1] - result.energies[0]) / (
            result.n_energies - 1
        )
        f.write(
            f"#en {result.n_energies} {result.energies[0]*HART:.6f} {result.energies[result.n_energies-1]*HART:.6f} {energy_step*HART:.6f}\n"
        )
        f.write(f"#bn {result.n_beams}\n")
        for i in range(result.n_beams):
            f.write(
                f"#bi {i} {result.beam_index1[i]:.6f} {result.beam_index2[i]:.6f} {result.beam_set[i]}\n"
            )
        for i in range(result.n_energies):
            f.write(
                " ".join(
                    [f"{result.energies[i]*HART:.2f}"]
                    + [
                        f"{result.iv_curves[i*result.n_beams+j]:.6e}"
                        for j in range(result.n_beams)
                    ]
                )
            )
            f.write("\n")


def leed(
    parameters_file: str = typer.Option(  # noqa: B008
        "leed.inp", "--input", "-i", help="Input file with parameters"
    ),
    bulk_file: str = typer.Option(None, "--bulk", "-b", help="Bulk file"),  # noqa: B008
    output_file: str = typer.Option(  # noqa: B008
        "leed.out", "--output", "-o", help="Output file"
    ),  # noqa: B008
):
    """Leed CLI."""
    result = call_cleed(parameters_file, bulk_file)

    print_cleed_results(result, output_file)


def cli():
    """Leed CLI."""
    typer.run(leed)


if __name__ == "__main__":
    cli()
