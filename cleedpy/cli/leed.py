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


def print_cleed_results(result):
    """Function that prints the results of cleed.
        The format looks as follows
    # ####################################### #
    #            output from CLEED            #
    # ####################################### #
    #vn 0.0 (test version GH/11.08.95)
    #ts Mon Jul 29 14:47:35 2024
    #
    #en 108 70.000000 498.100000 4.000000
    #bn 55
    #bi 0 0.000000 0.000000 0
    #bi 1 -1.000000 0.000000 0
    #bi 2 -1.000000 1.000000 0
    #bi 3 0.000000 -1.000000 0
    ...
    70.00 4.364253e-03 5.741584e-03 7.114631e-03 7.115966e-03 5.744055e-03 5.741942e-03 7.113053e-03 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00

    """
    with open("leed.out", "w") as f:
        f.write("# ####################################### #\n")
        f.write("#            output from CLEED            #\n")
        f.write("# ####################################### #\n")
        f.write(f"#vn {__version__}\n")
        f.write(f"#ts {datetime.datetime.now().ctime()}\n")
        f.write("#\n")
        energy_step = (result.energies[result.n_energies - 1] - result.energies[0]) / (
            result.n_energies - 1
        )
        f.write(
            f"#en {result.n_energies} {result.energies[0]*HART:.6f} {result.energies[result.n_energies-1]*HART:.6f} {energy_step*HART:.6f}\n"
        )
        f.write(f"#bn {result.n_beams}\n")
        for i in range(result.n_beams):
            f.write(
                f"#bi {i} {result.beam_x[i]:.6f} {result.beam_y[i]:.6f} {result.beam_set[i]}\n"
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
            f.write(" \n")


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
    result = call_cleed(parameters_file, bulk_file, output_file)

    print_cleed_results(result)


def cli():
    """Leed CLI."""
    typer.run(leed)


if __name__ == "__main__":
    cli()
