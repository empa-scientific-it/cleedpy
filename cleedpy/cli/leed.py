from pathlib import Path

import typer

from ..config import load_parameters
from ..interface.cleed import call_cleed
from ..physics.constants import BOHR_TO_ANGSTROM
from ..preprocess import extract_bulk_parameters, transform_cells


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
    print(result)


def cli():
    """Leed CLI."""
    typer.run(leed)


if __name__ == "__main__":
    cli()
