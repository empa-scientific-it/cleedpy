import numpy as np
import typer

from ..rfactor import compute_rfactor


def rfactor(
    experimental_iv: str = typer.Option(  # noqa: B008
        "experimental.txt", "--experimental", "-e", help="Experimental IV curves"
    ),
    theoretical_iv: str = typer.Option(  # noqa: B008
        "theoretical.txt", "--theoretical", "-t", help="Theoretical IV curves"
    ),
    rfactor_type: str = typer.Option(  # noqa: B008
        "r2", "--rfactor", "-r", help="Rfactor type: r2 or pendry"
    ),
):
    """Command line interface for the rfactor tool"""
    exper_iv = np.loadtxt(experimental_iv)
    theor_iv = np.loadtxt(theoretical_iv)

    r = compute_rfactor(
        theoretical_iv=theor_iv, experimental_iv=exper_iv, rfactor_type=rfactor_type
    )

    print(f"Rfactor={r}")


def cli():
    """Leed CLI."""
    typer.run(rfactor)


if __name__ == "__main__":
    cli()
