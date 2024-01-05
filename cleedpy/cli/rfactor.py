import click
import yaml
import numpy as np
from .. import rfactor

@click.command("cli")
@click.option("--input", "-i", help="Input file", required=True)
def cli(input):
    """Command line interface for the rfactor tool"""
    with open(input, "r") as f:
        data = yaml.safe_load(f)

    y_true = np.array(data["y_true"])
    y_pred = np.array(data["y_pred"])

    r_factor = rfactor.mean_square_error(y_true=y_true, y_pred=y_pred)

    print(f"R factor: {r_factor}")


if __name__ == "__main__":
    cli()