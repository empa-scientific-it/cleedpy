import click
import numpy as np
import yaml

from .. import rfactor


@click.command("cli")
@click.option("--config", "-c", help="Config file", required=True)
def cli(config):
    """Command line interface for the rfactor tool"""
    with open(config) as f:
        data = yaml.safe_load(f)

    y_true = np.array(data["y_true"])
    y_pred = np.array(data["y_pred"])

    r_factor = rfactor.mean_square_error(y_true=y_true, y_pred=y_pred)

    print(f"R factor: {r_factor}")


if __name__ == "__main__":
    cli()
