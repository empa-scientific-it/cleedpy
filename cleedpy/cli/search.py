import click
import numpy as np
import yaml

from .. import search


@click.command("cli")
@click.option("--config", "-c", help="Config file", required=True)
def cli(config):
    """Command line interface for the search tool"""
    with open(config) as fobj:
        data = yaml.safe_load(fobj)

    x_init = np.array(data["x_init"])

    def f(x):
        return (x[0] - 2) ** 2 + (x[1] - 3) ** 2

    res = search.simplex(f, x_init)
    print(f"Solution: {res}")


if __name__ == "__main__":
    cli()
