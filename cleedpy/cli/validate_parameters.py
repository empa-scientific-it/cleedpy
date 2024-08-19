from pathlib import Path

import typer
import yaml

from ..config import InputParameters


def validate(path: Path):
    """Validate the input parameters"""
    with open(path) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
        print(data)
        InputParameters.model_validate(data)
        obj = InputParameters(**data)
        print(obj)


if __name__ == "__main__":
    typer.run(validate)
