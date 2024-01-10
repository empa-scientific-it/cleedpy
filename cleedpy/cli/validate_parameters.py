from pathlib import Path

import typer
import yaml

from ..config.bulk import NonGeometricalParameters


def validate(path: Path):
    """Validate bulk parameters"""
    with open(path) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
        print(data)
        NonGeometricalParameters.model_validate(data)
        obj = NonGeometricalParameters(**data)
        print(obj)


if __name__ == "__main__":
    typer.run(validate)
