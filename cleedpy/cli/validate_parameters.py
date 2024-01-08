import typer
import pydantic
from pathlib import Path 
import json
import yaml

from cleedpy.config.bulk import NonGeometricalParameters

def validate(path: Path):
    """Validate bulk parameters"""
    with open(path, "r") as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
        print(data)
        NonGeometricalParameters.model_validate(data)
        obj = NonGeometricalParameters(**data)
        print(obj.bulk_layers[0].vibrational_displacement)

    


if __name__ == "__main__":
    typer.run(validate)