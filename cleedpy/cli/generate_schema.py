import json
from pathlib import Path

import pydantic
import typer

from ..config.bulk import NonGeometricalParameters


def generate_schema(path: Path):
    """Generate schema for bulk calculations"""
    schema = NonGeometricalParameters.model_json_schema()

    with open(path, "w") as f:
        f.write(json.dumps(schema, indent=4))


if __name__ == "__main__":
    typer.run(generate_schema)
