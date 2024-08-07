import json
from pathlib import Path

import typer

from ..config import InputParameters


def generate_schema(path: Path):
    """Generate schema for the LEED calculations"""
    schema = InputParameters.model_json_schema()

    with open(path, "w") as f:
        f.write(json.dumps(schema, indent=4))


if __name__ == "__main__":
    typer.run(generate_schema)
