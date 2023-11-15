import click
import yaml

@click.command("cli")
@click.option("--input", "-i", help="Input file")
def cli(input):
    """Command line interface for rfactor"""
    with open(input, "r") as f:
        data = yaml.safe_load(f)
    
    print(data)

if __name__ == "__main__":
    cli()