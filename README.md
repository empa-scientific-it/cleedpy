# cleedpy

Python port of the CLEED code written by [Georg Held](https://github.com/GeorgHeld).

## Installation

The easiest way to install the package is via pip:

```bash
pip install cleedpy
```

## Usage

The cleedpy package provides a command line interface (CLI) to run LEED calculations.
Those include: `rfactor`, `search` and `leed` sub-programs.
Each program can be called with the `cleedpy-` prefix, e.g. `cleedpy-leed`:

```bash
cleedpy-leed -i input.yml -e experiment.txt -o search.out -p PHASE
```
To learn more about the options of each program, use the `-h` flag:

```bash
cleedpy-leed --help
```

For example runs please see the [examples folder](https://github.com/empa-scientific-it/cleedpy/tree/main/examples).

## Documentation

The documentation is available at the [Wiki page](https://github.com/empa-scientific-it/cleedpy/wiki) of the repository.


## License
MIT
