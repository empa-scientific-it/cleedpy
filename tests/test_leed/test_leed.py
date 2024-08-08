from pathlib import Path

import numpy as np
import pytest

from cleedpy.interface.cleed import call_cleed
from cleedpy.physics.constants import HART


@pytest.mark.parametrize(
    "folder",
    [
        "../../examples/cleedpy_example_intermediate/",
    ],
)
def test_leed(folder):
    # Script directory with pathlib

    script_dir = Path(__file__).resolve().parent
    parameter_file = script_dir / folder / "leed.inp"
    phase_shift = script_dir / folder / "PHASE"
    result = call_cleed(str(parameter_file), str(parameter_file), str(phase_shift))

    assert result.n_beams == 55
    assert result.n_energies == 108

    # Read beams.txt file using numpy. The file contains 3 colums: 1st beam index (float), 2nd beam index (float), and beam set (int)
    beams = np.loadtxt(script_dir / folder / "beams.txt", dtype=float)

    # Compare the beam indexes from the file with the result.beam_index
    assert np.allclose(
        [result.beam_index1[i] for i in range(result.n_beams)], beams[:, 0]
    )
    assert np.allclose(
        [result.beam_index2[i] for i in range(result.n_beams)], beams[:, 1]
    )
    assert np.allclose([result.beam_set[i] for i in range(result.n_beams)], beams[:, 2])

    # Read iv_curves.txt file using numpy. The file contains 56 columns, the first column is the energy and the rest are the iv curves
    iv_curves = np.loadtxt(script_dir / folder / "iv_curves.txt", dtype=float).reshape(
        result.n_energies, result.n_beams + 1
    )
    assert np.allclose(
        iv_curves[:, 0], [result.energies[i] * HART for i in range(result.n_energies)]
    )
    iv_curves = iv_curves[:, 1:]
    for i in range(result.n_energies):
        assert np.allclose(
            iv_curves[i],
            [result.iv_curves[i * result.n_beams + j] for j in range(result.n_beams)],
        )
