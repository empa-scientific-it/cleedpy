import numpy as np
import pytest

from cleedpy import search


@pytest.mark.parametrize(
    "f, x0, expected",
    [
        (lambda x: x[0] ** 2, [2], [0]),
        (lambda x: (x[0] - 2) ** 2 + (x[1] - 3) ** 2, [0, 0], [2, 3]),
    ],
)
def test_simplex(f, x0, expected):
    res = search.simplex(f, x0)
    assert np.isclose(res, np.array(expected), atol=1e-5).all()
