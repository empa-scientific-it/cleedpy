import math

import numpy as np
import pytest

from cleedpy import rfactor


@pytest.mark.parametrize(
    "y_true, y_pred, expected",
    [
        ([1, 2, 3], [1, 2, 3], 0),
        ([1, 2, 3], [1.1, 2.2, 3.3], 0.04666666666666666),
        ([1.0, 2.0, 3.0], [1.0, 2.0, 5.0], 1.3333333333333333),
    ],
)
def test_mean_squer_error(y_true, y_pred, expected):
    assert math.isclose(
        expected, rfactor.mean_square_error(np.array(y_true), np.array(y_pred))
    )
