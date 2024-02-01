import pytest
import math
import numpy as np
from cleedpy.rfactor import r2_factor, rp_factor
from tests.curves_helper import curve_A, curve_B, curve_C


@pytest.mark.parametrize(
    "the_curve, exp_curve, expected_r",
    [
        (curve_A(), curve_B(), 1.8688830931),
        (curve_B(), curve_C(), 0.0001429368),
    ]
)
def test_r2_factor(the_curve, exp_curve, expected_r):
    assert math.isclose(
        expected_r, r2_factor(np.array(the_curve), np.array(exp_curve)), abs_tol=5
    )


@pytest.mark.parametrize(
    "the_curve, exp_curve, expected_r",
    [
        (curve_A(), curve_B(), 1.4898282448),
        (curve_B(), curve_C(), 0),
    ]
)
def test_rp_factor(the_curve, exp_curve, expected_r):
    assert math.isclose(
        expected_r, rp_factor(np.array(the_curve), np.array(exp_curve)), abs_tol=5
    )
