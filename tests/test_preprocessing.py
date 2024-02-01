import pytest
import math
import numpy as np
import matplotlib.pyplot as plt
from cleedpy.preprocessing import lorentzian_smoothing, preprocessing_loop
from tests.curves_helper import curve_A, curve_B, curve_A_smoothed


@pytest.mark.parametrize(
    "curve, vi, expected",
    [
        (np.array(curve_A()), 4, np.array(curve_A_smoothed())),
    ]
)
def test_lorentzian_smoothing(curve, vi, expected):
    L_curve = lorentzian_smoothing(curve, vi)
    
    x_values = curve[:, 0]
    y_values = curve[:, 1]
    plt.plot(x_values, y_values, marker='o', color='g', label='Input: Initial curve')

    x_values = L_curve[:, 0]
    y_values = L_curve[:, 1]
    plt.plot(x_values, y_values, marker='x', color='r', label='Output: Smoothed curve')

    y_min = 0
    y_max = 0.0110
    y_step = 0.0005
    plt.ylim(y_min, y_max)
    plt.yticks(np.arange(y_min, y_max + y_step, y_step))
    plt.grid()
    plt.legend()
    plt.savefig("test_lorentzian_smoothing_output.png")

    assert np.allclose(
        expected, L_curve
    )


@pytest.mark.parametrize(
    "the_curve, exp_curve, shift, r_factor, vi, expected",
    [
        ([curve_A(), curve_B()], [curve_A(), curve_B()], 1, "r2_factor", 4, 0),
    ]
)
def test_preprocessing_loop(the_curve, exp_curve, shift, r_factor, vi, expected):
    assert math.isclose(
        expected, preprocessing_loop(np.array(the_curve), np.array(exp_curve), shift, r_factor, vi), abs_tol=5
    )