import matplotlib.pyplot as plt
import numpy as np
import pytest

from cleedpy.preprocessing import lorentzian_smoothing
from tests.curves_helper import curve_a, curve_a_smoothed


@pytest.mark.parametrize(
    "curve, vi, expected",
    [
        (np.array(curve_a()), 4, np.array(curve_a_smoothed())),
    ],
)
def test_lorentzian_smoothing(curve, vi, expected):
    l_curve = lorentzian_smoothing(curve, vi)

    x_values = curve[:, 0]
    y_values = curve[:, 1]
    plt.plot(x_values, y_values, marker="o", color="g", label="Input: Initial curve")

    x_values = l_curve[:, 0]
    y_values = l_curve[:, 1]
    plt.plot(x_values, y_values, marker="x", color="r", label="Output: Smoothed curve")

    y_min = 0
    y_max = 0.0110
    y_step = 0.0005
    plt.ylim(y_min, y_max)
    plt.yticks(np.arange(y_min, y_max + y_step, y_step))
    plt.grid()
    plt.legend()
    plt.savefig("test_lorentzian_smoothing_output.png")

    assert np.allclose(expected, l_curve)
