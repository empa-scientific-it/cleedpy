import matplotlib.pyplot as plt
import numpy as np
import pytest

from cleedpy.preprocessing import lorentzian_smoothing


def curve_a():
    return [
        [70, 0.0037557780],
        [74, 0.0013476560],
        [78, 0.0010809600],
        [82, 0.0009844219],
        [86, 0.0009538341],
        [90, 0.0011735780],
        [94, 0.0020866750],
        [98, 0.0040076710],
        [102, 0.0049233240],
        [106, 0.0051959950],
        [110, 0.0044225620],
        [114, 0.0020602790],
        [118, 0.0005740963],
        [122, 0.0005428890],
        [126, 0.0002391057],
    ]


def curve_a_smoothed():
    return [
        [70, 0.00258495],
        [74, 0.00182937],
        [78, 0.00145542],
        [82, 0.00135855],
        [86, 0.00144012],
        [90, 0.00176726],
        [94, 0.00244699],
        [98, 0.00335730],
        [102, 0.00395957],
        [106, 0.00404653],
        [110, 0.00348618],
        [114, 0.00238778],
        [118, 0.00145559],
        [122, 0.00101030],
        [126, 0.00079836],
    ]


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
