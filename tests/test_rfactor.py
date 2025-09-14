import math

import matplotlib.pyplot as plt
import numpy as np
import pytest

from cleedpy import rfactor as rf


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
    l_curve = rf.lorentzian_smoothing(curve, vi)

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


def test_r2_factor():
    """Test for r2_factor with different curves."""
    x = np.arange(50, 300, 0.1)
    y = np.sin(x * np.pi / 20) ** 2 + 0.1  # The 0.1 here is to avoid zero values

    # Same curves should give 0.0
    assert math.isclose(
        rf.r2_factor(np.column_stack([x, y]), np.column_stack([x, y])),
        0.0,
        abs_tol=1e-8,
    )

    # Multiplying by a constant should not change the r2 factor
    y_times_constant = y * 5
    assert math.isclose(
        rf.r2_factor(np.column_stack([x, y]), np.column_stack([x, y_times_constant])),
        0,
        abs_tol=1e-8,
    )

    # Shifting x by 10 should put the maximum of one curve at the minimum of the other, so the r2 factor should be close to 2
    y_shifted = np.sin((x + 10) * np.pi / 20) ** 2 + 0.1
    assert math.isclose(
        rf.r2_factor(np.column_stack([x, y]), np.column_stack([x, y_shifted])),
        2,
        abs_tol=1e-2,
    )

    # A gradual increase of shift from 0 to 10 should lead to a smooth increase of the R factor.
    r_previous = 0.0
    for shift in np.arange(0.1, 10.1, 0.1):
        y_shifted = np.sin((x + shift) * np.pi / 20) ** 2 + 0.1
        r_current = rf.r2_factor(
            np.column_stack([x, y]), np.column_stack([x, y_shifted])
        )
        assert r_current > r_previous, f"R factor should increase with shift: {shift}"
        r_previous = r_current

    # Introducing modulation. A gradual increase of alpha from 0 to 1 should lead to a smooth increase of the R factor.
    r_previous = 0.0
    for alpha in np.arange(0.1, 1.1, 0.1):
        y_modulated = (
            np.sin(x * np.pi / 20) ** 2 * (1 + alpha * np.sin(x * np.pi / 100)) + 0.1
        )
        r_current = rf.r2_factor(
            np.column_stack([x, y]), np.column_stack([x, y_modulated])
        )
        assert (
            r_current > r_previous
        ), f"R factor should increase with modulation: {alpha}"
        r_previous = r_current


def test_rp_factor():
    """Test for rp_factor with different curves."""
    x = np.arange(50, 300, 0.1)
    y = np.sin(x * np.pi / 20) ** 2 + 0.1  # The 0.1 here is to avoid zero values

    # Same curves should give 0.0
    assert math.isclose(
        rf.rp_factor(np.column_stack([x, y]), np.column_stack([x, y])),
        0.0,
        abs_tol=1e-8,
    )

    # Multiplying by a constant should not change the rp factor
    y_times_constant = y * 5
    assert math.isclose(
        rf.rp_factor(np.column_stack([x, y]), np.column_stack([x, y_times_constant])),
        0,
        abs_tol=1e-8,
    )

    # Shifting x by 10 should put the maximum of one curve at the minimum of the other, so the rp factor should be close to 2
    y_shifted = np.sin((x + 10) * np.pi / 20) ** 2 + 0.1
    assert math.isclose(
        rf.rp_factor(np.column_stack([x, y]), np.column_stack([x, y_shifted])),
        2,
        abs_tol=1e-1,
    )

    # A gradual increase of shift from 0 to 10 should lead to a smooth increase of the R factor.
    r_previous = 0.0
    for shift in np.arange(0.1, 10.1, 0.1):
        y_shifted = np.sin((x + shift) * np.pi / 20) ** 2 + 0.1
        r_current = rf.rp_factor(
            np.column_stack([x, y]), np.column_stack([x, y_shifted])
        )
        assert r_current > r_previous, f"Rp factor should increase with shift: {shift}"
        r_previous = r_current

    # Introducing modulation. A gradual increase of alpha from 0 to 1 should keep the R factor close to 0.
    for alpha in np.arange(0.1, 1.1, 0.1):
        y_modulated = (
            np.sin(x * np.pi / 20) ** 2 * (1 + alpha * np.sin(x * np.pi / 100)) + 0.1
        )
        assert math.isclose(
            rf.rp_factor(np.column_stack([x, y]), np.column_stack([x, y_modulated])),
            0,
            abs_tol=1e-1,
        ), f"Rp factor should be close to 0 with modulation: {alpha}"


def test_find_common_x_axis():
    """Test for find_common_x_axis with different curves."""

    # Test with overlapping ranges
    x1 = np.array([1, 2, 3, 4, 5])
    x2 = np.array([3, 4, 5, 6, 7])
    expected = np.array([3, 4, 5])
    result = rf.find_common_x_axis(reference_grid=x1, other_grid=x2)
    assert np.array_equal(result, expected), f"Expected {expected}, but got {result}"

    # Test with no common points
    x1 = np.array([1, 2])
    x2 = np.array([3, 4])
    expected = np.array([])
    result = rf.find_common_x_axis(reference_grid=x1, other_grid=x2)
    assert np.array_equal(result, expected), f"Expected {expected}, but got {result}"

    # Test with overlapping ranges but different points
    x1 = np.array([1, 2, 3, 4, 5])
    x2 = np.array([2.4, 3.5, 4.6, 5.7, 6.8])
    expected = np.array([3, 4, 5])
    result = rf.find_common_x_axis(reference_grid=x1, other_grid=x2)
    assert np.array_equal(result, expected), f"Expected {expected}, but got {result}"
