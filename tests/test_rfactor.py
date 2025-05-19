import math

import numpy as np

from cleedpy.rfactor import r2_factor, rp_factor


def test_r2_factor_other():
    """Test for r2_factor with different curves."""
    x = np.arange(50, 300, 0.1)
    y = np.sin(x * np.pi / 20) ** 2 + 0.1  # The 0.1 here is to avoid zero values

    # Same curves should give 0.0
    assert math.isclose(
        r2_factor(np.column_stack([x, y]), np.column_stack([x, y])), 0.0, abs_tol=1e-8
    )

    # Multiplying by a constant should not change the r2 factor
    y_times_constant = y * 5
    assert math.isclose(
        r2_factor(np.column_stack([x, y]), np.column_stack([x, y_times_constant])),
        0,
        abs_tol=1e-8,
    )

    # Shifting x by 10 should put the maximum of one curve at the minimum of the other, so the r2 factor should be close to 2
    y_shifted = np.sin((x + 10) * np.pi / 20) ** 2 + 0.1
    assert math.isclose(
        r2_factor(np.column_stack([x, y]), np.column_stack([x, y_shifted])),
        2,
        abs_tol=1e-2,
    )

    # A gradual increase of shift from 0 to 10 should lead to a smooth increase of the R factor.
    r_previous = 0.0
    for shift in np.arange(0.1, 10.1, 0.1):
        y_shifted = np.sin((x + shift) * np.pi / 20) ** 2 + 0.1
        r_current = r2_factor(np.column_stack([x, y]), np.column_stack([x, y_shifted]))
        assert r_current > r_previous, f"R factor should increase with shift: {shift}"
        r_previous = r_current

    # Introducing modulation. A gradual increase of alpha from 0 to 1 should lead to a smooth increase of the R factor.
    r_previous = 0.0
    for alpha in np.arange(0.1, 1.1, 0.1):
        y_modulated = (
            np.sin(x * np.pi / 20) ** 2 * (1 + alpha * np.sin(x * np.pi / 100)) + 0.1
        )
        r_current = r2_factor(
            np.column_stack([x, y]), np.column_stack([x, y_modulated])
        )
        assert (
            r_current > r_previous
        ), f"R factor should increase with modulation: {alpha}"
        r_previous = r_current


def test_rp_factor_other():
    """Test for rp_factor with different curves."""
    x = np.arange(50, 300, 0.1)
    y = np.sin(x * np.pi / 20) ** 2 + 0.1  # The 0.1 here is to avoid zero values

    # Same curves should give 0.0
    assert math.isclose(
        rp_factor(np.column_stack([x, y]), np.column_stack([x, y])), 0.0, abs_tol=1e-8
    )

    # Multiplying by a constant should not change the rp factor
    y_times_constant = y * 5
    assert math.isclose(
        rp_factor(np.column_stack([x, y]), np.column_stack([x, y_times_constant])),
        0,
        abs_tol=1e-8,
    )

    # Shifting x by 10 should put the maximum of one curve at the minimum of the other, so the rp factor should be close to 2
    y_shifted = np.sin((x + 10) * np.pi / 20) ** 2 + 0.1
    assert math.isclose(
        rp_factor(np.column_stack([x, y]), np.column_stack([x, y_shifted])),
        2,
        abs_tol=1e-1,
    )

    # A gradual increase of shift from 0 to 10 should lead to a smooth increase of the R factor.
    r_previous = 0.0
    for shift in np.arange(0.1, 10.1, 0.1):
        y_shifted = np.sin((x + shift) * np.pi / 20) ** 2 + 0.1
        r_current = rp_factor(np.column_stack([x, y]), np.column_stack([x, y_shifted]))
        assert r_current > r_previous, f"Rp factor should increase with shift: {shift}"
        r_previous = r_current

    # Introducing modulation. A gradual increase of alpha from 0 to 1 should keep the R factor close to 0.
    for alpha in np.arange(0.1, 1.1, 0.1):
        y_modulated = (
            np.sin(x * np.pi / 20) ** 2 * (1 + alpha * np.sin(x * np.pi / 100)) + 0.1
        )
        assert math.isclose(
            rp_factor(np.column_stack([x, y]), np.column_stack([x, y_modulated])),
            0,
            abs_tol=1e-1,
        ), f"Rp factor should be close to 0 with modulation: {alpha}"
