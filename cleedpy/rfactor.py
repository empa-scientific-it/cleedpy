import numpy as np
from scipy.interpolate import CubicSpline

from . import config

GEO_PREFACTOR = (
    5.0  # To scale the geometrical rfactor to be comparable with the iv rfactor.
)


def lorentzian_smoothing(curve, vi=4.0):
    """
    Draws a Lorentzian curve in the same grid as the given curve.
    Performs the convolution of the two curves.
    The given curve is represented by an array of (e,i) points, where e = energy and i = intensity.

    Inputs:
        - curve: numpy array of [e,i] arrays
        - vi: predefined constant representing the width

    Output:
        The result of the convolution: a new array of [e,i] arrays
    """

    energies = curve[:, 0]
    intensities = curve[:, 1]
    result = []

    for j in range(energies.size):
        c = 0
        l_value = 0
        for i in range(energies.size):
            c += vi**2 / ((energies[i] - energies[j]) ** 2 + vi**2)
            l_value += (
                intensities[i] * vi**2 / ((energies[i] - energies[j]) ** 2 + vi**2)
            )
        result.append(l_value / c)

    return np.array(list(zip(energies, result)))


def r2_factor(experimental_curve, theoretical_curve):
    """
    Calculates R2-factor of the curve of a specific spot (x,y) of the experiment.
    A curve is represented by an array of (e,i) points, where e = energy and i = intensity.

    Inputs:
        - theoretical curve: numpy array of [e,i] arrays
        - experimental curve: numpy array of [e,i] arrays

    Design:
        R2 = sqrt{ S(it - c*ie)^2 / S(it - it_avg)^2 }

        where:
            c = sqrt( S|it|^2 / S|ie|^ 2)
            it_avg = (S it)/ dE
    """

    # Extract the it, ie values for theoretical and experimental intensity
    it = theoretical_curve[:, 1]
    ie = experimental_curve[:, 1]

    # Calculate normalization factor c and it_avg
    c = np.sqrt(np.sum(it**2) / np.sum(ie**2))
    it_avg = np.sum(it) / it.size

    # Calculate the numerator and denominator of R2
    numerator = np.sum((it - c * ie) ** 2)
    denominator = np.sum((it - it_avg) ** 2)

    # Calculate R2
    r2 = np.sqrt(numerator / denominator)

    # [TODO] error handling: may return NaN
    return r2


def rp_factor(experimental_curve, theoretical_curve):
    """
    Calculates Pendry's R-factor of the curve of a specific spot (x,y) of the experiment.
    A curve is represented by an array of (e,i) points, where e = energy and i = intensity.

    Inputs:
        - theoretical curve: numpy array of [e,i] arrays
        - experimental curve: numpy array of [e,i] arrays

    Design:
        Rp = S(ye - yt)^2 / S(ye^2 + yt^2)
    """

    # Extract the it, ie values for theoretical and experimental intensity
    it = theoretical_curve[:, 1]
    ie = experimental_curve[:, 1]

    # Extract the et, ee values for theoretical and experimental energy
    et = theoretical_curve[:, 0]
    ee = experimental_curve[:, 0]

    # Calculate theoretical and experimental energy steps
    step_et = energy_step(et)
    step_ee = energy_step(ee)

    # Calculate Y for theoretical and experimental intensity
    yt = y_function(it, step_et)
    ye = y_function(ie, step_ee)

    # Calculate the numerator and denominator of Rp
    numerator = np.sum((yt - ye) ** 2 * step_et)
    denominator = np.sum(ye**2 * step_ee) + np.sum(yt**2 * step_et)

    # Calculate Rp
    rp = numerator / denominator

    # [TODO] error handling: may return NaN
    return rp


def y_function(intensities, energy_step):
    """
    Calculates y for a given curve.

    Inputs:
        - I: array of intensity values of the curve
        - E: array of energy values of the curve

    Design:
        Y = L / (1 + L^2 * vi^2)

        where:
            L = (I[i] - I[i-1]) / (energy_step * 0.5 * (I[i] + I[i-1]))
    """

    # [TODO] constants should be defined elsewhere
    vi = 4

    y_values = (intensities[1:] - intensities[:-1]) / (
        energy_step * 0.5 * (intensities[1:] + intensities[:-1])
    )
    return y_values / (1 + y_values**2 * vi**2)


def energy_step(energies):
    """
    Calculates the pairwise energy step. Returns an array.

    Inputs:
        - E: array of energy values of the curve

    Design:
        step = E[i] - E[i-1]
    """

    return energies[1:] - energies[:-1]


def split_in_pairs(exp_iv, theo_iv):
    """This function takes a complete datasets with a experimental and theoretical iv curves splits them them in pairs corresponding to the same index."""

    beam_indices = np.unique(np.vstack([exp_iv, theo_iv])[:, [0, 1]], axis=0)

    for indx1, indx2 in beam_indices:
        exp_curve = exp_iv[(exp_iv[:, 0] == indx1) & (exp_iv[:, 1] == indx2)][
            :, [-2, -1]
        ]
        theo_curve = theo_iv[(theo_iv[:, 0] == indx1) & (theo_iv[:, 1] == indx2)][
            :, [-2, -1]
        ]
        yield exp_curve, theo_curve


def find_common_x_axis(reference_grid, other_grid):
    min_x = max(np.min(reference_grid), np.min(other_grid))
    max_x = min(np.max(reference_grid), np.max(other_grid))

    reference_grid = reference_grid[
        (reference_grid >= min_x) & (reference_grid <= max_x)
    ]
    return reference_grid


def compute_rfactor(experimental_iv, theoretical_iv, shift=0.0, rfactor_type="r2"):

    r_tot = 0.0
    rfactor = {
        "r2": r2_factor,
        "pendry": rp_factor,
    }

    delta_e = []  # To store the energy range of each curve.

    # Looping over the pairs of experimental and theoretical curves corresponding to the same index.
    for exp_curve, theo_curve in split_in_pairs(experimental_iv, theoretical_iv):
        if len(exp_curve) < 1 or len(theo_curve) < 1:
            continue

        # Perform Lorentzian smoothing.
        exp_curve = lorentzian_smoothing(exp_curve)
        theo_curve = lorentzian_smoothing(theo_curve)

        # Shifting the theoretical curve.
        theo_curve[:, 0] += shift

        # Finding common x axis.
        common_x = find_common_x_axis(exp_curve[:, 0], theo_curve[:, 0])

        # If no common points - skip the curve.
        if len(common_x) == 0:
            continue

        exp_spline = CubicSpline(x=exp_curve[:, 0], y=exp_curve[:, 1])(common_x)
        theo_spline = CubicSpline(x=theo_curve[:, 0], y=theo_curve[:, 1])(common_x)

        # Uncomment for testing.
        # import matplotlib.pyplot as plt  # noqa: E800
        # plt.plot(theo_curve[:,0], theo_curve[:,1], 'o', label='Data Points')  # noqa: E800
        # plt.plot(common_x, theo_spline, label='Cubic Spline')  # noqa: E800
        # #plt.plot(exp_curve_before_smooth[:,0], exp_curve_before_smooth[:,1], label='Before smooth')  # noqa: E800
        # plt.legend()  # noqa: E800
        # plt.show()  # noqa: E800

        r = rfactor[rfactor_type](
            np.column_stack([common_x, exp_spline]),
            np.column_stack([common_x, theo_spline]),
        )

        delta_e.append(common_x[-1] - common_x[0])
        r_tot += r * (common_x[-1] - common_x[0])

    return r_tot / sum(delta_e)


def compute_geometrical_rfactor(config: config.InputParameters):
    """Compute R-factor for a given geometry configuration"""

    structure = config.get_ase_structure()
    distance_matrix = structure.get_all_distances(mic=True)

    # Create a matrix of minimum allowed distances based on atomic radii.
    atomic_radii = np.array(
        [config.minimum_radius[atom] for atom in structure.get_chemical_symbols()]
    )
    min_distance_matrix = atomic_radii[:, None] + atomic_radii[None, :]

    # Prepare a mask to set the diagonal and lower triangular to False
    # to avoid self-comparison and double counting.
    mask = np.triu(np.ones(distance_matrix.shape, dtype=bool), k=1)

    # Check for overlaps
    overlap_matrix = (distance_matrix < min_distance_matrix) & mask

    # Sum distances that overlap.
    total_overlap = (
        (min_distance_matrix - distance_matrix) * overlap_matrix * GEO_PREFACTOR
    )
    return np.sum(total_overlap**2)
