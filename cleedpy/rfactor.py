import numpy as np


def r2_factor(theoretical_curve, experimental_curve):
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

    # [TODO] (not sure) Use the length of the shortest curve
    min_length = min(len(theoretical_curve), len(experimental_curve))
    theoretical_curve = theoretical_curve[:min_length]
    experimental_curve = experimental_curve[:min_length]

    # Extract the it, ie values for theoretical and experimental intensity
    it = theoretical_curve[:, 1]
    ie = experimental_curve[:, 1]

    # Calculate normalization factor c and it_avg
    c = np.sqrt(np.sum(it**2) / np.sum(ie**2))
    it_avg = np.sum(it) / it.size  # dE = number of energy steps

    # Calculate the numerator and denominator of R2
    numerator = np.sum((it - c * ie) ** 2)
    denominator = np.sum((it - it_avg) ** 2)

    # Calculate R2
    r2 = np.sqrt(numerator / denominator)

    # [TODO] error handling: may return NaN
    return r2


def rp_factor(theoretical_curve, experimental_curve):
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
