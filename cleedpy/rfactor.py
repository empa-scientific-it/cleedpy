import numpy as np


def r2_factor(theoretical_curve, experimental_curve):
    """
    Calculates R2-factor of the curve of a specific spot (x,y) of the experiment.
    A curve is represented by an array of (e,i) points, where e = energy and i = intensity.

    Inputs:
        - theoretical curve: numpy array of [e,i] arrays
        - experimental curve: numpy array of [e,i] arrays

    Design:
        R2 = sqrt{ S(It - c*Ie)^2 / S(It - It_avg)^2 }

        where:
            c = sqrt( S|It|^2 / S|Ie|^ 2)
            It_avg = (S It)/ dE
    """

    # [TODO] (not sure) Use the length of the shortest curve
    min_length = min(len(theoretical_curve), len(experimental_curve))
    theoretical_curve = theoretical_curve[:min_length]
    experimental_curve = experimental_curve[:min_length]

    # Extract the It, Ie values for theoretical and experimental intensity
    It = theoretical_curve[:, 1]
    Ie = experimental_curve[:, 1]

    # Calculate normalization factor c and It_avg
    c = np.sqrt(np.sum(It**2) / np.sum(Ie**2))
    It_avg = np.sum(It) / It.size  # dE = number of energy steps

    # Calculate the numerator and denominator of R2
    numerator = np.sum((It - c * Ie) ** 2)
    denominator = np.sum((It - It_avg) ** 2)

    # Calculate R2
    R2 = np.sqrt(numerator / denominator)

    # [TODO] error handling: may return NaN
    return R2


def rp_factor(theoretical_curve, experimental_curve):
    """
    Calculates Pendry's R-factor of the curve of a specific spot (x,y) of the experiment.
    A curve is represented by an array of (e,i) points, where e = energy and i = intensity.

    Inputs:
        - theoretical curve: numpy array of [e,i] arrays
        - experimental curve: numpy array of [e,i] arrays

    Design:
        Rp = S(Ye - Yt)^2 / S(Ye^2 + Yt^2)
    """

    # Extract the It, Ie values for theoretical and experimental intensity
    It = theoretical_curve[:, 1]
    Ie = experimental_curve[:, 1]

    # Extract the Et, Ee values for theoretical and experimental energy
    Et = theoretical_curve[:, 0]
    Ee = experimental_curve[:, 0]

    # Calculate theoretical and experimental energy steps
    step_Et = energy_step(Et)
    step_Ee = energy_step(Ee)

    # Calculate Y for theoretical and experimental intensity
    Yt = Y_function(It, step_Et)
    Ye = Y_function(Ie, step_Ee)

    # Calculate the numerator and denominator of Rp
    numerator = np.sum((Yt - Ye) ** 2 * step_Et)
    denominator = np.sum(Ye**2 * step_Ee) + np.sum(Yt**2 * step_Et)

    # Calculate Rp
    Rp = numerator / denominator

    # [TODO] error handling: may return NaN
    return Rp


def Y_function(I, energy_step):
    """
    Calculates Y for a given curve.

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

    L = (I[1:] - I[:-1]) / (energy_step * 0.5 * (I[1:] + I[:-1]))
    Y = L / (1 + L**2 * vi**2)

    return Y


def energy_step(E):
    """
    Calculates the pairwise energy step. Returns an array.

    Inputs:
        - E: array of energy values of the curve

    Design:
        step = E[i] - E[i-1]
    """

    return E[1:] - E[:-1]
