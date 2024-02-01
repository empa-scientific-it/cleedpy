import numpy as np 
from cleedpy.rfactor import r2_factor, rp_factor
from scipy.interpolate import CubicHermiteSpline


def lorentzian_smoothing(curve, vi):
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

    E = curve[:, 0]
    I = curve[:, 1]
    L = []

    for j in range(E.size):
        c = 0
        L_value = 0
        for i in range(E.size):
            c += vi**2 / ((E[i] - E[j])**2 + vi**2)
            L_value += (I[i] * vi**2 / ((E[i] - E[j])**2 + vi**2))
        L.append(L_value / c)

    return np.array(list(zip(E, L)))


def preprocessing_loop(theoretical_curves, experimental_curves, shift, r_factor, vi):
    """
        Performs all the preprocessing steps and R-factor calculations, required before the Search.
        A curve is represented by an array of (e,i) points, where e = energy and i = intensity.

        Inputs:
            - theoretical_curves: numpy array of arrays of [e,i] arrays
            - experimental_curves: numpy array of arrays of [e,i] arrays
                where:
                    - number of theoretical curves = number of experimental curves
            - shift: provided by user, used to edit the energy value per step
            - r_factor: provided by user, name of the preferred R-factor method

        Design:
            1. Lorentzian smoothing of all curves
            2. For each energy point:
                a. Shift theoretical curve
                b. Find boundaries of overlap
                c. Filter both curves to keep only the elements in the overlap
                d. For each experimental curve:
                    i. Spline = interpolate it to the theoretical grid
                    ii. Calculate R-factor
                e. Calculate average R-factor of all curves in this shift
            3. Find the min of the average values

        Output:
            The optimal R-factor value per curve.
    """

    for i in range(len(theoretical_curves)):
        theoretical_curves[i] = lorentzian_smoothing(theoretical_curves[i], vi)
        experimental_curves[i] = lorentzian_smoothing(experimental_curves[i], vi)

    energy_shift = [shift, 0]
    shifted_theoretical_curves = theoretical_curves + energy_shift
    print(shifted_theoretical_curves)

    # Find the boundaries of the overlap after the shift
    min_bound = np.maximum(shifted_theoretical_curves[0, 0], experimental_curves[0, 0])[0]
    max_bound = np.minimum(shifted_theoretical_curves[0, -1], experimental_curves[0, -1])[0]
    print(min_bound, max_bound)

    # Filter both curves based on the boundaries
    filtered_theoretical_curves = shifted_theoretical_curves[shifted_theoretical_curves[:, 0] >= min_bound]
    filtered_theoretical_curves = filtered_theoretical_curves[filtered_theoretical_curves[:, 0] <= max_bound]

    filtered_experimental_curves = experimental_curves[experimental_curves[:, 0] >= min_bound]
    filtered_experimental_curves = filtered_experimental_curves[filtered_experimental_curves[:, 0] <= max_bound]
    print(filtered_theoretical_curves, filtered_experimental_curves)

    # Spline
    for i in range(len(filtered_experimental_curves)):
        E_values = filtered_experimental_curves[i][:, 0]
        I_values = filtered_experimental_curves[i][:, 1]
        theoretical_E_grid = filtered_theoretical_curves[i][:, 0]
        CubicHermiteSpline(E_values, I_values, theoretical_E_grid)


    print(globals().get(r_factor))
    return
