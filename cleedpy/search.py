from scipy import optimize


def simplex(f, x0):
    """Minimize a function using the Nelder-Mead simplex algorithm.

    :param f: Function to minimize.
    :param x0: Initial guess.
    :return: Optimal value.
    """
    res = optimize.minimize(f, x0, method="Nelder-Mead", tol=1e-7)
    return res.x
