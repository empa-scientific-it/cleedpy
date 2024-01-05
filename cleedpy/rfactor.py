import numpy as np

def mean_square_error(y_true, y_pred):
    """Mean Square Error"""
    return np.mean(np.square(y_true - y_pred))
    