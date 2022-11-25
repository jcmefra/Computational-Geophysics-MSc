import numpy as np

def vec2mat(a, new_shape):

    padding = (new_shape - np.vstack((new_shape, a.shape))).T.tolist()

    return np.pad(a, padding, mode='constant')