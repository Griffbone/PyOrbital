import numpy as np


def rotmat(t, axis):
    """ Function to create a rotation matrix for a rotation about a given axis
        :param t: rotation angle theta (rad)
        :param axis: axis of rotation ('X', 'Y', or 'Z')

        :return mat: rotation matrix
    """

    if axis.upper() == 'X':
        return np.array([
            [1, 0, 0],
            [0, np.cos(t), -np.sin(t)],
            [0, np.sin(t), np.cos(t)]
        ])

    elif axis.upper() == 'Y':
        return np.array([
            [np.cos(t), 0, np.sin(t)],
            [0, 1, 0],
            [-np.sin(t), 0, np.cos(t)]
        ])

    elif axis.upper() == 'Z':
        return np.array([
            [np.cos(t), -np.sin(t), 0],
            [np.sin(t), np.cos(t), 0],
            [0, 0, 1]
        ])
