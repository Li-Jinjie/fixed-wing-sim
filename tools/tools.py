#!/usr/bin/env python
'''
@Author : Li-Jinjie
@Date    : 2020/12/30
'''
import numpy as np


def Quaternion2Euler(quaternion):
    '''
    conversion from quaternions to euler angles
    Args:
        quaternion: ndarray, four numbers of quaternions: e0, e1, e2, e3

    Returns:
        euler angles: phi, theta, psi
    '''
    # phi, theta, psi = Quaternion2Euler(self._state[6:10])
    e0, e1, e2, e3 = quaternion[0], quaternion[1], quaternion[2], quaternion[3]
    phi = np.arctan2(2 * (e0 * e1 + e2 * e3), (e0 ** 2 + e3 ** 2 - e1 ** 2 - e2 ** 2))
    theta = np.arcsin(2 * (e0 * e2 - e1 * e3))
    psi = np.arctan2(2 * (e0 * e3 + e1 * e2), (e0 ** 2 + e1 ** 2 - e2 ** 2 - e3 ** 2))

    return phi, theta, psi
