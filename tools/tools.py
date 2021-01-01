#!/usr/bin/env python
'''
@Author : Li-Jinjie
@Date    : 2020/12/30
'''
import numpy as np


def Quaternion2Euler(quaternion):
    '''
    Conversion from quaternions to euler angles, page 259
    Args:
        quaternion: ndarray, four numbers of quaternions: e0, e1, e2, e3

    Returns:
        euler angles: phi, theta, psi
    '''
    e0, e1, e2, e3 = quaternion[0], quaternion[1], quaternion[2], quaternion[3]
    phi = np.arctan2(2 * (e0 * e1 + e2 * e3), (e0 ** 2 + e3 ** 2 - e1 ** 2 - e2 ** 2))
    theta = np.arcsin(2 * (e0 * e2 - e1 * e3))
    psi = np.arctan2(2 * (e0 * e3 + e1 * e2), (e0 ** 2 + e1 ** 2 - e2 ** 2 - e3 ** 2))

    return phi, theta, psi


def Euler2Quaternion(phi, theta, psi):
    '''
    Conversion from euler angles to quaternions, page 259
    Args:
        Euler angles: phi, theta, psi

    Returns:
        quaternions e
    '''
    s_phi_2 = np.sin(phi / 2)
    c_phi_2 = np.cos(phi / 2)
    s_theta_2 = np.sin(theta / 2)
    c_theta_2 = np.cos(theta / 2)
    s_psi_2 = np.sin(psi / 2)
    c_psi_2 = np.cos(psi / 2)

    e = np.zeros(4)
    e[0] = c_psi_2 * c_theta_2 * c_phi_2 + s_psi_2 * s_theta_2 * s_phi_2
    e[1] = c_psi_2 * c_theta_2 * s_phi_2 - s_psi_2 * s_theta_2 * c_phi_2
    e[2] = c_psi_2 * s_theta_2 * c_phi_2 + s_psi_2 * c_theta_2 * s_phi_2
    e[3] = s_psi_2 * c_theta_2 * c_phi_2 - c_psi_2 * s_theta_2 * s_phi_2

    return e


if __name__ == "__main__":
    phi = 120 * np.pi / 180
    theta = 10 * np.pi / 180
    psi = -11 * np.pi / 180
    e = Euler2Quaternion(phi, theta, psi)
    print("quaternions is ", e)

    phi_new, theta_new, psi_new = Quaternion2Euler(e)
    print("the Euler angle is ", phi_new * 180 / np.pi, theta_new * 180 / np.pi, psi_new * 180 / np.pi)

    # After test, there is small quantization error in the conversion.
