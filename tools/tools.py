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

    return phi.item(), theta.item(), psi.item()


def Euler2Quaternion(phi, theta, psi):
    '''
    Conversion from euler angles to quaternions, in page 259
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


def Euler2Rotation(phi, theta, psi):
    """
    Converts euler angles to rotation matrix (R_b^i, i.e., body to inertial)
    """
    # only call sin and cos once for each angle to speed up rendering
    c_phi = np.cos(phi)
    s_phi = np.sin(phi)
    c_theta = np.cos(theta)
    s_theta = np.sin(theta)
    c_psi = np.cos(psi)
    s_psi = np.sin(psi)

    R_roll = np.array([[1, 0, 0],
                       [0, c_phi, s_phi],
                       [0, -s_phi, c_phi]], dtype=object)
    R_pitch = np.array([[c_theta, 0, -s_theta],
                        [0, 1, 0],
                        [s_theta, 0, c_theta]], dtype=object)
    R_yaw = np.array([[c_psi, s_psi, 0],
                      [-s_psi, c_psi, 0],
                      [0, 0, 1]], dtype=object)
    R = R_roll @ R_pitch @ R_yaw  # inertial to body (Equation 2.4 in book)
    return R.T  # transpose to return body to inertial


if __name__ == "__main__":
    phi = 120 * np.pi / 180
    theta = 10 * np.pi / 180
    psi = -11 * np.pi / 180
    e = Euler2Quaternion(phi, theta, psi)
    print("The quaternions is ", e)

    phi_new, theta_new, psi_new = Quaternion2Euler(e)
    print("The euler angle is ", phi_new * 180 / np.pi, theta_new * 180 / np.pi, psi_new * 180 / np.pi)

    R = Euler2Rotation(phi, theta, psi)
    print("The rotation matrix is ", R)

    # After test, there is small quantization error in the conversion.
