"""
compute_ss_model
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        2/4/2019 - RWB
"""
import sys

sys.path.append('..')
import numpy as np
from scipy.optimize import minimize
from tools.tools import Euler2Quaternion, Quaternion2Euler
from tools.transfer_function import transfer_function
import parameters.aerosonde_parameters as MAV
from parameters.simulation_parameters import ts_simulation as Ts


def compute_tf_model(mav, trim_state, trim_input):
    # trim values

    return T_phi_delta_a, T_chi_phi, T_theta_delta_e, T_h_theta, T_h_Va, T_Va_delta_t, T_Va_theta, T_beta_delta_r


def compute_ss_model(mav, trim_state, trim_input):
    # x_euler = [pn, pe, pd, u, v, w, phi, theta, psi, p, q, r].T
    x_euler = euler_state(trim_state)
    pn = x_euler.item(0)
    pe = x_euler.item(1)
    pd = x_euler.item(2)
    u = x_euler.item(3)
    v = x_euler.item(4)
    w = x_euler.item(5)
    phi = x_euler.item(6)
    theta = x_euler.item(7)
    psi = x_euler.item(8)
    p = x_euler.item(9)
    q = x_euler.item(10)
    r = x_euler.item(11)

    delta_a = trim_input.item(0)
    delta_e = trim_input.item(1)
    delta_r = trim_input.item(2)
    delta_t = trim_input.item(3)

    Va = mav._Va
    alpha = mav._alpha
    beta = mav._beta

    # m: External moment applied to the airframe about the body frame y-axis.
    rVS_2 = MAV.rho * (Va ** 2) * MAV.S_wing / 2
    m = rVS_2 * MAV.c * (MAV.C_m_0 + MAV.C_m_alpha * alpha +
                         MAV.C_m_q * MAV.c / (2 * Va) * q + MAV.C_m_delta_e * delta_e)

    # Lateral Model Coefficients
    Y_v = MAV.rho * MAV.S_wing * MAV.b * v / (4 * m * Va) * (MAV.C_Y_p * p + MAV.C_Y_r * r) \
          + (MAV.rho * MAV.S_wing * v / MAV.mass) * (
                  MAV.C_Y_0 + MAV.C_Y_beta * beta + MAV.C_Y_delta_a * delta_a + MAV.C_Y_delta_r * delta_r) \
          + (MAV.rho * MAV.S_wing * MAV.C_Y_beta / (2 * MAV.mass)) * np.sqrt(u ** 2 + w ** 2)
    Y_p = w + (MAV.rho * Va * MAV.S_wing * MAV.b / (4 * m)) * MAV.C_Y_p
    Y_r = -u + (MAV.rho * Va * MAV.S_wing * MAV.b / (4 * m)) * MAV.C_Y_r
    Y_delta_a = (MAV.rho * (Va ** 2) * MAV.S_wing / (2 * m)) * MAV.C_Y_delta_a
    Y_delta_r = (MAV.rho * (Va ** 2) * MAV.S_wing / (2 * m)) * MAV.C_Y_delta_r
    L_v = MAV.rho * MAV.S_wing * (MAV.b ** 2) * v / (4 * Va) * (MAV.C_p_p * p + MAV.C_p_r * r) \
          + MAV.rho * MAV.S_wing * MAV.b * v * (
                  MAV.C_p_0 + MAV.C_p_beta * beta + MAV.C_p_delta_a * delta_a + MAV.C_p_delta_r * delta_r) \
          + (MAV.rho * MAV.S_wing * MAV.b * MAV.C_p_beta / 2) * np.sqrt(u ** 2 + w ** 2)
    L_p = MAV.gamma1 * q + (MAV.rho * Va * MAV.S_wing * (MAV.b ** 2) / 4) * MAV.C_p_p
    L_r = -MAV.gamma2 * q + (MAV.rho * Va * MAV.S_wing * (MAV.b ** 2) / 4) * MAV.C_p_r
    L_delta_a = (MAV.rho * (Va ** 2) * MAV.S_wing * MAV.b / 2) * MAV.C_p_delta_a
    L_delta_r = (MAV.rho * (Va ** 2) * MAV.S_wing * MAV.b / 2) * MAV.C_p_delta_r
    N_v = MAV.rho * MAV.S_wing * (MAV.b ** 2) * v / (4 * Va) * (MAV.C_r_p * p + MAV.C_r_r * r) \
          + MAV.rho * MAV.S_wing * MAV.b * v * (
                  MAV.C_r_0 + MAV.C_r_beta * beta + MAV.C_r_delta_a * delta_a + MAV.C_r_delta_r * delta_r) \
          + (MAV.rho * MAV.S_wing * MAV.b * MAV.C_r_beta / 2) * np.sqrt(u ** 2 + w ** 2)
    N_p = MAV.gamma7 * q + (MAV.rho * Va * MAV.S_wing * (MAV.b ** 2) / 4) * MAV.C_r_p
    N_r = -MAV.gamma1 * q + (MAV.rho * Va * MAV.S_wing * (MAV.b ** 2) / 4) * MAV.C_r_r
    N_delta_a = (MAV.rho * (Va ** 2) * MAV.S_wing * MAV.b / 2) * MAV.C_r_delta_a
    N_delta_r = (MAV.rho * (Va ** 2) * MAV.S_wing * MAV.b / 2) * MAV.C_r_delta_r

    # Lateral State_space Model
    A_lat = np.array([[Y_v, Y_p, Y_r, MAV.gravity * np.cos(theta) * np.cos(phi), 0],
                      [L_v, L_p, L_r, 0, 0],
                      [N_v, N_p, N_r, 0, 0],
                      [0, 1, np.cos(phi) * np.tan(theta),
                       q * np.cos(phi) * np.tan(theta) - r * np.sin(phi) * np.tan(theta), 0],
                      [0, 0, np.cos(phi) / np.cos(theta),
                       p * np.cos(phi) / np.cos(theta) - r * np.sin(phi) / np.cos(theta), 0]])
    B_lat = np.array([[Y_delta_a, Y_delta_r],
                      [L_delta_a, L_delta_r],
                      [N_delta_a, N_delta_r],
                      [0, 0],
                      [0, 0]])

    # Longitudinal Model Coefficients
    # Cannot find C_X_0, C_X_alpha .etc in the aerosonde_parameters.py and the book.
    # X_u = (u * MAV.rho * MAV.S_wing / MAV.mass) * ()

    # Longitudinal State_space Model
    A_lon = None
    B_lon = None

    x_quat = quaternion_state(x_euler)
    return A_lon, B_lon, A_lat, B_lat


def euler_state(x_quat):
    # convert state x with attitude represented by quaternion
    # to x_euler with attitude represented by Euler angles
    # x_quat = [pn, pe, pd, u, v, w, e0, e1, e2, e3, p, q, r].T

    phi, theta, psi = Quaternion2Euler(x_quat[6:10])
    x_euler = np.zeros([12, 1])
    x_euler[0:6] = x_quat[0:6]
    x_euler[6:9] = np.array([[phi, theta, psi]]).T
    x_euler[9:] = x_quat[10:]

    return x_euler


def quaternion_state(x_euler):
    # convert state x_euler with attitude represented by Euler angles
    # to x_quat with attitude represented by quaternions
    # x_euler = [pn, pe, pd, u, v, w, phi, theta, psi, p, q, r].T

    e = Euler2Quaternion(x_euler[6][0], x_euler[7][0], x_euler[8][0])
    x_quat = np.zeros([13, 1])
    x_quat[0:6] = x_euler[0:6]
    x_quat[6:10] = np.array([e]).T
    x_quat[10:] = x_euler[9:]

    return x_quat


def f_euler(mav, x_euler, input):
    # return 12x1 dynamics (as if state were Euler state)
    # compute f at euler_state
    return f_euler_


def df_dx(mav, x_euler, input):
    # take partial of f_euler with respect to x_euler
    return A


def df_du(mav, x_euler, delta):
    # take partial of f_euler with respect to delta
    return B


def dT_dVa(mav, Va, delta_t):
    # returns the derivative of motor thrust with respect to Va
    return dThrust


def dT_ddelta_t(mav, Va, delta_t):
    # returns the derivative of motor thrust with respect to delta_t
    return dThrust
