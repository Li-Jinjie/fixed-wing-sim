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
from tools.tools import euler_2_quaternion, quaternion_2_euler, euler_2_rotation
from tools.transfer_function import TransferFunction
import parameters.aerosonde_parameters as MAV
from parameters.simulation_parameters import ts_simulation as Ts


def compute_tf_model(mav, trim_state, trim_input):
    # trim values
    Va = mav._Va
    u, v, w = mav._state[3:6]
    Vg = np.sqrt(u ** 2 + v ** 2 + w ** 2)
    # Follow the equations in page 68.

    # Lateral Transfer Functions:
    # 1. Roll Angle
    rV2Sb_2 = (1 / 2) * MAV.rho * (Va ** 2) * MAV.S_wing * MAV.b
    a_phi_1 = -rV2Sb_2 * MAV.C_p_p * MAV.b / (2 * Va)
    a_phi_2 = rV2Sb_2 * MAV.C_p_delta_a
    T_phi_delta_a = TransferFunction(np.array([[a_phi_2]]), np.array([[1, a_phi_1, 0]]), Ts)

    # 2. Course and Heading
    T_chi_phi = TransferFunction(np.array([[MAV.gravity / Vg]]), np.array([[1, 0]]), Ts)

    # 3. Sideslip
    rVS_2m = (MAV.rho * Va * MAV.S_wing / (2 * MAV.mass))
    a_beta_1 = - rVS_2m * MAV.C_Y_beta
    a_beta_2 = rVS_2m * MAV.C_Y_delta_r
    T_beta_delta_r = TransferFunction(np.array([[a_beta_2]]), np.array([[1, a_beta_1]]), Ts)

    # Longitudinal Transfer Functions:
    # 1. Pitch Angle
    rV2cS_2Jy = (MAV.rho * (Va ** 2) * MAV.c * MAV.S_wing / (2 * MAV.Jy))
    a_theta_1 = - rV2cS_2Jy * MAV.C_m_q * (MAV.c / (2 * Va))
    a_theta_2 = - rV2cS_2Jy * MAV.C_m_alpha
    a_theta_3 = rV2cS_2Jy * MAV.C_m_delta_e
    T_theta_delta_e = TransferFunction(np.array([[a_theta_3]]), np.array([[1, a_theta_1, a_theta_2]]), Ts)

    # 2. Altitude
    T_h_theta = TransferFunction(np.array([[Va]]), np.array([[1, 0]]), Ts)
    _, theta, _ = quaternion_2_euler(mav._state[6:10])
    T_h_Va = TransferFunction(np.array([[theta]]), np.array([[1, 0]]), Ts)

    # 3. Airspeed
    # Pay Attention: Use trim condition!!!
    u_trim, v_trim, w_trim = trim_state.item(3), trim_state.item(4), trim_state.item(5)
    Va_trim = u_trim  # u in body frame, no wind
    delta_e_trim = trim_input[1].item()
    delta_t_trim = trim_input[3].item()
    phi_trim, theta_trim, psi_trim = quaternion_2_euler(trim_state[6:10])
    alpha_trim = np.arctan2(w_trim, u_trim)
    R_trim = euler_2_rotation(phi_trim, theta_trim, psi_trim)  # R: body to inertial
    [[_], [v_i_trim], [w_i_trim]] = R_trim @ trim_state[3:6]  # in the inertial frame
    chi_trim = np.arctan2(v_i_trim, w_i_trim)  # -pi to pi

    a_V_1 = (MAV.rho * Va_trim * MAV.S_wing / MAV.mass) * \
            (MAV.C_D_0 + MAV.C_D_alpha * alpha_trim + MAV.C_D_delta_e * delta_e_trim) + \
            (MAV.rho * MAV.S_prop / MAV.mass) * MAV.C_prop * Va_trim
    a_V_2 = (MAV.rho * MAV.S_prop / MAV.mass) * MAV.C_prop * (MAV.k_motor ** 2) * delta_t_trim
    a_V_3 = MAV.gravity * np.cos(theta_trim - chi_trim)
    T_Va_delta_t = TransferFunction(np.array([[a_V_2]]), np.array([[1, a_V_1]]), Ts)
    T_Va_theta = TransferFunction(np.array([[-a_V_3]]), np.array([[1, a_V_1]]), Ts)

    return T_phi_delta_a, T_chi_phi, T_beta_delta_r, T_theta_delta_e, T_h_theta, T_h_Va, T_Va_delta_t, T_Va_theta


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

    # [e, a, r, t]
    delta_e = trim_input.item(0)
    delta_a = trim_input.item(1)
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

    # TODO: Longitudinal Model Coefficients. The linearize process should be calculated using  control.iosys.linearize()
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

    phi, theta, psi = quaternion_2_euler(x_quat[6:10])
    x_euler = np.zeros([12, 1])
    x_euler[0:6] = x_quat[0:6]
    x_euler[6:9] = np.array([[phi, theta, psi]]).T
    x_euler[9:] = x_quat[10:]

    return x_euler


def quaternion_state(x_euler):
    # convert state x_euler with attitude represented by Euler angles
    # to x_quat with attitude represented by quaternions
    # x_euler = [pn, pe, pd, u, v, w, phi, theta, psi, p, q, r].T

    e = euler_2_quaternion(x_euler[6][0], x_euler[7][0], x_euler[8][0])
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


# Read the slides of chapter 5 on: https://uavbook.byu.edu/lib/exe/fetch.php?media=lecture:chap5.pdf.
# No need in calculating the transfer functions here.
def dT_dVa(mav, Va_trim):
    # returns the derivative of motor thrust with respect to Va
    dThrust = - MAV.rho * MAV.S_prop * MAV.C_prop * Va_trim
    return dThrust


def dT_ddelta_t(mav, delta_t_trim):
    # returns the derivative of motor thrust with respect to delta_t
    dThrust = MAV.rho * MAV.S_prop * MAV.C_prop * (MAV.k_motor ** 2) * delta_t_trim
    return dThrust
