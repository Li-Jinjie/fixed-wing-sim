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
from tools.rotations import euler_2_quaternion, quaternion_2_euler
import parameters.aerosonde_parameters as MAV
from parameters.simulation_parameters import ts_simulation as Ts
from message_types.msg_delta import MsgDelta


def compute_model(mav, trim_state, trim_input):
    A_lon, B_lon, A_lat, B_lat = compute_ss_model(mav, trim_state, trim_input)
    Va_trim, alpha_trim, theta_trim, a_phi1, a_phi2, a_theta1, a_theta2, a_theta3, \
    a_V1, a_V2, a_V3 = compute_tf_model(mav, trim_state, trim_input)

    # write transfer function gains to file
    file = open('model_coef.py', 'w')
    file.write('import numpy as np\n')
    file.write('x_trim = np.array([[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]]).T\n' %
               (trim_state.item(0), trim_state.item(1), trim_state.item(2), trim_state.item(3),
                trim_state.item(4), trim_state.item(5), trim_state.item(6), trim_state.item(7),
                trim_state.item(8), trim_state.item(9), trim_state.item(10), trim_state.item(11),
                trim_state.item(12)))
    file.write('u_trim = np.array([[%f, %f, %f, %f]]).T\n' %
               (trim_input.elevator, trim_input.aileron, trim_input.rudder, trim_input.throttle))
    file.write('Va_trim = %f\n' % Va_trim)
    file.write('alpha_trim = %f\n' % alpha_trim)
    file.write('theta_trim = %f\n' % theta_trim)
    file.write('a_phi1 = %f\n' % a_phi1)
    file.write('a_phi2 = %f\n' % a_phi2)
    file.write('a_theta1 = %f\n' % a_theta1)
    file.write('a_theta2 = %f\n' % a_theta2)
    file.write('a_theta3 = %f\n' % a_theta3)
    file.write('a_V1 = %f\n' % a_V1)
    file.write('a_V2 = %f\n' % a_V2)
    file.write('a_V3 = %f\n' % a_V3)
    file.write('A_lon = np.array([\n    [%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f]])\n' %
               (A_lon[0][0], A_lon[0][1], A_lon[0][2], A_lon[0][3], A_lon[0][4],
                A_lon[1][0], A_lon[1][1], A_lon[1][2], A_lon[1][3], A_lon[1][4],
                A_lon[2][0], A_lon[2][1], A_lon[2][2], A_lon[2][3], A_lon[2][4],
                A_lon[3][0], A_lon[3][1], A_lon[3][2], A_lon[3][3], A_lon[3][4],
                A_lon[4][0], A_lon[4][1], A_lon[4][2], A_lon[4][3], A_lon[4][4]))
    file.write('B_lon = np.array([\n    [%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f]])\n' %
               (B_lon[0][0], B_lon[0][1],
                B_lon[1][0], B_lon[1][1],
                B_lon[2][0], B_lon[2][1],
                B_lon[3][0], B_lon[3][1],
                B_lon[4][0], B_lon[4][1],))
    file.write('A_lat = np.array([\n    [%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f],\n    '
               '[%f, %f, %f, %f, %f]])\n' %
               (A_lat[0][0], A_lat[0][1], A_lat[0][2], A_lat[0][3], A_lat[0][4],
                A_lat[1][0], A_lat[1][1], A_lat[1][2], A_lat[1][3], A_lat[1][4],
                A_lat[2][0], A_lat[2][1], A_lat[2][2], A_lat[2][3], A_lat[2][4],
                A_lat[3][0], A_lat[3][1], A_lat[3][2], A_lat[3][3], A_lat[3][4],
                A_lat[4][0], A_lat[4][1], A_lat[4][2], A_lat[4][3], A_lat[4][4]))
    file.write('B_lat = np.array([\n    [%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f],\n    '
               '[%f, %f]])\n' %
               (B_lat[0][0], B_lat[0][1],
                B_lat[1][0], B_lat[1][1],
                B_lat[2][0], B_lat[2][1],
                B_lat[3][0], B_lat[3][1],
                B_lat[4][0], B_lat[4][1],))
    file.write('Ts = %f\n' % Ts)
    file.close()


def compute_tf_model(mav, trim_state, trim_input):
    # please refer to page 26 of the addendum by McLain
    # trim values
    mav._state = trim_state
    mav._update_velocity_data()
    Va_trim = mav._Va
    alpha_trim = mav._alpha
    phi, theta_trim, psi = quaternion_2_euler(trim_state[6:10])

    # define transfer function constants
    rV2Sb_2 = (1. / 2.) * MAV.rho * (Va_trim ** 2) * MAV.S_wing * MAV.b
    a_phi1 = -rV2Sb_2 * MAV.C_p_p * MAV.b / (2. * Va_trim)
    a_phi2 = rV2Sb_2 * MAV.C_p_delta_a

    rV2cS_2Jy = (MAV.rho * (Va_trim ** 2) * MAV.c * MAV.S_wing / (2. * MAV.Jy))
    a_theta1 = - rV2cS_2Jy * MAV.C_m_q * (MAV.c / (2. * Va_trim))
    a_theta2 = - rV2cS_2Jy * MAV.C_m_alpha
    a_theta3 = rV2cS_2Jy * MAV.C_m_delta_e

    # Compute transfer function coefficients using new propulsion model
    a_V1 = (MAV.rho * Va_trim * MAV.S_wing / MAV.mass) * \
           (MAV.C_D_0 + MAV.C_D_alpha * alpha_trim + MAV.C_D_delta_e * trim_input.elevator) - \
           (1. / MAV.mass) * dT_dVa(mav, Va_trim, trim_input.throttle)
    a_V2 = (1. / MAV.mass) * dT_ddelta_t(mav, Va_trim, trim_input.throttle)
    a_V3 = MAV.gravity * np.cos(theta_trim - alpha_trim)

    return Va_trim, alpha_trim, theta_trim, a_phi1, a_phi2, a_theta1, a_theta2, a_theta3, a_V1, a_V2, a_V3


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

    delta_e = trim_input.elevator
    delta_a = trim_input.aileron
    delta_r = trim_input.rudder
    delta_t = trim_input.throttle

    Va = mav._Va
    alpha = mav._alpha
    beta = mav._beta

    # m: External moment applied to the airframe about the body frame y-axis.
    rVS_2 = MAV.rho * (Va ** 2) * MAV.S_wing / 2.
    m = rVS_2 * MAV.c * (MAV.C_m_0 + MAV.C_m_alpha * alpha +
                         MAV.C_m_q * MAV.c / (2. * Va) * q + MAV.C_m_delta_e * delta_e)

    # Lateral Model Coefficients
    Y_v = MAV.rho * MAV.S_wing * MAV.b * v / (4. * m * Va) * (MAV.C_Y_p * p + MAV.C_Y_r * r) \
          + (MAV.rho * MAV.S_wing * v / MAV.mass) * (
                  MAV.C_Y_0 + MAV.C_Y_beta * beta + MAV.C_Y_delta_a * delta_a + MAV.C_Y_delta_r * delta_r) \
          + (MAV.rho * MAV.S_wing * MAV.C_Y_beta / (2. * MAV.mass)) * np.sqrt(u ** 2 + w ** 2)
    Y_p = w + (MAV.rho * Va * MAV.S_wing * MAV.b / (4. * m)) * MAV.C_Y_p
    Y_r = -u + (MAV.rho * Va * MAV.S_wing * MAV.b / (4. * m)) * MAV.C_Y_r
    Y_delta_a = (MAV.rho * (Va ** 2) * MAV.S_wing / (2. * MAV.mass)) * MAV.C_Y_delta_a
    Y_delta_r = (MAV.rho * (Va ** 2) * MAV.S_wing / (2. * MAV.mass)) * MAV.C_Y_delta_r
    L_v = MAV.rho * MAV.S_wing * (MAV.b ** 2) * v / (4. * Va) * (MAV.C_p_p * p + MAV.C_p_r * r) \
          + MAV.rho * MAV.S_wing * MAV.b * v * (
                  MAV.C_p_0 + MAV.C_p_beta * beta + MAV.C_p_delta_a * delta_a + MAV.C_p_delta_r * delta_r) \
          + (MAV.rho * MAV.S_wing * MAV.b * MAV.C_p_beta / 2.) * np.sqrt(u ** 2 + w ** 2)
    L_p = MAV.gamma1 * q + (MAV.rho * Va * MAV.S_wing * (MAV.b ** 2) / 4.) * MAV.C_p_p
    L_r = -MAV.gamma2 * q + (MAV.rho * Va * MAV.S_wing * (MAV.b ** 2) / 4.) * MAV.C_p_r
    L_delta_a = (MAV.rho * (Va ** 2) * MAV.S_wing * MAV.b / 2.) * MAV.C_p_delta_a
    L_delta_r = (MAV.rho * (Va ** 2) * MAV.S_wing * MAV.b / 2.) * MAV.C_p_delta_r
    N_v = MAV.rho * MAV.S_wing * (MAV.b ** 2) * v / (4. * Va) * (MAV.C_r_p * p + MAV.C_r_r * r) \
          + MAV.rho * MAV.S_wing * MAV.b * v * (
                  MAV.C_r_0 + MAV.C_r_beta * beta + MAV.C_r_delta_a * delta_a + MAV.C_r_delta_r * delta_r) \
          + (MAV.rho * MAV.S_wing * MAV.b * MAV.C_r_beta / 2.) * np.sqrt(u ** 2 + w ** 2)
    N_p = MAV.gamma7 * q + (MAV.rho * Va * MAV.S_wing * (MAV.b ** 2) / 4.) * MAV.C_r_p
    N_r = -MAV.gamma1 * q + (MAV.rho * Va * MAV.S_wing * (MAV.b ** 2) / 4.) * MAV.C_r_r
    N_delta_a = (MAV.rho * (Va ** 2) * MAV.S_wing * MAV.b / 2.) * MAV.C_r_delta_a
    N_delta_r = (MAV.rho * (Va ** 2) * MAV.S_wing * MAV.b / 2.) * MAV.C_r_delta_r

    A_lat_14 = MAV.gravity * np.cos(theta) * np.cos(phi)
    A_lat_43 = np.cos(phi) * np.tan(theta)
    A_lat_44 = q * np.cos(phi) * np.tan(theta) - r * np.sin(phi) * np.tan(theta)
    A_lat_53 = np.cos(phi) / np.cos(theta)
    A_lat_54 = p * np.cos(phi) / np.cos(theta) - r * np.sin(phi) / np.cos(theta)

    # Lateral State_space Model
    # extract lateral states (v, p, r, phi, psi)
    A_lat = np.array([[Y_v, Y_p, Y_r, A_lat_14, 0.],
                      [L_v, L_p, L_r, 0., 0.],
                      [N_v, N_p, N_r, 0., 0.],
                      [0., 1., A_lat_43, A_lat_44, 0.],
                      [0., 0., A_lat_53, A_lat_54, 0.]])
    B_lat = np.array([[Y_delta_a, Y_delta_r],
                      [L_delta_a, L_delta_r],
                      [N_delta_a, N_delta_r],
                      [0., 0.],
                      [0., 0.]])

    # Longitudinal Model Coefficients
    # refer to page 49 on uavbook
    [[C_X_0, C_X_q, C_X_alpha, C_X_delta_e],
     [C_Z_0, C_Z_q, C_Z_alpha, C_Z_delta_e]] = np.array([[np.cos(alpha), -np.sin(alpha)],
                                                         [np.sin(alpha), np.cos(alpha)]]) @ \
                                               np.array([[-MAV.C_D_0, -MAV.C_D_q, -MAV.C_D_alpha, -MAV.C_D_delta_e],
                                                         [-MAV.C_L_0, -MAV.C_L_q, -MAV.C_L_alpha, -MAV.C_L_delta_e]])

    dT_du = dT_dVa(mav, Va, delta_t) * u / Va  # I calculate this. --LJJ
    dT_dw = dT_dVa(mav, Va, delta_t) * w / Va
    X_u = (u * MAV.rho * MAV.S_wing / MAV.mass) * (C_X_0 + C_X_alpha * alpha + C_X_delta_e * delta_e) \
          - (MAV.rho * MAV.S_wing * C_X_alpha * w / (2. * MAV.mass)) \
          + (MAV.rho * MAV.S_wing * MAV.c * C_X_q * u * q / (4. * MAV.mass * Va)) + dT_du

    X_w = -q + w * MAV.rho * MAV.S_wing / MAV.mass * (C_X_0 + C_X_alpha * alpha + C_X_delta_e * delta_e) \
          + (MAV.rho * MAV.S_wing * C_X_alpha * u / (2. * MAV.mass)) \
          + (MAV.rho * MAV.S_wing * MAV.c * C_X_q * w * q / (4. * MAV.mass * Va)) + dT_dw
    X_q = -w + (MAV.rho * Va * MAV.S_wing * C_X_q * MAV.c / (4. * MAV.mass))
    X_delta_e = (MAV.rho * (Va ** 2) * MAV.S_wing * C_X_delta_e / (2. * MAV.mass))
    X_delta_t = dT_ddelta_t(mav, Va, delta_t)
    Z_u = q + (u * MAV.rho * MAV.S_wing / MAV.mass) * (C_Z_0 + C_Z_alpha * alpha + C_Z_delta_e * delta_e) \
          - (MAV.rho * MAV.S_wing * C_Z_alpha * w / (2. * MAV.mass)) \
          + (u * MAV.rho * MAV.S_wing * C_Z_q * MAV.c * q / (4. * MAV.mass * Va))
    Z_w = (w * MAV.rho * MAV.S_wing / MAV.mass) * (C_Z_0 + C_Z_alpha * alpha + C_Z_delta_e * delta_e) \
          + (MAV.rho * MAV.S_wing * C_Z_alpha * u / (2. * MAV.mass)) \
          + (w * MAV.rho * MAV.S_wing * C_Z_q * MAV.c * q / (4. * MAV.mass * Va))
    Z_q = u + (MAV.rho * Va * MAV.S_wing * C_Z_q * MAV.c / (4. * MAV.mass))
    Z_delta_e = (MAV.rho * (Va ** 2) * MAV.S_wing * C_Z_delta_e / (2. * MAV.mass))
    M_u = u * MAV.rho * MAV.S_wing * MAV.c / MAV.Jy * (MAV.C_m_0 + MAV.C_m_alpha * alpha + MAV.C_m_delta_e * delta_e) \
          - (MAV.rho * MAV.S_wing * MAV.c * MAV.C_m_alpha * w / (2. * MAV.Jy)) \
          + (MAV.rho * MAV.S_wing * (MAV.c ** 2.) * MAV.C_m_q * q * u / (4. * MAV.Jy * Va))
    M_w = (w * MAV.rho * MAV.S_wing * MAV.c / MAV.Jy) * (MAV.C_m_0 + MAV.C_m_alpha * alpha + MAV.C_m_delta_e * delta_e) \
          + (MAV.rho * MAV.S_wing * MAV.c * MAV.C_m_alpha * u / (2. * MAV.Jy)) \
          + (MAV.rho * MAV.S_wing * (MAV.c ** 2.) * MAV.C_m_q * q * w / (4. * MAV.Jy * Va))
    M_q = (MAV.rho * Va * MAV.S_wing * (MAV.c ** 2.) * MAV.C_m_q / (4. * MAV.Jy))
    M_delta_e = (MAV.rho * (Va ** 2.) * MAV.S_wing * MAV.c * MAV.C_m_delta_e / (2. * MAV.Jy))

    A_lon_12 = X_w * Va * np.cos(alpha)
    A_lon_14 = -MAV.gravity * np.cos(theta)
    A_lon_21 = Z_u / (Va * np.cos(alpha))
    A_lon_23 = Z_q / (Va * np.cos(alpha))
    A_lon_24 = -MAV.gravity * np.sin(theta) / (Va * np.cos(alpha))
    A_lon_32 = M_w * Va * np.cos(alpha)
    A_lon_52 = -Va * np.cos(theta) * np.cos(alpha)
    A_lon_54 = u * np.cos(theta) + w * np.sin(theta)

    B_lon_21 = Z_delta_e / (Va * np.cos(alpha))

    # extract longitudinal states (u, w, q, theta, pd) and change pd to h
    A_lon = np.array([[X_u, A_lon_12, X_q, A_lon_14, 0.],
                      [A_lon_21, Z_w, A_lon_23, A_lon_24, 0.],
                      [M_u, A_lon_32, M_q, 0., 0.],
                      [0., 0., 1., 0., 0.],
                      [np.sin(theta), A_lon_52, 0., A_lon_54, 0.]])
    B_lon = np.array([[X_delta_e, X_delta_t],
                      [B_lon_21, 0.],
                      [M_delta_e, 0],
                      [0., 0.],
                      [0., 0.]])

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


# def f_euler(mav, x_euler, delta):
#     # return 12x1 dynamics (as if state were Euler state)
#     # compute f at euler_state
#
#     f_euler_ =
#     return f_euler_
#
#
# def df_dx(mav, x_euler, delta):
#     # take partial of f_euler with respect to x_euler
#     A =
#     return A
#
#
# def df_du(mav, x_euler, delta):
#     # take partial of f_euler with respect to input
#     B =
#     return B


def dT_dVa(mav, Va, delta_t):
    # returns the derivative of motor thrust with respect to Va
    eps = 1e-5
    T_eps, Q_eps = mav._motor_thrust_torque(Va + eps, delta_t)
    T, Q = mav._motor_thrust_torque(Va, delta_t)
    return (T_eps - T) / eps


def dT_ddelta_t(mav, Va, delta_t):
    # returns the derivative of motor thrust with respect to delta_t
    eps = 1e-5
    T_eps, Q_eps = mav._motor_thrust_torque(Va, delta_t + eps)
    T, Q = mav._motor_thrust_torque(Va, delta_t)
    return (T_eps - T) / eps
