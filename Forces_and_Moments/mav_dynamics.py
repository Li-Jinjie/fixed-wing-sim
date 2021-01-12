"""
mav_dynamics
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state
    
"""
import sys

sys.path.append('..')
import numpy as np

# load message types
from message_types.msg_state import msg_state

import parameters.aerosonde_parameters as MAV
from tools.tools import *


class mav_dynamics:
    def __init__(self, Ts):
        self._ts_simulation = Ts
        # set initial states based on parameter file
        # _state is the 13x1 internal state of the aircraft that is being propagated:
        # _state = [pn, pe, pd, u, v, w, e0, e1, e2, e3, p, q, r]
        # We will also need a variety of other elements that are functions of the _state and the wind.
        # self.true_state is a 19x1 vector that is estimated and used by the autopilot to control the aircraft:
        # true_state = [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        self._state = np.array([[MAV.pn0],  # (0)
                                [MAV.pe0],  # (1)
                                [MAV.pd0],  # (2)
                                [MAV.u0],  # (3)
                                [MAV.v0],  # (4)
                                [MAV.w0],  # (5)
                                [MAV.e0],  # (6)
                                [MAV.e1],  # (7)
                                [MAV.e2],  # (8)
                                [MAV.e3],  # (9)
                                [MAV.p0],  # (10)
                                [MAV.q0],  # (11)
                                [MAV.r0]])  # (12)
        # store wind data for fast recall since it is used at various points in simulation
        self._wind = np.array([[0.], [0.], [0.]])  # wind in NED frame in meters/sec
        self._update_velocity_data()
        # store forces to avoid recalculation in the sensors function
        self._forces = np.array([[0.], [0.], [0.]])
        self._Va = MAV.u0
        self._alpha = 0
        self._beta = 0
        # initialize true_state message
        self.msg_true_state = msg_state()

    ###################################
    # public functions
    def update_state(self, delta, wind):
        '''
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        '''
        # get forces and moments acting on rigid bod
        forces_moments = self._forces_moments(delta)

        # Integrate ODE using Runge-Kutta RK4 algorithm
        # TODO: The formulas of k2 and k3 are different from the definitions in wikipedia.
        time_step = self._ts_simulation
        k1 = self._derivatives(self._state, forces_moments)
        k2 = self._derivatives(self._state + time_step / 2. * k1, forces_moments)
        k3 = self._derivatives(self._state + time_step / 2. * k2, forces_moments)
        k4 = self._derivatives(self._state + time_step * k3, forces_moments)
        self._state += time_step / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

        # normalize the quaternion
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        normE = np.sqrt(e0 ** 2 + e1 ** 2 + e2 ** 2 + e3 ** 2)
        self._state[6][0] = self._state.item(6) / normE
        self._state[7][0] = self._state.item(7) / normE
        self._state[8][0] = self._state.item(8) / normE
        self._state[9][0] = self._state.item(9) / normE

        # update the airspeed, angle of attack, and side slip angles using new state
        self._update_velocity_data(wind)

        # update the message class for the true state
        self._update_msg_true_state()

    ###################################
    # private functions
    def _derivatives(self, state, forces_moments):
        """
        for the dynamics xdot = f(x, u), returns f(x, u)
        """
        # extract the states
        pn = state.item(0)
        pe = state.item(1)
        pd = state.item(2)
        u = state.item(3)
        v = state.item(4)
        w = state.item(5)
        e0 = state.item(6)
        e1 = state.item(7)
        e2 = state.item(8)
        e3 = state.item(9)
        p = state.item(10)
        q = state.item(11)
        r = state.item(12)
        #   extract forces/moments
        fx = forces_moments.item(0)
        fy = forces_moments.item(1)
        fz = forces_moments.item(2)
        l = forces_moments.item(3)
        m = forces_moments.item(4)
        n = forces_moments.item(5)

        # From page 256, Appendix B.2
        # B.1 - B.4
        # position kinematics
        pn_dot = (e1 ** 2 + e0 ** 2 - e2 ** 2 - e3 ** 2) * u + 2 * (e1 * e2 - e3 * e0) * v + 2 * (e1 * e3 + e2 * e0) * w
        pe_dot = 2 * (e1 * e2 + e3 * e0) * u + (e2 ** 2 + e0 ** 2 - e1 ** 2 - e3 ** 2) * v + 2 * (e2 * e3 - e1 * e0) * w
        pd_dot = 2 * (e1 * e3 - e2 * e0) * u + 2 * (e2 * e3 + e1 * e0) * v + (e3 ** 2 + e0 ** 2 - e1 ** 2 - e2 ** 2) * w

        # position dynamics
        u_dot = r * v - q * w + fx / MAV.mass
        v_dot = p * w - r * u + fy / MAV.mass
        w_dot = q * u - p * v + fz / MAV.mass

        # rotational kinematics
        e0_dot = (0 - p * e1 - q * e2 - r * e3) / 2
        e1_dot = (p * e0 + 0 + r * e2 - q * e3) / 2
        e2_dot = (q * e0 - r * e1 + 0 + p * e3) / 2
        e3_dot = (r * e0 + q * e1 - p * e2 + 0) / 2

        # rotatonal dynamics
        p_dot = MAV.gamma1 * p * q - MAV.gamma2 * q * r + MAV.gamma3 * l + MAV.gamma4 * n
        q_dot = MAV.gamma5 * p * r - MAV.gamma6 * (p ** 2 - r ** 2) + m / MAV.Jy
        r_dot = MAV.gamma7 * p * q - MAV.gamma1 * q * r + MAV.gamma4 * l + MAV.gamma8 * n

        # collect the derivative of the states
        x_dot = np.array([[pn_dot, pe_dot, pd_dot, u_dot, v_dot, w_dot,
                           e0_dot, e1_dot, e2_dot, e3_dot, p_dot, q_dot, r_dot]]).T
        return x_dot

    def _update_velocity_data(self, wind=np.zeros((6, 1))):
        # the formula in page 57
        if wind[0] != 0:  # if activate the wind simulation
            # TODO: coordinate transformation using quaternion
            phi, theta, psi = Quaternion2Euler(self._state[6:10])
            R = Euler2Rotation(phi, theta, psi)  # R: body to inertial  R.T: inertial to body
            V_w = R.T @ wind[0:3] + wind[3:]  # in the body frame
        else:
            V_w = np.zeros((3, 1))

        u_r = self._state[3][0] - V_w[0][0]
        v_r = self._state[4][0] - V_w[1][0]
        w_r = self._state[5][0] - V_w[2][0]

        # compute airspeed
        self._Va = np.sqrt(u_r ** 2 + v_r ** 2 + w_r ** 2)
        # compute angle of attack
        self._alpha = np.arctan2(w_r, u_r)
        # compute sideslip angle
        self._beta = np.arcsin(v_r / self._Va)

    def _forces_moments(self, delta):
        """
        return the forces on the UAV based on the state, wind, and control surfaces
        :param delta: np.matrix(delta_a, delta_e, delta_r, delta_t)
        :return: Forces_and_Moments on the UAV np.matrix(Fx, Fy, Fz, Ml, Mn, Mm)
        """
        delta_a = delta[0]
        delta_e = delta[1]
        delta_r = delta[2]
        delta_t = delta[3]

        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)

        p = self._state.item(10)
        q = self._state.item(11)
        r = self._state.item(12)

        # ---- calculate C_L, in page 63 ----
        term_1 = np.exp(-MAV.M * (self._alpha - MAV.alpha0))
        term_2 = np.exp(MAV.M * (self._alpha + MAV.alpha0))
        sigma = (1 + term_1 + term_2) / ((1 + term_1) * (1 + term_2))
        C_L = (1 - sigma) * (MAV.C_L_0 + MAV.C_L_alpha * self._alpha) + sigma * (
                2 * np.sign(self._alpha) * (np.sin(self._alpha) ** 2) * np.cos(self._alpha))

        # C_L = MAV.C_L_0 + MAV.C_L_alpha * self._alpha   # linear model

        # ---- calculate C_D, in page 63 ----
        C_D = MAV.C_D_p + (MAV.C_L_0 + MAV.C_L_alpha * self._alpha) ** 2 / (np.pi * MAV.e * MAV.AR)
        # The e here is the Oswald efficiency factor, NOT the euler number.

        # C_D = MAV.C_D_0 + MAV.C_D_alpha * self._alpha   # linear model

        # ---- calculate C_X(alpha).etc, in page 62 ----
        c_alpha = np.cos(self._alpha)
        s_alpha = np.sin(self._alpha)
        C_X = - C_D * c_alpha + C_L * s_alpha
        C_X_q = - MAV.C_D_q * c_alpha + MAV.C_L_q * s_alpha
        C_X_delta_e = -MAV.C_D_delta_e * c_alpha + MAV.C_L_delta_e * s_alpha
        C_Z = - C_D * s_alpha - C_L * c_alpha
        C_Z_q = - MAV.C_D_q * s_alpha - MAV.C_L_q * c_alpha
        C_Z_delta_e = - MAV.C_D_delta_e * s_alpha - MAV.C_L_delta_e * c_alpha

        # The first term is gravity. Please read page 257.
        # TODO: check
        rVS_2 = MAV.rho * (self._Va ** 2) * MAV.S_wing / 2

        fx = MAV.mass * MAV.gravity * 2 * (e1 * e3 - e2 * e0) \
             + rVS_2 * (C_X + C_X_q * MAV.c / (2 * self._Va) * q + C_X_delta_e * delta_e) \
             + 1 / 2 * MAV.rho * MAV.S_prop * MAV.C_prop * ((MAV.k_motor * delta_t) ** 2 - self._Va ** 2)

        fy = MAV.mass * MAV.gravity * 2 * (e2 * e3 + e1 * e0) \
             + rVS_2 * (MAV.C_Y_0 + MAV.C_Y_beta * self._beta + MAV.C_Y_p * MAV.b / (2 * self._Va) * p +
                        MAV.C_Y_r * MAV.b / (2 * self._Va) * r + MAV.C_Y_delta_a * delta_a + MAV.C_Y_delta_r * delta_r)

        fz = MAV.mass * MAV.gravity * (e3 ** 2 + e0 ** 2 - e1 ** 2 - e2 ** 2) \
             + rVS_2 * (C_Z + C_Z_q * MAV.c / (2 * self._Va) * q + C_Z_delta_e * delta_e)

        self._forces[0] = fx
        self._forces[1] = fy
        self._forces[2] = fz

        # Pay attention: to distinguish between L(lift) and l(the torque along x axis), denote l as ell.
        Mx = rVS_2 * MAV.b * \
             (MAV.C_ell_0 + MAV.C_ell_beta * self._beta + MAV.C_ell_p * MAV.b / (2 * self._Va) * p +
              MAV.C_ell_r * MAV.b / (2 * self._Va) * r + MAV.C_ell_delta_a * delta_a + MAV.C_ell_delta_r * delta_r) \
             - MAV.kTp * (MAV.kOmega * delta_t) ** 2

        My = rVS_2 * MAV.c * (MAV.C_m_0 + MAV.C_m_alpha * self._alpha +
                              MAV.C_m_q * MAV.c / (2 * self._Va) * q + MAV.C_m_delta_e * delta_e)

        Mz = rVS_2 * MAV.b * \
             (MAV.C_n_0 + MAV.C_n_beta * self._beta + MAV.C_n_p * MAV.b / (2 * self._Va) * p +
              MAV.C_n_r * MAV.b / (2 * self._Va) * r + MAV.C_n_delta_a * delta_a + MAV.C_n_delta_r * delta_r)

        return np.array([[fx, fy, fz, Mx, My, Mz]]).T

    def _update_msg_true_state(self):
        # update the class structure for the true state:
        #   [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        self.msg_true_state.pn = self._state.item(0)
        self.msg_true_state.pe = self._state.item(1)
        self.msg_true_state.h = -self._state.item(2)
        self.msg_true_state.Va = self._Va
        self.msg_true_state.alpha = self._alpha
        self.msg_true_state.beta = self._beta
        self.msg_true_state.phi = phi
        self.msg_true_state.theta = theta
        self.msg_true_state.psi = psi
        # self.msg_true_state.Vg =
        # self.msg_true_state.gamma =
        # self.msg_true_state.chi =
        self.msg_true_state.p = self._state.item(10)
        self.msg_true_state.q = self._state.item(11)
        self.msg_true_state.r = self._state.item(12)
        self.msg_true_state.wn = self._wind.item(0)
        self.msg_true_state.we = self._wind.item(1)
