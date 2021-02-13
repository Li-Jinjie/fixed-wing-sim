"""
mavDynamics 
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state
    
mavsim_python
    - Beard & McLain, PUP, 2012
    - Update history:  
        2/24/2020 - RWB
"""
import sys

sys.path.append('..')
import numpy as np

# load message types
from message_types.msg_state import MsgState
from message_types.msg_sensors import MsgSensors
from message_types.msg_delta import MsgDelta

import parameters.aerosonde_parameters as MAV
import parameters.sensor_parameters as SENSOR
from tools.rotations import quaternion_2_rotation, quaternion_2_euler, euler_2_rotation


class MavDynamics:
    def __init__(self, Ts):
        self._ts_simulation = Ts
        # set initial states based on parameter file
        # _state is the 13x1 internal state of the aircraft that is being propagated:
        # _state = [pn, pe, pd, u, v, w, e0, e1, e2, e3, p, q, r]
        # We will also need a variety of other elements that are functions of the _state and the wind.
        # self.true_state is a 19x1 vector that is estimated and used by the autopilot to control the aircraft:
        # true_state = [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        self._state = np.array([[MAV.north0],  # (0)
                                [MAV.east0],  # (1)
                                [MAV.down0],  # (2)
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
        # store forces to avoid recalculation in the sensors function
        self._forces = np.array([[0.], [0.], [0.]])
        self._Va = MAV.u0
        self._alpha = 0
        self._beta = 0
        # initialize true_state message
        self.true_state = MsgState()
        # initialize the sensors message
        self._sensors = MsgSensors()
        # random walk parameters for GPS
        self._gps_eta_n = 0.
        self._gps_eta_e = 0.
        self._gps_eta_h = 0.
        # timer so that gps only updates every ts_gps seconds
        self._t_gps = 999.  # large value ensures gps updates at initial time.
        # update velocity data and forces and moments
        self._update_velocity_data()
        self._forces_moments(delta=MsgDelta())

    ###################################
    # public functions
    def update(self, delta, wind):
        '''
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        '''
        # get forces and moments acting on rigid bod
        forces_moments = self._forces_moments(delta)

        # Integrate ODE using Runge-Kutta RK4 algorithm
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
        self._update_true_state()

    def sensors(self):
        "Return value of sensors on MAV: gyros, accels, absolute_pressure, dynamic_pressure, GPS"
        # true value
        north = self._state.item(0)
        east = self._state.item(1)
        h = - self._state.item(2)  # h = -pd
        p = self._state.item(10)
        q = self._state.item(11)
        r = self._state.item(12)
        phi, theta, psi = quaternion_2_euler(self._state[6:10])
        # pdot =

        # simulate rate gyros(units are rad / sec)
        self._sensors.gyro_x = p + np.random.normal(loc=0, scale=SENSOR.gyro_sigma) + SENSOR.gyro_x_bias
        self._sensors.gyro_y = q + np.random.normal(loc=0, scale=SENSOR.gyro_sigma) + SENSOR.gyro_y_bias
        self._sensors.gyro_z = r + np.random.normal(loc=0, scale=SENSOR.gyro_sigma) + SENSOR.gyro_z_bias

        # simulate accelerometers(units of g)
        # self._forces: x,y,z
        self._sensors.accel_x = self._forces.item(0) / MAV.mass + MAV.gravity * np.sin(theta) + \
                                np.random.normal(loc=0, scale=SENSOR.accel_sigma)
        self._sensors.accel_y = self._forces.item(1) / MAV.mass - MAV.gravity * np.cos(theta) * np.sin(phi) + \
                                np.random.normal(loc=0, scale=SENSOR.accel_sigma)
        self._sensors.accel_z = self._forces.item(2) / MAV.mass - MAV.gravity * np.cos(theta) * np.cos(phi) + \
                                np.random.normal(loc=0, scale=SENSOR.accel_sigma)
        # simulate magnetometers
        # magnetic field in provo has magnetic declination of 12.5 degrees
        # and magnetic inclination of 66 degrees
        # Beijing: declination: -9.0 degrees, inclination: 57.5 degrees.
        # From https://www.magnetic-declination.com/China/Beijing/388314.html
        R_mag = euler_2_rotation(0, np.radians(57.5), np.radians(-9.0))
        # magnetic field in inertial frame: unit vector
        mag_inertial = R_mag @ np.array([[1], [0], [0]])
        R = euler_2_rotation(phi, theta, psi)  # body to inertial, need to inverse.
        # magnetic field in body frame: unit vector
        mag_body = R.T @ mag_inertial
        self._sensors.mag_x = mag_body.item(0) + SENSOR.mag_beta + np.random.normal(loc=0, scale=SENSOR.mag_sigma)
        self._sensors.mag_y = mag_body.item(1) + SENSOR.mag_beta + np.random.normal(loc=0, scale=SENSOR.mag_sigma)
        self._sensors.mag_z = mag_body.item(2) + SENSOR.mag_beta + np.random.normal(loc=0, scale=SENSOR.mag_sigma)

        # simulate pressure sensors
        self._sensors.abs_pressure = MAV.rho * MAV.gravity * (h - 0) + \
                                     np.random.normal(loc=0, scale=SENSOR.abs_pres_sigma)
        self._sensors.diff_pressure = MAV.rho * (self._Va ** 2) / 2 + \
                                      np.random.normal(loc=0, scale=SENSOR.diff_pres_sigma)
        # simulate GPS sensor
        if self._t_gps >= SENSOR.ts_gps:
            coe_1 = np.exp(-SENSOR.gps_k * SENSOR.ts_gps)
            self._gps_eta_n = coe_1 * self._gps_eta_n + np.random.normal(loc=0, scale=SENSOR.gps_n_sigma)
            self._gps_eta_e = coe_1 * self._gps_eta_e + np.random.normal(loc=0, scale=SENSOR.gps_e_sigma)
            self._gps_eta_h = coe_1 * self._gps_eta_h + np.random.normal(loc=0, scale=SENSOR.gps_h_sigma)
            self._sensors.gps_n = north + self._gps_eta_n
            self._sensors.gps_e = east + self._gps_eta_e
            self._sensors.gps_h = h + self._gps_eta_h
            # wind: N E D
            self._sensors.gps_Vg = (np.sqrt(((self._Va * np.cos(psi) + self._wind[0]) ** 2)
                                            + ((self._Va * np.sin(psi) + self._wind[1]) ** 2)) +
                                    np.random.normal(loc=0, scale=SENSOR.gps_Vg_sigma)).item()
            self._sensors.gps_course = \
                (np.arctan2(self._Va * np.sin(psi) + self._wind[1], self._Va * np.cos(psi) + self._wind[0]) +
                 np.random.normal(loc=0, scale=SENSOR.gps_course_sigma)).item()
            self._t_gps = 0.
        else:
            self._t_gps += self._ts_simulation
        return self._sensors

    def external_set_state(self, new_state):
        self._state = new_state

    ###################################
    # private functions
    def _derivatives(self, state, forces_moments):
        """
        for the dynamics xdot = f(x, u), returns f(x, u)
        """
        # extract the states
        # north = state.item(0)
        # east = state.item(1)
        # down = state.item(2)
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
        north_dot = (e1 ** 2 + e0 ** 2 - e2 ** 2 - e3 ** 2) * u + 2 * (e1 * e2 - e3 * e0) * v + \
                    2 * (e1 * e3 + e2 * e0) * w
        east_dot = 2 * (e1 * e2 + e3 * e0) * u + (e2 ** 2 + e0 ** 2 - e1 ** 2 - e3 ** 2) * v + \
                   2 * (e2 * e3 - e1 * e0) * w
        down_dot = 2 * (e1 * e3 - e2 * e0) * u + 2 * (e2 * e3 + e1 * e0) * v + \
                   (e3 ** 2 + e0 ** 2 - e1 ** 2 - e2 ** 2) * w

        # position dynamics
        u_dot = r * v - q * w + fx / MAV.mass
        v_dot = p * w - r * u + fy / MAV.mass
        w_dot = q * u - p * v + fz / MAV.mass

        # rotational kinematics
        e0_dot = (0 - p * e1 - q * e2 - r * e3) / 2
        e1_dot = (p * e0 + 0 + r * e2 - q * e3) / 2
        e2_dot = (q * e0 - r * e1 + 0 + p * e3) / 2
        e3_dot = (r * e0 + q * e1 - p * e2 + 0) / 2

        # rotational dynamics
        p_dot = MAV.gamma1 * p * q - MAV.gamma2 * q * r + MAV.gamma3 * l + MAV.gamma4 * n
        q_dot = MAV.gamma5 * p * r - MAV.gamma6 * (p ** 2 - r ** 2) + m / MAV.Jy
        r_dot = MAV.gamma7 * p * q - MAV.gamma1 * q * r + MAV.gamma4 * l + MAV.gamma8 * n

        # collect the derivative of the states
        x_dot = np.array([[north_dot, east_dot, down_dot, u_dot, v_dot, w_dot,
                           e0_dot, e1_dot, e2_dot, e3_dot, p_dot, q_dot, r_dot]]).T
        return x_dot

    def _update_velocity_data(self, wind=np.zeros((6, 1))):
        # the formula in page 57
        if wind[0][0] != 0:  # if activate the wind simulation
            steady_state = wind[0:3]
            gust = wind[3:6]
            # convert wind vector from world to body frame and add gust
            # TODO: directly coordinate transformation using quaternion
            R = quaternion_2_rotation(self._state[6:10])  # R: body to inertial  R.T: inertial to body
            V_wind_body_frame = R.T @ steady_state + gust  # in the body frame
            self._wind = R @ V_wind_body_frame  # wind in the inertial frame

        else:
            V_wind_body_frame = np.zeros((3, 1))

        # velocity vector relative to the airmass
        u_r = self._state[3][0] - V_wind_body_frame[0][0]
        v_r = self._state[4][0] - V_wind_body_frame[1][0]
        w_r = self._state[5][0] - V_wind_body_frame[2][0]
        # compute airspeed
        self._Va = np.sqrt(u_r ** 2 + v_r ** 2 + w_r ** 2)
        # compute angle of attack
        if u_r == 0:
            self._alpha = np.pi / 2  # check this line
        else:
            self._alpha = np.arctan2(w_r, u_r)
            # compute sideslip angle
        if self._Va == 0:
            self._beta = 0  # check this line
        else:
            self._beta = np.arcsin(v_r / self._Va)

    def _forces_moments(self, delta):
        """
        return the forces on the UAV based on the state, wind, and control surfaces
        :param delta: MsgDelta()
        :return: Forces and Moments on the UAV np.matrix(Fx, Fy, Fz, Ml, Mn, Mm)
        """
        phi, theta, psi = quaternion_2_euler(self._state[6:10])

        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)

        p = self._state.item(10)
        q = self._state.item(11)
        r = self._state.item(12)

        # compute gravitaional forces in body frame
        f_g_x = MAV.mass * MAV.gravity * 2 * (e1 * e3 - e2 * e0)
        f_g_y = MAV.mass * MAV.gravity * 2 * (e2 * e3 + e1 * e0)
        f_g_z = MAV.mass * MAV.gravity * (e3 ** 2 + e0 ** 2 - e1 ** 2 - e2 ** 2)

        # compute Lift and Drag coefficients
        # ---- calculate C_L, in page 63 ----
        # C_L = MAV.C_L_0 + MAV.C_L_alpha * self._alpha   # linear model
        term_1 = np.exp(-MAV.M * (self._alpha - MAV.alpha0))
        term_2 = np.exp(MAV.M * (self._alpha + MAV.alpha0))
        sigma = (1 + term_1 + term_2) / ((1 + term_1) * (1 + term_2))
        C_L = (1 - sigma) * (MAV.C_L_0 + MAV.C_L_alpha * self._alpha) + sigma * (
                2 * np.sign(self._alpha) * (np.sin(self._alpha) ** 2) * np.cos(self._alpha))

        # ---- calculate C_D, in page 63 ----
        # C_D = MAV.C_D_0 + MAV.C_D_alpha * self._alpha   # linear model

        # PAY ATTENTION: The e here is the Oswald efficiency factor, NOT the euler number.
        C_D = MAV.C_D_p + (MAV.C_L_0 + MAV.C_L_alpha * self._alpha) ** 2 / (np.pi * MAV.e_os * MAV.AR)

        # compute Lift and Drag Forces
        rVS_2 = MAV.rho * (self._Va ** 2) * MAV.S_wing / 2
        F_lift = rVS_2 * (C_L + MAV.C_L_q * MAV.c / (2 * self._Va) * q + MAV.C_L_delta_e * delta.elevator)
        F_drag = rVS_2 * (C_D + MAV.C_D_q * MAV.c / (2 * self._Va) * q + MAV.C_D_delta_e * delta.elevator)

        # compute propeller thrust and torque
        thrust_prop, torque_prop = self._motor_thrust_torque(self._Va, delta.throttle)

        # compute longitudinal forces in body frame
        c_alpha = np.cos(self._alpha)
        s_alpha = np.sin(self._alpha)
        fx = f_g_x + (c_alpha * -F_drag + s_alpha * F_lift) + thrust_prop
        fz = f_g_z + (s_alpha * - F_drag + c_alpha * -F_lift)

        # compute lateral forces in body frame
        fy = f_g_y + rVS_2 * (MAV.C_Y_0 + MAV.C_Y_beta * self._beta + MAV.C_Y_p * MAV.b / (2 * self._Va) * p +
                              MAV.C_Y_r * MAV.b / (2 * self._Va) * r + MAV.C_Y_delta_a * delta.aileron +
                              MAV.C_Y_delta_r * delta.rudder)

        # compute longitudinal torque in body frame
        My = rVS_2 * MAV.c * (MAV.C_m_0 + MAV.C_m_alpha * self._alpha +
                              MAV.C_m_q * MAV.c / (2 * self._Va) * q + MAV.C_m_delta_e * delta.elevator)
        # compute lateral torques in body frame
        Mx = rVS_2 * MAV.b * (MAV.C_ell_0 + MAV.C_ell_beta * self._beta + MAV.C_ell_p * MAV.b / (2 * self._Va) * p +
                              MAV.C_ell_r * MAV.b / (2 * self._Va) * r + MAV.C_ell_delta_a * delta.aileron +
                              MAV.C_ell_delta_r * delta.rudder) + torque_prop
        Mz = rVS_2 * MAV.b * (MAV.C_n_0 + MAV.C_n_beta * self._beta + MAV.C_n_p * MAV.b / (2 * self._Va) * p +
                              MAV.C_n_r * MAV.b / (2 * self._Va) * r + MAV.C_n_delta_a * delta.aileron +
                              MAV.C_n_delta_r * delta.rudder)

        self._forces[0] = fx
        self._forces[1] = fy
        self._forces[2] = fz
        return np.array([[fx, fy, fz, Mx, My, Mz]]).T

    def _motor_thrust_torque(self, Va, delta_t):
        # compute thrust and torque due to propeller  (See addendum by McLain, in page 8)
        # map delta_t throttle command(0 to 1) into motor input voltage, according to equation 4.6
        V_in = MAV.V_max * delta_t

        # Angular speed of propeller, according to equation 4.5
        a = (MAV.rho * (MAV.D_prop ** 5) / ((2 * np.pi) ** 2)) * MAV.C_Q0
        b = (MAV.rho * (MAV.D_prop ** 4) / (2 * np.pi)) * MAV.C_Q1 * Va + (MAV.K_Q ** 2) / MAV.R_motor
        c = MAV.rho * (MAV.D_prop ** 3) * MAV.C_Q2 * (Va ** 2) - MAV.K_Q / MAV.R_motor * V_in + MAV.K_Q * MAV.i0
        Omega_p = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)

        # thrust and torque due to propeller
        thrust_prop = ((MAV.rho * (MAV.D_prop ** 4) * MAV.C_T0) / (4 * (np.pi ** 2))) * (Omega_p ** 2) + \
                      ((MAV.rho * (MAV.D_prop ** 3) * MAV.C_T1 * Va) / (2 * np.pi)) * Omega_p + \
                      MAV.rho * (MAV.D_prop ** 2) * MAV.C_T2 * (Va ** 2)
        torque_prop = ((MAV.rho * (MAV.D_prop ** 5) * MAV.C_Q0) / (4 * (np.pi ** 2))) * (Omega_p ** 2) + \
                      ((MAV.rho * (MAV.D_prop ** 4) * MAV.C_Q1 * Va) / (2 * np.pi)) * Omega_p + \
                      MAV.rho * (MAV.D_prop ** 3) * MAV.C_Q2 * (Va ** 2)

        # The original motor model, deprecated in addendum.
        # thrust_prop = 1 / 2 * MAV.rho * MAV.S_prop * MAV.C_prop * ((MAV.k_motor * delta_t) ** 2 - Va ** 2)
        # torque_prop = - MAV.kTp * (MAV.kOmega * delta_t) ** 2

        return thrust_prop, torque_prop

    def _update_true_state(self):
        # update the class structure for the true state:
        #   [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        phi, theta, psi = quaternion_2_euler(self._state[6:10])
        pdot = quaternion_2_rotation(self._state[6:10]) @ self._state[3:6]
        self.true_state.north = self._state.item(0)
        self.true_state.east = self._state.item(1)
        self.true_state.altitude = -self._state.item(2)
        self.true_state.Va = self._Va
        self.true_state.alpha = self._alpha
        self.true_state.beta = self._beta
        self.true_state.phi = phi
        self.true_state.theta = theta
        self.true_state.psi = psi
        self.true_state.Vg = np.linalg.norm(pdot)
        self.true_state.gamma = np.arcsin(pdot.item(2) / self.true_state.Vg)
        self.true_state.chi = np.arctan2(pdot.item(1), pdot.item(0))
        self.true_state.p = self._state.item(10)
        self.true_state.q = self._state.item(11)
        self.true_state.r = self._state.item(12)
        self.true_state.wn = self._wind.item(0)
        self.true_state.we = self._wind.item(1)
