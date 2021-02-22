"""
observer
    - Beard & McLain, PUP, 2012
    - Last Update:
        3/2/2019 - RWB
        refer to: https://github.com/eyler94/EE674LQR/blob/master/chap8/observer.py，
                  https://github.com/betaBison/EC-EN-674-Flight-Dynamics-Controls/blob/master/chap8/observer_ekf.py
"""
import sys
import numpy as np
from scipy import stats

sys.path.append('..')
import parameters.control_parameters as CTRL
import parameters.simulation_parameters as SIM
import parameters.sensor_parameters as SENSOR
import parameters.aerosonde_parameters as MAV
from tools.wrap import wrap
from message_types.msg_state import MsgState


class Observer:
    def __init__(self, ts_control):
        # initialized estimated state message
        self.estimated_state = MsgState()
        # use alpha filters to low pass filter gyros and accels
        self.lpf_gyro_x = AlphaFilter(alpha=0.5)
        self.lpf_gyro_y = AlphaFilter(alpha=0.2)
        self.lpf_gyro_z = AlphaFilter(alpha=0.8)
        self.lpf_accel_x = AlphaFilter(alpha=0.5)
        self.lpf_accel_y = AlphaFilter(alpha=0.5)
        self.lpf_accel_z = AlphaFilter(alpha=0.5)
        # use alpha filters to low pass filter absolute and differential pressure
        self.lpf_abs = AlphaFilter(alpha=0.9, y0=MAV.rho * MAV.gravity * -MAV.down0)
        self.lpf_diff = AlphaFilter(alpha=0.5, y0=(MAV.rho * (MAV.Va0 ** 2.)) / 2.)
        # ekf for phi and theta
        self.attitude_ekf = EkfAttitude()
        # ekf for pn, pe, Vg, chi, wn, we, psi
        self.position_ekf = EkfPosition()

    def update(self, measurement):
        # estimates for p, q, r are low pass filter of gyro minus bias estimate
        self.estimated_state.p = self.lpf_gyro_x.update(measurement.gyro_x) - SENSOR.gyro_x_bias
        self.estimated_state.q = self.lpf_gyro_y.update(measurement.gyro_y) - SENSOR.gyro_y_bias
        self.estimated_state.r = self.lpf_gyro_z.update(measurement.gyro_z) - SENSOR.gyro_z_bias

        # invert sensor model to get altitude and airspeed
        self.estimated_state.altitude = self.lpf_abs.update(measurement.abs_pressure) / (MAV.rho * MAV.gravity)
        self.estimated_state.Va = np.sqrt(2. * self.lpf_diff.update(measurement.diff_pressure) / MAV.rho)

        # estimate phi and theta with simple ekf
        self.attitude_ekf.update(measurement, self.estimated_state)

        # estimate pn, pe, Vg, chi, wn, we, psi
        self.position_ekf.update(measurement, self.estimated_state)

        # not estimating these
        self.estimated_state.alpha = self.estimated_state.theta
        self.estimated_state.beta = 0.0
        self.estimated_state.bx = 0.0
        self.estimated_state.by = 0.0
        self.estimated_state.bz = 0.0
        return self.estimated_state


class AlphaFilter:
    # alpha filter implements a simple low pass filter
    # y[k] = alpha * y[k-1] + (1-alpha) * u[k]
    def __init__(self, alpha=0.5, y0=0.0):
        self.alpha = alpha  # filter parameter
        self.y = y0  # initial condition

    def update(self, u):
        self.y = self.alpha * self.y + (1 - self.alpha) * u
        return self.y


class EkfAttitude:
    # implement continous-discrete EKF to estimate roll and pitch angles
    # Refer to "8.11.4 Algorithm for Direct Kalman Filter" on the supplement materials, page 90.

    def __init__(self):
        self.Q = 1e-10 * np.eye(2)  # represent modelling error
        self.Q_gyro = np.array([[SENSOR.gyro_sigma ** 2, 0., 0.],
                                [0., SENSOR.gyro_sigma ** 2, 0.],
                                [0., 0., SENSOR.gyro_sigma ** 2]])
        self.R_accel = np.array([[SENSOR.accel_sigma ** 2, 0., 0.],
                                 [0., SENSOR.accel_sigma ** 2, 0.],
                                 [0., 0., SENSOR.accel_sigma ** 2]])
        self.N = 4  # number of prediction step per sample
        self.xhat = np.array([[MAV.phi0],
                              [MAV.theta0]])  # initial state: phi, theta
        self.P = np.eye(2) * 0.1  # Covariance of the estimation error
        self.Tp = (SIM.ts_control / self.N)
        self.gate_threshold = stats.chi2.isf(q=0.01, df=3)  # 99%; measurement equations are 3-dimension

    def update(self, measurement, state):
        self.propagate_model(measurement, state)
        self.measurement_update(measurement, state)
        state.phi = self.xhat.item(0)
        state.theta = self.xhat.item(1)

    def f(self, x, measurement, state):
        # system dynamics for propagation model: xdot = f(x, u)
        # refer to page 157 of uavbook
        p = state.p  # u
        q = state.q
        r = state.r
        phi = x.item(0)
        theta = x.item(1)

        G = np.array([[1., np.sin(phi) * np.tan(theta), np.cos(phi) * np.tan(theta)],
                      [0., np.cos(phi), -np.sin(phi)]])
        f_ = G @ np.array([[p, q, r]]).T
        return f_

    def h(self, x, measurement, state):
        # measurement model y
        # refer to page 157
        p = state.p  # u
        q = state.q
        r = state.r
        Va = state.Va
        phi = x.item(0)
        theta = x.item(1)

        h_ = np.array([[q * Va * np.sin(theta) + MAV.gravity * np.sin(theta)],
                       [r * Va * np.cos(theta) - p * Va * np.sin(theta) - MAV.gravity * np.cos(theta) * np.sin(phi)],
                       [-q * Va * np.cos(theta) - MAV.gravity * np.cos(theta) * np.cos(phi)]])
        return h_

    def propagate_model(self, measurement, state):
        # model propagation
        # refer to page 73 of uavbook_supplement
        for i in range(0, self.N):
            # propagate model
            Tp = self.Tp
            self.xhat = self.xhat + Tp * self.f(self.xhat, measurement, state)
            # compute Jacobian
            A = jacobian(self.f, self.xhat, measurement, state)
            # compute G matrix for gyro noise
            phi = self.xhat.item(0)
            theta = self.xhat.item(1)
            G = np.array([[1., np.sin(phi) * np.tan(theta), np.cos(phi) * np.tan(theta)],
                          [0., np.cos(phi), -np.sin(phi)]])
            # convert to discrete time models
            A_d = np.eye(A.shape[0]) + A * Tp + A @ A * (Tp ** 2)
            # update P with discrete time model
            self.P = A_d @ self.P @ A_d.T + (Tp ** 2) * (G @ self.Q_gyro @ G.T + self.Q)

    def measurement_update(self, measurement, state):
        # measurement updates
        h = self.h(self.xhat, measurement, state)
        C = jacobian(self.h, self.xhat, measurement, state)
        y = np.array([[measurement.accel_x, measurement.accel_y, measurement.accel_z]]).T
        S_inv = np.linalg.inv(self.R_accel + C @ self.P @ C.T)
        if ((y - h).T @ S_inv @ (y - h)) < self.gate_threshold:  # refer to page 74, 8.8 in uav_book supplement
            L = self.P @ C.T @ S_inv
            tmp = np.eye(L.shape[0]) - L @ C
            self.P = tmp @ self.P @ tmp.T + L @ self.R_accel @ L.T
            self.xhat = self.xhat + L @ (y - h)
            # print('updating')


class EkfPosition:
    # implement continous-discrete EKF to estimate pn, pe, Vg, chi, wn, we, psi
    def __init__(self):
        self.Q = np.diag([1e-10,
                          1e-10,
                          1e-1,
                          1e-5,
                          5e-1,
                          5e-2,
                          1e-1])
        self.R_gps = np.array([[SENSOR.gps_n_sigma ** 2, 0., 0., 0.],
                               [0., SENSOR.gps_e_sigma ** 2, 0., 0.],
                               [0., 0., SENSOR.gps_Vg_sigma ** 2, 0.],
                               [0., 0., 0., SENSOR.gps_course_sigma ** 2]])  # n, e, Vg, chi
        self.R_pseudo = np.diag([0.01, 0.01])  #
        self.N = 4  # number of prediction step per sample
        self.Tp = (SIM.ts_control / self.N)
        self.xhat = np.array([[MAV.north0],
                              [MAV.east0],
                              [MAV.Va0],
                              [0.],
                              [0.],
                              [0.],
                              [0.]])
        self.P = np.eye(7) * 0.1
        self.gps_n_old = 9999
        self.gps_e_old = 9999
        self.gps_Vg_old = 9999
        self.gps_course_old = 9999
        self.pseudo_threshold = stats.chi2.isf(q=0.01, df=2)
        self.gps_threshold = np.inf  # don't gate GPS

    def update(self, measurement, state):
        self.propagate_model(measurement, state)
        self.measurement_update(measurement, state)
        state.north = self.xhat.item(0)
        state.east = self.xhat.item(1)
        state.Vg = self.xhat.item(2)
        state.chi = self.xhat.item(3)
        state.wn = self.xhat.item(4)
        state.we = self.xhat.item(5)
        state.psi = self.xhat.item(6)

    def f(self, x, measurement, state):
        # system dynamics for propagation model: xdot = f(x, u)
        # refer to page 159 of uav_book
        Vg = x.item(2)
        chi = x.item(3)
        wn = x.item(4)
        we = x.item(5)
        psi = x.item(6)
        psidot = state.q * np.sin(state.phi) / np.cos(state.theta) + state.r * np.cos(state.phi) / np.cos(state.theta)
        # refer to page 75 of uav_book supplement
        Vgdot = state.Va * psidot * (we * np.cos(psi) - wn * np.sin(psi)) / Vg
        # Vgdot = ((state.Va * np.cos(psi) + wn) * (-state.Va * psidot * np.sin(psi)) +
        #          (state.Va * np.sin(psi) + we) * (state.Va * psidot * np.cos(psi))) / Vg  # Original version
        f_ = np.array([[Vg * np.cos(chi)],
                       [Vg * np.sin(chi)],
                       [Vgdot],
                       [(MAV.gravity / Vg) * np.tan(state.phi) * np.cos(chi - psi)],
                       [0],
                       [0],
                       [psidot]])
        return f_

    def h_gps(self, x, measurement, state):
        # measurement model for gps measurements
        pn = x.item(0)
        pe = x.item(1)
        Vg = x.item(2)
        chi = x.item(3)
        h_ = np.array([[pn],
                       [pe],
                       [Vg],
                       [chi]])
        return h_

    def h_pseudo(self, x, measurement, state):
        # measurement model for wind triangale pseudo measurement
        Vg = x.item(2)
        chi = x.item(3)
        wn = x.item(4)
        we = x.item(5)
        psi = x.item(6)
        h_ = np.array([[state.Va * np.cos(psi) + wn - Vg * np.cos(chi)],
                       [state.Va * np.sin(psi) + we - Vg * np.sin(chi)]])
        return h_

    def propagate_model(self, measurement, state):
        # model propagation
        for i in range(0, self.N):
            # propagate model
            Tp = self.Tp
            self.xhat = self.xhat + Tp * self.f(self.xhat, measurement, state)
            # compute Jacobian
            A = jacobian(self.f, self.xhat, measurement, state)
            # convert to discrete time models
            A_d = np.eye(A.shape[0]) + A * Tp + A @ A * (Tp ** 2)
            # update P with discrete time model
            self.P = A_d @ self.P @ A_d.T + (Tp ** 2) * self.Q

    def measurement_update(self, measurement, state):
        # always update based on wind triangle pseudu measurement
        h = self.h_pseudo(self.xhat, measurement, state)
        C = jacobian(self.h_pseudo, self.xhat, measurement, state)
        y = np.array([[0, 0]]).T
        S_inv = np.linalg.inv(self.R_pseudo + C @ self.P @ C.T)
        if (y - h).T @ S_inv @ (y - h) < self.pseudo_threshold:
            L = self.P @ C.T @ S_inv
            tmp = np.eye(L.shape[0]) - L @ C
            self.P = tmp @ self.P @ tmp.T + L @ self.R_pseudo @ L.T
            self.xhat = self.xhat + L @ (y - h)

            # only update GPS when one of the signals changes
        if (measurement.gps_n != self.gps_n_old) \
                or (measurement.gps_e != self.gps_e_old) \
                or (measurement.gps_Vg != self.gps_Vg_old) \
                or (measurement.gps_course != self.gps_course_old):

            h = self.h_gps(self.xhat, measurement, state)
            C = jacobian(self.h_gps, self.xhat, measurement, state)
            y_chi = wrap(measurement.gps_course, h[3, 0])
            y = np.array([[measurement.gps_n,
                           measurement.gps_e,
                           measurement.gps_Vg,
                           y_chi]]).T
            S_inv = np.linalg.inv(self.R_gps + C @ self.P @ C.T)
            if (y - h).T @ S_inv @ (y - h) < self.gps_threshold:
                L = self.P @ C.T @ S_inv
                tmp = np.eye(L.shape[0]) - L @ C
                self.P = tmp @ self.P @ tmp.T + L @ self.R_gps @ L.T
                self.xhat = self.xhat + L @ (y - h)

                # update stored GPS signals
            self.gps_n_old = measurement.gps_n
            self.gps_e_old = measurement.gps_e
            self.gps_Vg_old = measurement.gps_Vg
            self.gps_course_old = measurement.gps_course


def jacobian(fun, x, measurement, state):
    # compute jacobian of fun with respect to x
    f = fun(x, measurement, state)
    m = f.shape[0]
    n = x.shape[0]
    eps = 0.0001  # deviation
    J = np.zeros((m, n))
    for i in range(0, n):
        x_eps = np.copy(x)
        x_eps[i][0] += eps
        f_eps = fun(x_eps, measurement, state)
        df = (f_eps - f) / eps
        J[:, i] = df[:, 0]
    return J
