"""
observer
    - Beard & McLain, PUP, 2012
    - Last Update:
        3/4/2019 - RWB
"""
import sys
import numpy as np
from scipy import stats

sys.path.append('..')
import parameters.control_parameters as CTRL
import parameters.simulation_parameters as SIM
import parameters.sensor_parameters as SENSOR
import parameters.aerosonde_parameters as MAV
from tools.rotations import euler_2_rotation
from tools.wrap import wrap

from message_types.msg_state import MsgState

# TODO: finish this file in the future !!!

class Observer:
    def __init__(self, ts_control):
        # initialized estimated state message
        self.ekf = ekfFullState()

    def update(self, measurement):
        estimated_state = self.ekf.update(measurement)
        return estimated_state


class ekfFullState:
    # implement continous-discrete EKF to estimate full state
    def __init__(self):
        self.Q = np.diag([1e-10,  # 0 pn
                          1e-10,  # 1 pe
                          1e-1,  # 2 pd
                          1e-5,  # 3 u
                          5e-1,  # 4 v
                          5e-2,  # 5 w
                          1e-1,  # 6 phi
                          1e-10,  # 7 theta
                          1e-10,  # 8 psi
                          1e-10,  # 9 bx
                          1e-10,  # 10 by
                          1e-10,  # 11 bz
                          1e-10,  # 12 wn
                          1e-10])  # 13 we
        self.Q_gyro = np.array([[SENSOR.gyro_sigma ** 2, 0., 0.],
                                [0., SENSOR.gyro_sigma ** 2, 0.],
                                [0., 0., SENSOR.gyro_sigma ** 2]])
        self.Q_accel = np.array([[SENSOR.accel_sigma ** 2, 0., 0.],
                                 [0., SENSOR.accel_sigma ** 2, 0.],
                                 [0., 0., SENSOR.accel_sigma ** 2]])
        self.R_analog =
        self.R_gps = np.array([[SENSOR.gps_n_sigma ** 2, 0., 0., 0.],
                               [0., SENSOR.gps_e_sigma ** 2, 0., 0.],
                               [0., 0., SENSOR.gps_Vg_sigma ** 2, 0.],
                               [0., 0., 0., SENSOR.gps_course_sigma ** 2]])  # n, e, Vg, chi
        self.N = 4  # number of prediction step per sample
        self.Tp = (SIM.ts_control / self.N)
        self.xhat = np.array([[MAV.north0],  # (0)  pn
                              [MAV.east0],  # (1)  pe
                              [MAV.down0],  # (2)  pd
                              [MAV.u0],  # (3)  u
                              [MAV.v0],  # (4)  v
                              [MAV.w0],  # (5)  w
                              [MAV.phi0],  # (6)  phi
                              [MAV.theta0],  # (7)  theta
                              [MAV.psi0],  # (8)  psi
                              [0],  # (9)  bx
                              [0],  # (10)  by
                              [0],  # (11)  bz
                              [0],  # (12)  wn
                              [0]])  # (13)  we
        self.P = np.eye(14) * 0.01
        self.analog_threshold = stats.chi2.isf(q=0.01, df=)
        self.gps_n_old = 9999
        self.gps_e_old = 9999
        self.gps_Vg_old = 9999
        self.gps_course_old = 9999
        self.estimated_state = MsgState()

    def update(self, measurement):
        self.propagate_model(measurement)
        self.measurement_update(measurement)

        # write out estimate state
        self.estimated_state.pn = self.xhat.item(0)
        self.estimated_state.pe = self.xhat.item(1)
        self.estimated_state.h = -self.xhat.item(2)
        self.estimated_state.u = self.xhat.item(3)
        self.estimated_state.v = self.xhat.item(4)
        self.estimated_state.w = self.xhat.item(5)
        self.estimated_state.phi = self.xhat.item(6)
        self.estimated_state.theta = self.xhat.item(7)
        self.estimated_state.psi = self.xhat.item(8)
        self.estimated_state.bx = self.xhat.item(9)
        self.estimated_state.by = self.xhat.item(10)
        self.estimated_state.bz = self.xhat.item(11)
        self.estimated_state.wn = self.xhat.item(12)
        self.estimated_state.we = self.xhat.item(13)

        # estimate needed quantities that are not part of state
        R = euler_2_rotation(self.estimated_state.phi,
                             self.estimated_state.theta,
                             self.estimated_state.psi)
        vel_body = self.xhat[3:6]
        vel_world = R @ vel_body
        wind_world = np.vstack((self.xhat[12:], 0))
        wind_body = R.T @ wind_world

        # velocity vector relative to the airmass
        u_r = vel_body[0][0] - wind_body[0][0]
        v_r = vel_body[1][0] - wind_body[1][0]
        w_r = vel_body[2][0] - wind_body[2][0]

        # compute airspeed
        self.estimated_state.Va = np.sqrt(u_r ** 2 + v_r ** 2 + w_r ** 2)
        # compute angle of attack
        if u_r == 0:
            self.estimated_state.alpha = np.pi / 2  # check this line
        else:
            self.estimated_state.alpha = np.arctan2(w_r, u_r)

        # compute sideslip angle
        if self.estimated_state.Va == 0:
            self.estimated_state.beta = 0  # check this line
        else:
            self.estimated_state.beta = np.arcsin(v_r / self.estimated_state.Va)

        self.estimated_state.Vg = np.linalg.norm(vel_world)
        self.estimated_state.chi = np.arctan2(vel_world.item(1), vel_world.item(0))
        self.estimated_state.p = measurement.gyro_x - self.estimated_state.bx
        self.estimated_state.q = measurement.gyro_y - self.estimated_state.by
        self.estimated_state.r = measurement.gyro_z - self.estimated_state.bz
        return self.estimated_state

    def propagate_model(self, measurement):
        # model propagation
        for i in range(0, self.N):
            state = self.estimated_state
            phi = state.phi
            theta = state.theta
            psi = state.psi
            S = np.array([[1, np.sin(phi) * np.tan(theta), np.cos(phi) * np.tan(theta)],
                          [0, np.cos(phi), -np.sin(phi)],
                          [0, np.sin(phi) / np.cos(theta), np.cos(phi) / np.cos(theta)]])
            [[u], [v], [w]] = state.Va * np.array([[np.cos(state.alpha) * np.cos(state.beta)],
                                                   [np.sin(state.beta)],
                                                   [np.sin(state.alpha) * np.cos(state.beta)]])
            vel_cross = np.array([[0, -w, v],
                                  [w, 0, -u],
                                  [-v, u, 0]])
            # propagate model
            Tp = self.Tp
            self.xhat = self.xhat + Tp * self.f(self.xhat, measurement, state)
            # compute Jacobian
            A = jacobian(self.f, self.xhat, measurement, state)
            # convert to discrete time models
            A_d = np.eye(A.shape[0]) + A * Tp + A @ A * (Tp ** 2)
            Gg_d =
            Ga_d =
            # update P with discrete time model
            self.P =

    def f(self, x, measurement):
        # system dynamics for propagation model: xdot = f(x, u)
        pos_dot =
        vel_dot =
        Theta_dot =
        bias_dot =
        wind_dot =
        f_ = np.concatenate((pos_dot, vel_dot, Theta_dot, bias_dot, wind_dot), axis=0)
        return f_

    def measurement_update(self, measurement):
        # always update based on sensor measurements
        h = self.h_analog(self.xhat, measurement)
        C = jacobian(self.h_analog, self.xhat, measurement)
        y = np.array([[measurement.abs_pressure,
                       measurement.diff_pressure,
                       0.0]]).T
        S_inv =
        if False:  # (y - h).T @ S_inv @ (y - h) < self.analog_threshold:
            L =
            self.P =
            self.xhat =

            # only update GPS when one of the signals changes
        if (measurement.gps_n != self.gps_n_old) \
                or (measurement.gps_e != self.gps_e_old) \
                or (measurement.gps_Vg != self.gps_Vg_old) \
                or (measurement.gps_course != self.gps_course_old):

            h =  # self.h_gps()
            C =  # jacobian()
            y_chi =  # wrap()
            y = np.array([[measurement.gps_n, measurement.gps_e, measurement.gps_Vg, y_chi]]).T
            if np.linalg.norm(y - h) < 0.5:
                L =
                self.P =
                self.xhat =

                # update stored GPS signals
            self.gps_n_old = measurement.gps_n
            self.gps_e_old = measurement.gps_e
            self.gps_Vg_old = measurement.gps_Vg
            self.gps_course_old = measurement.gps_course

    def h_analog(self, x, measurement):
        # analog sensor measurements and pseudo measurements
        abs_pres =
        diff_pres =
        sideslip =
        h = np.array([[abs_pres, diff_pres, sideslip]]).T
        return h

    def h_gps(self, x, measurement):
        # measurement model for gps measurements
        north =
        east =
        Vg =
        chi =
        h = np.array([[north, east, Vg, chi]]).T
        return h


def jacobian(fun, x, measurement):
    # compute jacobian of fun with respect to x
    f = fun(x, measurement)
    m = f.shape[0]
    n = x.shape[0]
    eps = 0.01  # deviation
    J = np.zeros((m, n))
    for i in range(0, n):
        x_eps = np.copy(x)
        x_eps[i][0] += eps
        f_eps = fun(x_eps, measurement)
        df = (f_eps - f) / eps
        J[:, i] = df[:, 0]
    return J
