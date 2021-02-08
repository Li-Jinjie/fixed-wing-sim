"""
pid_control
    - Beard & McLain, PUP, 2012
    - Last Update:
        2/6/2019 - RWB
"""
import sys
import numpy as np

sys.path.append('..')


class PidControl:
    def __init__(self, kp=0.0, ki=0.0, kd=0.0, Ts=0.01, sigma=0.05, limit=1.0):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.Ts = Ts
        self.limit = limit
        self.integrator = 0.0
        self.error_dot = 0.0
        self.y_dot = 0
        self.error_delay_1 = 0.0  # delayed by one time step
        # self.error_dot_delay_1 = 0.0
        # gains for differentiator
        self.a1 = (2.0 * sigma - Ts) / (2.0 * sigma + Ts)
        self.a2 = 2.0 / (2.0 * sigma + Ts)

    def update(self, y_ref, y, reset_flag=False):
        if reset_flag is True:
            self.integrator = 0
            self.error_dot = 0
            self.error_delay_1 = 0

        # compute the current error
        error = y_ref - y
        #  update the integrator using trapezoidal rule
        self.integrator = self.integrator + (self.Ts / 2) * (error + self.error_delay_1)
        #  update the differentiator using trapezoidal rule
        self.error_dot = self.a1 * self.error_dot + self.a2 * (error - self.error_delay_1)
        # update the delayed variable
        self.error_delay_1 = error
        # PID control
        u_unsat = self.kp * error + self.ki * self.integrator + self.kd * self.error_dot
        # saturate PID control at limit
        u_sat = self._saturate(u_unsat, self.limit)

        # implement integrator anti-windup
        if np.abs(self.ki) > 0.0001:  # not work if ki is too small
            self.integrator = self.integrator + 1.0 / self.ki * (u_sat - u_unsat)

        return u_sat

    def update_with_rate(self, y_ref, y, y_dot, reset_flag=False):
        # read the uavbook_supplement materials on https://uavbook.byu.edu/lib/exe/fetch.php?media=uavbook_supplement.pdf
        if reset_flag is True:
            self.integrator = 0
            self.error_delay_1 = 0

        # compute the current error
        error = y_ref - y
        # update the integrator using trapezoidal rule
        self.integrator = self.integrator + (self.Ts / 2) * (error + self.error_delay_1)
        # update the delayed variable
        self.error_delay_1 = error

        # PID control
        u_unsat = self.kp * error + self.ki * self.integrator - self.kd * y_dot   # error_dot = 0 - y_dot
        # saturate PID control at limit
        u_sat = self._saturate(u_unsat, self.limit)

        # implement integrator anti-windup
        if np.abs(self.ki) > 0.0001:  # not work if ki is too small
            self.integrator = self.integrator + 1.0 / self.ki * (u_sat - u_unsat)

        return u_sat

    def _saturate(self, u, limit):
        # saturate u at +- self.limit
        if u >= limit:
            u_sat = limit
        elif u <= -limit:
            u_sat = -limit
        else:
            u_sat = u
        return u_sat


class PiControl:
    def __init__(self, kp=0.0, ki=0.0, Ts=0.01, limit=1.0):
        self.kp = kp
        self.ki = ki
        self.Ts = Ts
        self.limit = limit
        self.integrator = 0.0
        self.error_delay_1 = 0.0

    def update(self, y_ref, y):

        # compute the current error
        error = y_ref - y
        #  update the integrator using trapezoidal rule
        self.integrator = self.integrator + (self.Ts / 2) * (error + self.error_delay_1)
        # update the delayed variable
        self.error_delay_1 = error
        # PI control
        u_unsat = self.kp * error + self.ki * self.integrator
        # saturate PID control at limit
        u_sat = self._saturate(u_unsat, self.limit)
        # implement integrator anti-windup
        if np.abs(self.ki) > 0.0001:
            self.integrator = self.integrator + 1.0 / self.ki * (u_sat - u_unsat)
        return u_sat

    def _saturate(self, u, limit):
        # saturate u at +- self.limit
        if u >= limit:
            u_sat = limit
        elif u <= -limit:
            u_sat = -limit
        else:
            u_sat = u
        return u_sat


class PdControlWithRate:
    # PD control with rate information
    # u = kp*(yref-y) - kd*ydot
    def __init__(self, kp=0.0, kd=0.0, limit=1.0):
        self.kp = kp
        self.kd = kd
        self.limit = limit
        self.error_delay_1 = 0

    def update(self, y_ref, y, y_dot):
        # read the uavbook_supplement materials on https://uavbook.byu.edu/lib/exe/fetch.php?media=uavbook_supplement.pdf

        # compute the current error
        error = y_ref - y
        # update the delayed variable
        self.error_delay_1 = error

        # PID control
        u_unsat = self.kp * error - self.kd * y_dot
        # saturate PID control at limit
        u_sat = self._saturate(u_unsat, self.limit)
        return u_sat

    def _saturate(self, u, limit):
        # saturate u at +- self.limit
        if u >= limit:
            u_sat = limit
        elif u <= -limit:
            u_sat = -limit
        else:
            u_sat = u
        return u_sat
