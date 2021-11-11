import numpy as np
# from math import sin, cos  TODO: change np.sin to math.sin
import sys

sys.path.append('..')
from message_types.msg_autopilot import MsgAutopilot
from tools.wrap import wrap
import parameters.aerosonde_parameters as MAV


class PathFollower:
    def __init__(self):
        self.chi_inf = np.radians(75.)  # approach angle for large distance from straight-line path
        self.k_path = 0.05  # 0.05  # proportional gain for straight-line path following
        self.k_orbit = 10.0  # 10.0  # proportional gain for orbit following
        self.gravity = MAV.gravity
        self.autopilot_commands = MsgAutopilot()  # message sent to autopilot

    def update(self, path, state):
        if path.type == 'line':
            self._follow_straight_line(path, state)
        elif path.type == 'orbit':
            self._follow_orbit(path, state)
        return self.autopilot_commands

    def _follow_straight_line(self, path, state):
        q = path.line_direction
        qn = q.item(0)
        qe = q.item(1)
        qd = q.item(2)
        r = path.line_origin
        rn = r.item(0)
        re = r.item(1)
        rd = r.item(2)
        pn = state.north
        pe = state.east
        pd = -state.altitude

        # airspeed command
        self.autopilot_commands.airspeed_command = path.airspeed

        # course command
        chi_q = wrap(np.arctan2(qe, qn), state.chi)
        e_py = - np.sin(chi_q) * (pn - rn) + np.cos(chi_q) * (pe - re)
        self.autopilot_commands.course_command = chi_q - self.chi_inf * (2. / np.pi) * np.arctan(self.k_path * e_py)

        # altitude command
        # page 177 of uav_book
        k = np.array([[0, 0, 1]]).T
        tmp = np.cross(q.T, k.T).T
        n = tmp / np.linalg.norm(tmp)
        ep_i = np.array([[pn], [pe], [pd]]) - r
        s_i = ep_i - np.dot(ep_i.T, n) * n
        sn = s_i.item(0)
        se = s_i.item(1)
        # refer to uav_book supplement
        self.autopilot_commands.altitude_command = - rd - np.sqrt(sn ** 2 + se ** 2) * (qd / np.sqrt(qn ** 2 + qe ** 2))
        # feedforward roll angle for straight line is zero
        self.autopilot_commands.phi_feedforward = 0.

    def _follow_orbit(self, path, state):
        if path.orbit_direction == 'CW':
            direction = 1.0
        else:
            direction = -1.0
        c = path.orbit_center
        cn = c.item(0)
        ce = c.item(1)
        cd = c.item(2)
        pn = state.north
        pe = state.east
        wn = state.wn
        we = state.we
        wd = 0.

        # airspeed command
        self.autopilot_commands.airspeed_command = path.airspeed
        # distance from orbit center
        d = np.sqrt((pn - cn) ** 2 + (pe - ce) ** 2)
        # compute wrapped version of angular position on orbit
        phi_var = wrap(np.arctan2(pe - ce, pn - cn), state.chi)
        # compute normalized orbit error
        orbit_error = (d - path.orbit_radius) / path.orbit_radius
        # course command
        self.autopilot_commands.course_command = phi_var + \
                                                 direction * (np.pi / 2. + np.arctan(self.k_orbit * orbit_error))
        # altitude command
        self.autopilot_commands.altitude_command = -cd
        # roll feedforward command
        if orbit_error < 5:  # 10 when the uav is very far from the circle, it will fly toward the center with phi=0.
            # please read the page 96 of uav_book supplement
            term_1 = path.airspeed ** 2
            term_2 = (wn * np.sin(state.chi) - we * np.cos(state.chi)) ** 2 + wd ** 2
            if term_1 >= term_2:
                Vg = wn * np.cos(state.chi) + we * np.sin(state.chi) + np.sqrt(term_1 - term_2)
                term_all = Vg ** 2 / (MAV.gravity * path.orbit_radius * np.sqrt((term_1 - term_2) / (term_1 - wd ** 2)))
                phi_ff = direction * np.arctan(term_all)
                self.autopilot_commands.phi_feedforward = phi_ff
            else:
                print("The wind is too strong to follow the path!")
                self.autopilot_commands.phi_feedforward = 0

        else:
            self.autopilot_commands.phi_feedforward = 0
