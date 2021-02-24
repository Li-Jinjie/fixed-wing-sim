import numpy as np
import sys

sys.path.append('..')
# from Path_Manager.dubins_parameters import DubinsParameters
from message_types.msg_path import MsgPath


class PathManager:
    def __init__(self):
        # message sent to path follower
        self.path = MsgPath()
        # pointers to previous, current, and next waypoints
        self.ptr_previous = 0
        self.ptr_current = 1
        self.ptr_next = 2
        self.num_waypoints = 0
        self.halfspace_n = np.inf * np.ones((3, 1))
        self.halfspace_r = np.inf * np.ones((3, 1))
        # state of the manager state machine
        self.manager_state = 1
        self.manager_requests_waypoints = True
        # self.dubins_path = DubinsParameters()

    def update(self, waypoints, radius, state):
        if waypoints.num_waypoints == 0:
            self.manager_requests_waypoints = True
        if self.manager_requests_waypoints is True \
                and waypoints.flag_waypoints_changed is True:
            self.manager_requests_waypoints = False
        if waypoints.type == 'straight_line':
            self.line_manager(waypoints, state)
        elif waypoints.type == 'fillet':
            self.fillet_manager(waypoints, radius, state)
        elif waypoints.type == 'dubins':
            self.dubins_manager(waypoints, radius, state)
        else:
            print('Error in Path Manager: Undefined waypoint type.')
        return self.path

    def initialize_pointers(self):
        if self.num_waypoints >= 3:
            self.ptr_previous = 0
            self.ptr_current = 1
            self.ptr_next = 2
        else:
            print('Error Path Manager: need at least three waypoints')

    def increment_pointers(self):
        self.ptr_previous += 1
        self.ptr_current += 1
        self.ptr_next += 1

    def inHalfSpace(self, pos):
        if (pos - self.halfspace_r).T @ self.halfspace_n >= 0:
            return True
        else:
            return False

    def line_manager(self, waypoints, state):
        mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
        # if the waypoints have changed, update the waypoint pointer
        if waypoints.flag_waypoints_changed is True:
            self.initialize_pointers()
            self.construct_line(waypoints)
            waypoints.flag_waypoints_changed = False

        # state machine for line path
        if self.inHalfSpace(mav_pos) is True:
            self.increment_pointers()
            self.construct_line(waypoints)

    def construct_line(self, waypoints):
        w_previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous + 1]
        if self.ptr_current == 9999:
            w_current = None  # TODO: fix this
        else:
            w_current = waypoints.ned[:, self.ptr_current:self.ptr_current + 1]
        if self.ptr_next == 9999:
            w_next = None  # TODO: fix this
        else:
            w_next = waypoints.ned[:, self.ptr_next:self.ptr_next + 1]
        # update path variables
        # refer to page 189, algorithm 5, line 4-8
        q_previous = (w_current - w_previous) / np.linalg.norm(w_current - w_previous)
        q_current = (w_next - w_current) / np.linalg.norm(w_next - w_current)

        self.path.type = 'line'
        self.path.line_origin = w_previous
        self.path.line_direction = q_previous

        # update half plane variables
        self.halfspace_r = w_current
        self.halfspace_n = (q_previous + q_current) / np.linalg.norm(q_previous + q_current)

    def fillet_manager(self, waypoints, radius, state):
        mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
        # if the waypoints have changed, update the waypoint pointer
        if waypoints.flag_waypoints_changed is True:
            self.initialize_pointers()
            self.manager_state = 1
            self.construct_fillet_line(waypoints, radius)
            waypoints.flag_waypoints_changed = False

        # state machine for fillet path
        if self.inHalfSpace(mav_pos) is True:
            if self.manager_state == 1:
                self.construct_fillet_circle(waypoints, radius)
                self.manager_state = 2
            elif self.manager_state == 2:
                self.increment_pointers()
                self.construct_fillet_line(waypoints, radius)
                self.manager_state = 1

    def construct_fillet_line(self, waypoints, radius):
        w_previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous + 1]
        if self.ptr_current == 9999:
            w_current = None  # TODO: fix this
        else:
            w_current = waypoints.ned[:, self.ptr_current:self.ptr_current + 1]
        if self.ptr_next == 9999:
            w_next = None  # TODO: fix this
        else:
            w_next = waypoints.ned[:, self.ptr_next:self.ptr_next + 1]

        # update path variables
        q_previous = (w_current - w_previous) / np.linalg.norm(w_current - w_previous)
        q_current = (w_next - w_current) / np.linalg.norm(w_next - w_current)
        rhovar = np.arccos(-q_previous.T @ q_current)

        z = w_current - (radius / np.tan(rhovar / 2.)) * q_previous

        self.path.type = 'line'
        self.path.line_origin = w_previous
        self.path.line_direction = q_previous

        # update half plane variables
        self.halfspace_r = z
        self.halfspace_n = q_previous

    def construct_fillet_circle(self, waypoints, radius):
        w_previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous + 1]
        if self.ptr_current == 9999:
            w_current = None  # TODO: fix this
        else:
            w_current = waypoints.ned[:, self.ptr_current:self.ptr_current + 1]
        if self.ptr_next == 9999:
            w_next = None  # TODO: fix this
        else:
            w_next = waypoints.ned[:, self.ptr_next:self.ptr_next + 1]

        # update path variables
        q_previous = (w_current - w_previous) / np.linalg.norm(w_current - w_previous)
        q_current = (w_next - w_current) / np.linalg.norm(w_next - w_current)
        rho_var = np.arccos(-q_previous.T @ q_current)

        self.path.type = 'orbit'
        self.path.orbit_center = w_current - (radius / np.sin(rho_var / 2.)) * \
                                 (q_previous - q_current) / np.linalg.norm(q_previous - q_current)

        self.path.orbit_radius = radius
        direction = np.sign(q_previous.item(0) * q_current.item(1)
                            - q_previous.item(1) * q_current.item(0))
        if direction == 1:
            self.path.orbit_direction = 'CW'
        else:
            self.path.orbit_direction = 'CCW'

        # update half plane variables
        self.halfspace_r = w_current + (radius / np.tan(rho_var / 2.)) * q_current
        self.halfspace_n = q_current

#
# def dubins_manager(self, waypoints, radius, state):
#     mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
#     # if the waypoints have changed, update the waypoint pointer
#
#     # state machine for dubins path
#
# def construct_dubins_circle_start(self, waypoints, dubins_path):
#     pass
#     # update path variables
#
# def construct_dubins_line(self, waypoints, dubins_path):
#     pass
#     # update path variables
#
# def construct_dubins_circle_end(self, waypoints, dubins_path):
#     pass
#     # update path variables
