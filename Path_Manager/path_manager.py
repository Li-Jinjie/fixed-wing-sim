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
        self.halfspace_r = w_current
        q_previous = (w_current - w_previous) / np.linalg.norm(w_current - w_previous)
        q_current = (w_next - w_current) / np.linalg.norm(w_next - w_current)
        self.halfspace_n = (q_previous + q_current) / np.linalg.norm(q_previous + q_current)

        self.path.type = 'line'
        self.path.line_origin = w_previous
        self.path.line_direction = q_previous

    # def fillet_manager(self, waypoints, radius, state):
    #     mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
    #     # if the waypoints have changed, update the waypoint pointer
    #
    #     # state machine for fillet path
    #
    # def construct_fillet_line(self, waypoints, radius):
    #     previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous + 1]
    #     if self.ptr_current == 9999:
    #         current =
    #     else:
    #         current =
    #     if self.ptr_next == 9999:
    #         next =
    #     else:
    #         next =
    #         # update path variables
    #
    # def construct_fillet_circle(self, waypoints, radius):
    #     previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous + 1]
    #     if self.ptr_current == 9999:
    #         current =
    #     else:
    #         current =
    #     if self.ptr_next == 9999:
    #         next =
    #     else:
    #         next =
    #         # update path variables
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
