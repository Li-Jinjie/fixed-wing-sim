import numpy as np
import sys

sys.path.append('..')
from i_Path_Manager.dubins_parameters import DubinsParameters
from message_types.msg_path import MsgPath


# TODO: separate several path manager algorithm into different python files.
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
        self.dubins_path = DubinsParameters()
        # flag of end point
        self.flag_end_point = False

    def update(self, waypoints, radius, state):
        self.num_waypoints = waypoints.num_waypoints
        if waypoints.num_waypoints == 0:
            self.manager_requests_waypoints = True
        if self.manager_requests_waypoints is True \
                and waypoints.flag_waypoints_changed is True:
            self.manager_requests_waypoints = False
            self.flag_end_point = False
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
        self.path.plot_updated = False

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

            # if fly to the end node, planning again
            if self.flag_end_point:
                self.manager_requests_waypoints = True
                return

            # normal mode
            if self.manager_state == 1:
                self.construct_fillet_circle(waypoints, radius)
                self.manager_state = 2
            elif self.manager_state == 2:
                self.increment_pointers()
                self.construct_fillet_line(waypoints, radius)
                self.manager_state = 1

    def construct_fillet_line(self, waypoints, radius):
        flag_end_point = False

        w_previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous + 1]
        if self.ptr_current == 9999:
            w_current = None  # TODO: fix this
        else:
            w_current = waypoints.ned[:, self.ptr_current:self.ptr_current + 1]
        if self.ptr_next == waypoints.num_waypoints:  # the current point is the end point
            self.flag_end_point = True  # TODO: fix this
        else:
            w_next = waypoints.ned[:, self.ptr_next:self.ptr_next + 1]

        # update path variables

        q_previous = (w_current - w_previous) / np.linalg.norm(w_current - w_previous)
        self.path.type = 'line'
        self.path.line_origin = w_previous
        self.path.line_direction = q_previous
        self.path.plot_updated = False

        # update half plane variables
        # end_point, special condition
        if self.flag_end_point:
            self.halfspace_r = w_current
            self.halfspace_n = q_previous
            return

        q_current = (w_next - w_current) / np.linalg.norm(w_next - w_current)
        rhovar = np.arccos(-q_previous.T @ q_current)

        z = w_current - (radius / np.tan(rhovar / 2.)) * q_previous

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
        if direction == +1:
            self.path.orbit_direction = 'CW'
        else:
            self.path.orbit_direction = 'CCW'

        self.path.plot_updated = False

        # update half plane variables
        self.halfspace_r = w_current + (radius / np.tan(rho_var / 2.)) * q_current
        self.halfspace_n = q_current

    def dubins_manager(self, waypoints, radius, state):
        mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
        # if the waypoints have changed, update the waypoint pointer
        if waypoints.flag_waypoints_changed is True:
            self.initialize_pointers()
            self.construct_dubins_circle_start(waypoints, self.dubins_path, radius)
            self.manager_state = 1
            waypoints.flag_waypoints_changed = False

        # state machine for dubins path
        if self.inHalfSpace(mav_pos) is True:
            if self.manager_state == 1:
                self.halfspace_r = self.dubins_path.r1
                self.halfspace_n = self.dubins_path.n1
                self.manager_state = 2
            elif self.manager_state == 2:
                self.construct_dubins_line(waypoints, self.dubins_path)
                self.manager_state = 3
            elif self.manager_state == 3:
                self.construct_dubins_circle_end(waypoints, self.dubins_path, radius)
                self.manager_state = 4
            elif self.manager_state == 4:
                self.halfspace_r = self.dubins_path.r3
                self.halfspace_n = self.dubins_path.n3
                self.manager_state = 5
            elif self.manager_state == 5:
                self.increment_pointers()
                self.construct_dubins_circle_start(waypoints, self.dubins_path, radius)
                self.manager_state = 1

    def construct_dubins_circle_start(self, waypoints, dubins_path, radius):
        w_previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous + 1]
        if self.ptr_current == 9999:
            w_current = None  # TODO: fix this
        else:
            w_current = waypoints.ned[:, self.ptr_current:self.ptr_current + 1]

        # find Dubins Parameters
        dubins_path.update(w_previous, waypoints.course.item(self.ptr_previous),
                           w_current, waypoints.course.item(self.ptr_current), radius)

        # update path variables
        self.path.type = 'orbit'
        self.path.orbit_center = dubins_path.center_s
        self.path.orbit_radius = dubins_path.radius
        if dubins_path.dir_s == +1:
            self.path.orbit_direction = 'CW'
        else:
            self.path.orbit_direction = 'CCW'
        self.path.plot_updated = False

        # update half plane variables
        self.halfspace_r = dubins_path.r1
        self.halfspace_n = -dubins_path.n1

    def construct_dubins_line(self, waypoints, dubins_path):
        # update path variables
        self.path.type = 'line'
        self.path.line_origin = dubins_path.r1
        self.path.line_direction = dubins_path.n1
        self.path.plot_updated = False

        # update half plane variables
        self.halfspace_r = dubins_path.r2
        self.halfspace_n = dubins_path.n1

    def construct_dubins_circle_end(self, waypoints, dubins_path, radius):
        # update path variables
        self.path.type = 'orbit'
        self.path.orbit_center = dubins_path.center_e
        self.path.orbit_radius = dubins_path.radius
        if dubins_path.dir_e == +1:
            self.path.orbit_direction = 'CW'
        else:
            self.path.orbit_direction = 'CCW'
        self.path.plot_updated = False

        # update half plane variables
        self.halfspace_r = dubins_path.r3
        self.halfspace_n = -dubins_path.n3
