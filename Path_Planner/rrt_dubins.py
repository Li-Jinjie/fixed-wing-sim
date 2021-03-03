# rrt dubins path planner for mavsim_python
#
# mavsim_python
#     - Beard & McLain, PUP, 2012
#     - Last updated:
#         4/16/2019 - RWB
import numpy as np
from message_types.msg_waypoints import MsgWaypoints
from Path_Manager.draw_waypoints import DrawWaypoints
from Path_Planner.draw_map import DrawMap
from Path_Manager.dubins_parameters import DubinsParameters
import pyqtgraph as pg
import pyqtgraph.opengl as gl


# TODO: finish this

class RRTDubins:
    def __init__(self):
        self.segment_length = 300  # standard length of path segments
        self.plot_window = []
        self.plot_app = []
        self.dubins_path = DubinsParameters()

    def update(self, start_pose, end_pose, Va, world_map, radius):
        # generate tree
        tree = MsgWaypoints()
        tree.type = 'dubins'
        # add the start pose to the tree
        # add(ned_position, airspeed, course, cost, parent, flag_is_goal)
        tree.add(start_pose, Va, np.inf, 0, -1, 0)

        # check to see if start_pose connects directly to end_pose
        if distance(start_pose, end_pose) < self.segment_length and not collision(start_pose, end_pose, world_map):
            # no collision, connect directly
            tree.add(end_pose, Va, np.inf, distance(start_pose, end_pose), 0, 1)
        else:
            num_paths = 0
            while num_paths < 10:
                if self.extend_tree(tree, end_pose, Va, world_map) is True:
                    num_paths += 1

        # find path with minimum cost to end_node
        waypoints_not_smooth = find_minimum_path(tree, end_pose)
        waypoints = smooth_path(waypoints_not_smooth, world_map)

        self.plot_map(world_map, tree, waypoints_not_smooth, waypoints, radius)  # display all paths
        return waypoints

    def extend_tree(self, tree, end_pose, Va, world_map, radius):
        # refer to page 215 on uav_book
        flag_found = False
        while not flag_found:
            # extend tree by randomly selecting pose and extending tree toward that pose
            p = random_pose(world_map, pd=end_pose.item(2))

            # find closest configuration
            idx_min = 0
            len_min = np.inf
            for i in range(tree.num_waypoints):
                v = get_node(tree, i)
                len_tmp = np.linalg.norm(p - v)
                if len_tmp < len_min:
                    len_min = len_tmp
                    idx_min = i
            v_star = get_node(tree, idx_min)

            # plan path
            unit_vec = (p - v_star) / np.linalg.norm(p - v_star)
            unit_vec[2, 0] = 0  # no influence to pd
            v_plus = v_star + self.segment_length * unit_vec

            # check collision from v_star to v_plus
            if not collision(v_star, v_plus, world_map):
                course = np.arctan2((v_plus - v_star).item(1), (v_plus - v_star).item(0))
                # no collision, add this node
                tree.add(v_plus, Va, course, distance(v_star, v_plus), idx_min, 0)
                del course

                # if v_plus is added successfully, check collision from v_plus to end_point
                if not collision(v_plus, end_pose, world_map):
                    course = np.arctan2((end_pose - v_plus).item(1), (end_pose - v_plus).item(0))
                    # no collision, add end node directly
                    tree.add(end_pose, Va, course, distance(v_plus, end_pose), tree.num_waypoints - 1, 1)
                    del course

                    flag_found = True  # find a path!
                    break

        return flag_found

    def points_along_path(self, Del):
        # returns a list of points along the dubins path
        initialize_points = True
        # points along start circle
        theta_list = []
        if initialize_points:
            points = 
            initialize_points = False
        for angle in theta_list:
            points = 

        # points along straight line
        sig = 0
        while sig <= 1:
            points =
            sig = 

        # points along end circle
    
        theta_list = []
        for angle in theta_list:
            points = 
        return points

    def collision(self, start_pose, end_pose, world_map, radius):
        # check to see of path from start_pose to end_pose colliding with world_map
        collision_flag = 
        return collision_flag

    def plot_map(self, world_map, tree, waypoints, smoothed_waypoints, radius):
        scale = 4000
        # initialize Qt gui application and window
        self.plot_app = pg.QtGui.QApplication([])  # initialize QT
        self.plot_window = gl.GLViewWidget()  # initialize the view object
        self.plot_window.setWindowTitle('World Viewer')
        self.plot_window.setGeometry(0, 0, 1500, 1500)  # args: upper_left_x, upper_right_y, width, height
        grid = gl.GLGridItem() # make a grid to represent the ground
        grid.scale(scale/20, scale/20, scale/20) # set the size of the grid (distance between each line)
        self.plot_window.addItem(grid) # add grid to viewer
        self.plot_window.setCameraPosition(distance=scale, elevation=50, azimuth=-90)
        self.plot_window.setBackgroundColor('k')  # set background color to black
        self.plot_window.show()  # display configured window
        self.plot_window.raise_() # bring window to the front

        blue = np.array([[30, 144, 255, 255]])/255.
        red = np.array([[204, 0, 0]])/255.
        green = np.array([[0, 153, 51]])/255.
        DrawMap(world_map, self.plot_window)
        DrawWaypoints(waypoints, radius, blue, self.plot_window)
        DrawWaypoints(smoothed_waypoints, radius, red, self.plot_window)
        self.draw_tree(radius, green)
        # draw things to the screen
        self.plot_app.processEvents()

    def draw_tree(self, radius, color, window):
        R = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
        Del = 0.05
        for i in range(1, tree.num_waypoints):
            parent = int(tree.parent.item(i))
            self.dubins_path.update(column(tree.ned, parent), tree.course[parent],
                                    column(tree.ned, i), tree.course[i], radius)
            points = self.points_along_path(Del)
            points = points @ R.T
            tree_color = np.tile(color, (points.shape[0], 1))
            tree_plot_object = gl.GLLinePlotItem(pos=points,
                                                color=tree_color,
                                                width=2,
                                                antialias=True,
                                                mode='line_strip')
            self.plot_window.addItem(tree_plot_object)

def get_node(waypoints, ptr):
    if ptr > waypoints.num_waypoints - 1:
        raise ValueError("Index is out of the range of waypoints!")

    return waypoints.ned[:, ptr: ptr + 1]

def find_minimum_path(tree, end_pose):
    # obj: find the lowest cost path to the end node

    # find every end_node
    id_end_nodes = []
    for i in range(tree.num_waypoints):
        if tree.is_goal[i] == 1:
            id_end_nodes.append(i)

    # find minimum cost end node
    sum_cost_min = np.inf
    for id_end_node in id_end_nodes:
        sum_cost = 0
        sum_cost += tree.cost[id_end_node]
        idx = tree.parent[id_end_node]
        while idx != -1:  # not origin
            sum_cost += tree.cost[idx]
            idx = tree.parent[idx]

        if sum_cost < sum_cost_min:
            sum_cost_min = sum_cost
            idx_min = id_end_node

    # construct lowest cost path order
    path = []
    idx = idx_min
    while idx != -1:  # not origin's parent
        path.insert(0, idx)
        idx = tree.parent[idx]

    # construct waypoint path
    waypoints = MsgWaypoints()
    waypoints.type = tree.type
    parent = 0
    for idx in path:
        if idx == 0:
            waypoints.add(get_node(tree, idx), tree.airspeed[idx], np.inf, 0, -1, 0)  # first node, start pose
            continue
        ned_now = get_node(tree, idx)
        ned_pre = get_node(waypoints, parent)

        course = np.arctan2((ned_now - ned_pre).item(1), (ned_now - ned_pre).item(0))

        waypoints.add(ned_now, tree.airspeed[idx], course,
                      distance(ned_now, ned_pre), parent, tree.is_goal[idx])
        parent += 1

    return waypoints


def smooth_path(waypoints, world_map, radius):
    # smooth the waypoint path
    smooth = [0]  # add the first waypoint
    smooth_ptr = 0  # i
    origin_ptr = 1  # j

    while origin_ptr < waypoints.num_waypoints - 1:
        w_s = get_node(waypoints, smooth_ptr)
        w_plus = get_node(waypoints, origin_ptr + 1)
        if collision(w_s, w_plus, world_map):
            smooth.append(origin_ptr)  # add deconflicted node to smoothed path
            # add cost
            smooth_ptr = origin_ptr
        origin_ptr += 1
    smooth.append(origin_ptr)  # add end node

    # construct smooth waypoint path
    smooth_waypoints = MsgWaypoints()
    smooth_waypoints.type = waypoints.type
    parent = 0
    for idx in smooth:
        if idx == 0:
            smooth_waypoints.add(get_node(waypoints, idx), waypoints.airspeed[idx], np.inf, 0, -1,
                                 0)  # first node, start pose
            continue
        ned_now = get_node(waypoints, idx)
        ned_pre = get_node(smooth_waypoints, parent)

        course = np.arctan2((ned_now - ned_pre).item(1), (ned_now - ned_pre).item(0))

        smooth_waypoints.add(ned_now, waypoints.airspeed[idx], course,
                             distance(ned_now, ned_pre), parent, waypoints.is_goal[idx])
        parent += 1

    if smooth_waypoints.num_waypoints != waypoints.num_waypoints:
        smooth_waypoints = smooth_path(smooth_waypoints, world_map)  # recursion!

    return smooth_waypoints


def distance(start_pose, end_pose):
    # compute distance between start and end pose
    DubinsParameters()
    d = 
    return d


def collision(start_pose, end_pose, world_map):
    # check to see of path from start_pose to end_pose colliding with map
    N = np.floor(np.linalg.norm(end_pose - start_pose) / 10.)  # every 10 meters get a point
    points = points_along_path(start_pose, end_pose, N)

    redundancy = 10.  # meters, redundancy to prevent the uav from hitting the building.

    for point in points:
        for p_n in range(world_map.num_city_blocks):  # North
            for p_e in range(world_map.num_city_blocks):  # East
                safe_distance = world_map.building_width / 2. + redundancy
                if world_map.building_height[p_n, p_e] > np.abs(point.item(2)) and \
                        np.abs(point.item(0) - world_map.building_north.item(p_n)) <= safe_distance and \
                        np.abs(point.item(1) - world_map.building_north.item(p_e)) <= safe_distance:
                    collision_flag = True
                    point.item(0) - world_map.building_north.item(p_n)
                    return collision_flag

    collision_flag = False
    return collision_flag

# def height_above_ground(world_map, point):
#     # find the altitude of point above ground level
#     point_height =
#     if (d_n<world_map.building_width) and (d_e<world_map.building_width):
#         map_height =
#     else:
#         map_height =
#     h_agl =
#     return h_agl


def random_pose(world_map, pd):
    # generate a random pose
    pose = np.array([[world_map.city_width * np.random.random()],
                     [world_map.city_width * np.random.random()],
                     [pd]])
    return pose


def mod(x):
    # force x to be between 0 and 2*pi
    while x < 0:
        x += 2*np.pi
    while x > 2*np.pi:
        x -= 2*np.pi
    return x


def column(A, i):
    # extracts the ith column of A and return column vector
    tmp = A[:, i]
    col = tmp.reshape(A.shape[0], 1)
    return col