"""
mavsim_python
    - Chapter 12 assignment for Beard & McLain, PUP, 2012
    - Last Update:
        4/3/2019 - BGM
        2/27/2020 - RWB
"""
import sys

sys.path.append('..')
import numpy as np
import copy
import parameters.simulation_parameters as SIM
import parameters.planner_parameters as PLAN

from f_Sensors_MAV.mav_dynamics import MavDynamics
from c_Forces_and_Moments.wind_simulation import WindSimulation
from e_Autopilot.autopilot import Autopilot
from g_State_Estimation.observer import Observer
from h_Path_Follower.path_follower import PathFollower
from i_Path_Manager.path_manager import PathManager
from j_Path_Planner.path_planner import PathPlanner

from b_Kinematics_and_Dynamics.data_viewer import DataViewer
from j_Path_Planner.world_viewer import WorldViewer

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
world_view = WorldViewer()  # initialize the viewer
data_view = DataViewer()  # initialize view of data plots
if VIDEO is True:
    from a_Coordinate_Frames.video_writer import VideoWriter

    video = VideoWriter(video_name="chap12_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)
autopilot = Autopilot(SIM.ts_simulation)
initial_state = copy.deepcopy(mav.true_state)
observer = Observer(SIM.ts_simulation, initial_state)
path_follower = PathFollower()
path_manager = PathManager()
path_planner = PathPlanner()

from message_types.msg_world_map import MsgWorldMap

world_map = MsgWorldMap()

# initialize the simulation time
sim_time = SIM.start_time
plot_timer = 0

# main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:
    # -------observer-------------
    measurements = mav.sensors()  # get sensor measurements
    estimated_state = observer.update(measurements)  # estimate states from measurements

    # -------path planner - ----
    if path_manager.manager_requests_waypoints is True:
        # waypoints = path_planner.update(world_map, estimated_state, PLAN.R_min)
        waypoints = path_planner.update(world_map, mav.true_state, PLAN.R_min)  # for debugging

    # -------path manager-------------
    # path = path_manager.update(waypoints, PLAN.R_min, estimated_state)
    path = path_manager.update(waypoints, PLAN.R_min, mav.true_state)  # for debugging

    # -------path follower-------------
    # autopilot_commands = path_follower.update(path, estimated_state)
    autopilot_commands = path_follower.update(path, mav.true_state)  # for debugging

    # -------autopilot-------------
    # delta, commanded_state = autopilot.update(autopilot_commands, estimated_state)
    delta, commanded_state = autopilot.update(autopilot_commands, mav.true_state)  # for debugging

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    mav.update(delta, current_wind)  # propagate the MAV dynamics

    # -------update viewer-------------
    if plot_timer > SIM.ts_plotting:
        world_view.update(mav.true_state, path, waypoints, world_map)  # plot path and MAV
        data_view.update(mav.true_state,  # true states
                         estimated_state,  # estimated states
                         commanded_state,  # commanded states
                         delta,  # input to aircraft
                         SIM.ts_plotting)
        plot_timer = 0
    plot_timer += SIM.ts_simulation
    if VIDEO is True:
        video.update(sim_time)

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if VIDEO is True:
    video.close()
