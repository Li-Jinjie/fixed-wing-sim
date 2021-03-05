"""
mavsim_python
    - Chapter 11 assignment for Beard & McLain, PUP, 2012
    - Last Update:
        3/26/2019 - RWB
        2/27/2020 - RWB
"""
import sys

sys.path.append('..')
import numpy as np
import parameters.simulation_parameters as SIM
import parameters.planner_parameters as PLAN

from Kinematics_and_Dynamics.data_viewer import DataViewer
from Forces_and_Moments.wind_simulation import WindSimulation
from Autopilot.autopilot import Autopilot
from Sensors_MAV.mav_dynamics import MavDynamics
from State_Estimation.observer import Observer
from Path_Follower.path_follower import PathFollower
from Path_Manager.path_manager import PathManager
from Path_Manager.waypoint_viewer import WaypointViewer

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
waypoint_view = WaypointViewer()  # initialize the viewer
data_view = DataViewer()  # initialize view of data plots
if VIDEO is True:
    from Coordinate_Frames.video_writer import VideoWriter

    video = VideoWriter(video_name="chap11_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)
autopilot = Autopilot(SIM.ts_simulation)
observer = Observer(SIM.ts_simulation)
path_follower = PathFollower()
path_manager = PathManager()

# waypoint definition
from message_types.msg_waypoints import MsgWaypoints

waypoints = MsgWaypoints()
# waypoints.type = 'straight_line'
# waypoints.type = 'fillet'
waypoints.type = 'dubins'
Va = PLAN.Va0
waypoints.add(np.array([[0, 0, -100]]).T, Va, np.radians(0), np.inf, 0, 0)
waypoints.add(np.array([[1000, 0, -100]]).T, Va, np.radians(45), np.inf, 0, 0)
waypoints.add(np.array([[0, 1000, -100]]).T, Va, np.radians(45), np.inf, 0, 0)
waypoints.add(np.array([[1000, 1000, -100]]).T, Va, np.radians(-135), np.inf, 0, 0)
waypoints.add(np.array([[0, 0, -100]]).T, Va, np.radians(0), np.inf, 0, 0)
waypoints.add(np.array([[1000, 0, -100]]).T, Va, np.radians(45), np.inf, 0, 0)
waypoints.add(np.array([[1000, 1000, -100]]).T, Va, np.radians(-135), np.inf, 0, 0)
waypoints.add(np.array([[0, 1000, -100]]).T, Va, np.radians(45), np.inf, 0, 0)

# # parabolic
# waypoints.add(np.array([[100, 0, -100 - 95/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[200, 0, -100 - 180/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[300, 0, -100 - 255/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[400, 0, -100 - 320/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[500, 0, -100 - 375/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[600, 0, -100 - 420/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[700, 0, -100 - 455/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[800, 0, -100 - 480/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[900, 0, -100 - 495/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[1000, 0, -100 - 500/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[1100, 0, -100 - 495/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[1200, 0, -100 - 480/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[1300, 0, -100 - 455/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[1400, 0, -100 - 420/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[1500, 0, -100 - 375/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[1600, 0, -100 - 320/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[1700, 0, -100 - 255/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[1800, 0, -100 - 180/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[1900, 0, -100 - 95/2]]).T, Va, np.radians(0), np.inf, 0, 0)
# waypoints.add(np.array([[2000, 0, -100 - 0]]).T, Va, np.radians(0), np.inf, 0, 0)

# initialize the simulation time
sim_time = SIM.start_time
plot_timer = 0

# main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:
    # -------observer-------------
    measurements = mav.sensors()  # get sensor measurements
    # estimated_state = observer.update(measurements)  # estimate states from measurements
    estimated_state = mav.true_state  # for debugging

    # -------path manager-------------
    path = path_manager.update(waypoints, PLAN.R_min, estimated_state)

    # -------path follower-------------
    autopilot_commands = path_follower.update(path, estimated_state)

    # -------autopilot-------------
    delta, commanded_state = autopilot.update(autopilot_commands, estimated_state)

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    # current_wind = np.zeros([6, 1])  # for debugging
    mav.update(delta, current_wind)  # propagate the MAV dynamics

    # -------update viewer-------------
    if plot_timer > SIM.ts_plotting:
        waypoint_view.update(mav.true_state, path, waypoints)  # plot path and MAV
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
