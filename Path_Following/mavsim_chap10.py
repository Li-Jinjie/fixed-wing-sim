"""
mavsim_python
    - Chapter 10 assignment for Beard & McLain, PUP, 2012
    - Last Update:
        3/11/2019 - RWB
        2/27/2020 - RWB
"""
import sys

sys.path.append('..')
import numpy as np
import parameters.simulation_parameters as SIM

from Kinematics_and_Dynamics.data_viewer import DataViewer
from Forces_and_Moments.wind_simulation import WindSimulation
from Autopilot.autopilot import Autopilot
from Sensors_MAV.mav_dynamics import MavDynamics
from State_Estimation.observer import Observer
from Path_Following.path_follower import PathFollower
from Path_Following.path_viewer import PathViewer

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
path_view = PathViewer()  # initialize the viewer
data_view = DataViewer()  # initialize view of data plots
if VIDEO is True:
    from Coordinate_Frames.video_writer import VideoWriter

    video = VideoWriter(video_name="chap10_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)
autopilot = Autopilot(SIM.ts_simulation)
observer = Observer(SIM.ts_simulation)
path_follower = PathFollower()

# path definition
from message_types.msg_path import MsgPath

path = MsgPath()
# path.type = 'line'
path.type = 'orbit'
if path.type == 'line':
    path.line_origin = np.array([[0.0, 0.0, -100.0]]).T
    path.line_direction = np.array([[0.5, 1.0, -0.1]]).T
    path.line_direction = path.line_direction / np.linalg.norm(path.line_direction)
elif path.type == 'orbit':
    path.orbit_center = np.array([[0.0, 0.0, -100.0]]).T  # center of the orbit
    path.orbit_radius = 300.0  # radius of the orbit
    path.orbit_direction = 'CW'  # orbit direction: 'CW'==clockwise, 'CCW'==counter clockwise

# initialize the simulation time
sim_time = SIM.start_time
plot_timer = 0

# main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:
    # -------observer-------------
    measurements = mav.sensors()  # get sensor measurements
    estimated_state = observer.update(measurements)  # estimate states from measurements

    # -------path follower-------------
    # autopilot_commands = path_follower.update(path, estimated_state)
    autopilot_commands = path_follower.update(path, mav.true_state)  # for debugging

    # -------autopilot-------------
    # delta, commanded_state = autopilot.update(autopilot_commands, estimated_state)
    delta, commanded_state = autopilot.update(autopilot_commands, mav.true_state)  # for debugging

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    # current_wind = np.zeros([6, 1])   # for debugging
    mav.update(delta, current_wind)  # propagate the MAV dynamics

    # -------update viewer-------------
    if plot_timer > SIM.ts_plotting:
        path_view.update(mav.true_state, path)  # plot path and MAV
        data_view.update(mav.true_state,  # true states
                         estimated_state,  # estimated states
                         commanded_state,  # commanded states
                         delta,  # input to aircraft
                         SIM.ts_simulation)
        plot_timer = 0
    plot_timer += SIM.ts_simulation

    if VIDEO is True:
        video.update(sim_time)

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if VIDEO is True:
    video.close()
