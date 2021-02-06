"""
mavsimPy
    - Chapter 4 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        12/27/2018 - RWB
        1/17/2019 - RWB
"""
import sys

sys.path.append('..')
import numpy as np
import msvcrt
import parameters.simulation_parameters as SIM

from Coordinate_Frames.spacecraft_viewer import SpacecraftViewer
from Coordinate_Frames.video_writer import VideoWriter
from Kinematics_and_Dynamics.data_viewer import DataViewer
from Forces_and_Moments.mav_dynamics import MavDynamics
from Forces_and_Moments.wind_simulation import WindSimulation

# Whether or not to add wind simulation?
wind_flag = True

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
mav_view = SpacecraftViewer()  # initialize the mav viewer
data_view = DataViewer()  # initialize view of data plots
if VIDEO == True:
    video = VideoWriter(video_name="chap4_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)

# initialize control commands
delta_a = 0.0  # 0.001
delta_e = 0.0  # -0.2
delta_r = 0.003  # 0.005
delta_t = 1.3  # 0.5

# initialize the simulation time
sim_time = SIM.start_time

# main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:
    # -------set control surfaces-------------
    # TODO: Use the keyboard to input the command. Need to check the document of pyqtgraph.
    # if msvcrt.kbhit():
    #     ch = msvcrt.getch()
    #     if ch == 'j':
    #         delta_a = delta_a - 0.0005  # 0.001
    #     elif ch == 'l':
    #         delta_a = delta_a + 0.0005
    #     elif ch == 'a':
    #         delta_r = delta_r - 0.005  # 0.005
    #     elif ch == 'd':
    #         delta_r = delta_r + 0.005
    #     elif ch == 's':
    #         delta_e = delta_e - 0.005
    #     elif ch == 'w':
    #         delta_e = delta_e + 0.005

    delta = np.array([[delta_e, delta_a, delta_r, delta_t]]).T  # [e, a, r, t] transpose to make it a column vector

    # -------physical system-------------
    if wind_flag is True:
        current_wind = wind.update()  # get the new wind vector
    else:
        current_wind = np.zeros((6, 1))
    mav.update_state(delta, current_wind)  # propagate the MAV dynamics

    # -------update viewer-------------
    mav_view.update(mav.msg_true_state)  # plot body of MAV
    data_view.update(mav.msg_true_state,  # true states
                     mav.msg_true_state,  # estimated states
                     mav.msg_true_state,  # commanded states
                     SIM.ts_simulation)
    if VIDEO == True:
        video.update(sim_time)

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if VIDEO == True:
    video.close()
