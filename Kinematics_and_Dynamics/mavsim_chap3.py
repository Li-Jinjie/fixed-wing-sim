"""
mavsimPy
    - Chapter 3 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        12/18/2018 - RWB
        1/14/2019 - RWB
"""
import sys

sys.path.append('..')
import numpy as np
import parameters.simulation_parameters as SIM

from Coordinate_Frames.spacecraft_viewer import spacecraft_viewer
from Kinematics_and_Dynamics.data_viewer import data_viewer
from Kinematics_and_Dynamics.mav_dynamics import mav_dynamics

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
mav_view = spacecraft_viewer()  # initialize the mav viewer
data_view = data_viewer()  # initialize view of data plots
if VIDEO == True:
    from Coordinate_Frames.video_writer import video_writer

    video = video_writer(video_name="chap3_video.avi",
                         bounding_box=(0, 0, 1000, 1000),
                         output_rate=SIM.ts_video)

# initialize elements of the architecture
mav = mav_dynamics(SIM.ts_simulation)

# initialize the simulation time
sim_time = SIM.start_time

# main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:
    # -------vary forces and moments to check dynamics-------------
    fx = 10
    fy = 0  # 10
    fz = 0  # 10
    Mx = 0  # 0.1
    My = 0  # 0.1
    Mz = 0  # 0.1
    forces_moments = np.array([[fx, fy, fz, Mx, My, Mz]]).T

    # -------physical system-------------
    mav.update(forces_moments)  # propagate the MAV dynamics

    # -------update viewer-------------
    mav_view.update(mav.true_state)  # plot body of MAV
    data_view.update(mav.true_state,  # true states
                     mav.true_state,  # estimated states
                     mav.true_state,  # commanded states
                     SIM.ts_simulation)
    if VIDEO == True: video.update(sim_time)

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if VIDEO == True: video.close()
