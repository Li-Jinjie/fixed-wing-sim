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
import parameters.simulation_parameters as SIM

from Coordinate_Frames.mav_viewer import MavViewer
from Kinematics_and_Dynamics.data_viewer import DataViewer
from Forces_and_Moments.mav_dynamics import MavDynamics
from Forces_and_Moments.wind_simulation import WindSimulation
from message_types.msg_delta import MsgDelta

# Whether or not to add wind simulation?
wind_flag = True

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
mav_view = MavViewer()  # initialize the mav viewer
data_view = DataViewer()  # initialize view of data plots
if VIDEO is True:
    from Coordinate_Frames.video_writer import VideoWriter

    video = VideoWriter(video_name="chap4_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)
delta = MsgDelta()

# initialize the simulation time
sim_time = SIM.start_time
plot_time = sim_time

# main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:
    # -------set control surfaces-------------
    delta.elevator = -0.05
    delta.aileron = 0.002
    delta.rudder = -0.003
    delta.throttle = 0.7
    # transpose to make it a column vector

    # -------physical system-------------
    if wind_flag is True:
        current_wind = wind.update()  # get the new wind vector
    else:
        current_wind = np.zeros((6, 1))
    mav.update(delta, current_wind)  # propagate the MAV dynamics

    # -------update viewer-------------
    if sim_time - plot_time > SIM.ts_plotting:
        mav_view.update(mav.true_state)  # plot body of MAV
        plot_time = sim_time
    data_view.update(mav.true_state,  # true states
                     mav.true_state,  # estimated states
                     mav.true_state,  # commanded states
                     delta,  # inputs to aircraft
                     SIM.ts_simulation)
    if VIDEO is True:
        video.update(sim_time)

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if VIDEO is True:
    video.close()
