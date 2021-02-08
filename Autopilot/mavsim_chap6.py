"""
mavsim_python
    - Chapter 6 assignment for Beard & McLain, PUP, 2012
    - Last Update:
        2/5/2019 - RWB
"""
import sys

sys.path.append('..')
import numpy as np
import parameters.simulation_parameters as SIM

from Coordinate_Frames.spacecraft_viewer import SpacecraftViewer
from Kinematics_and_Dynamics.data_viewer import DataViewer
from Forces_and_Moments.mav_dynamics import MavDynamics
from Forces_and_Moments.wind_simulation import WindSimulation
from Autopilot.autopilot import AutoPilot
from tools.signals import Signals

# autopilot commands
from message_types.msg_autopilot import MsgAutopilot

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
mav_view = SpacecraftViewer()  # initialize the mav viewer
data_view = DataViewer()  # initialize view of data plots
if VIDEO is True:
    from Coordinate_Frames.video_writer import VideoWriter

    video = VideoWriter(video_name="chap6_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)
ctrl = AutoPilot(SIM.ts_simulation)

commands = MsgAutopilot()
Va_command = Signals(dc_offset=25.0, amplitude=3.0, start_time=2.0, frequency=0.01)
h_command = Signals(dc_offset=100.0, amplitude=10.0, start_time=0.0, frequency=0.02)
chi_command = Signals(dc_offset=np.radians(180), amplitude=np.radians(45), start_time=5.0, frequency=0.015)

# initialize the simulation time
sim_time = SIM.start_time

# main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:

    # -------controller-------------
    estimated_state = mav.msg_true_state  # uses true states in the control algorithm
    commands.airspeed_command = Va_command.square(sim_time)
    commands.course_command = chi_command.square(sim_time)
    commands.altitude_command = h_command.square(sim_time)
    delta, commanded_state = ctrl.update(commands, estimated_state)

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    mav.update_state(delta, current_wind)  # propagate the MAV dynamics

    # -------update viewer-------------
    mav_view.update(mav.msg_true_state)  # plot body of MAV
    data_view.update(mav.msg_true_state,  # true states
                     mav.msg_true_state,  # estimated states
                     commanded_state,  # commanded states
                     delta,
                     SIM.ts_simulation)
    if VIDEO == True: video.update(sim_time)

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if VIDEO == True: video.close()
