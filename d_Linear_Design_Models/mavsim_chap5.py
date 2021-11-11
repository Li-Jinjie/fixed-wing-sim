"""
mavsim_python
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Last Update:
        2/2/2019 - RWB
"""
import sys

sys.path.append('..')
import numpy as np
import parameters.simulation_parameters as SIM

from a_Coordinate_Frames.mav_viewer import MavViewer
from b_Kinematics_and_Dynamics.data_viewer import DataViewer
from c_Forces_and_Moments.mav_dynamics import MavDynamics
from c_Forces_and_Moments.wind_simulation import WindSimulation
from d_Linear_Design_Models.trim import compute_trim
from d_Linear_Design_Models.compute_models import compute_model
from tools.signals import Signals

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
mav_view = MavViewer()  # initialize the mav viewer
data_view = DataViewer()  # initialize view of data plots
if VIDEO is True:
    from a_Coordinate_Frames.video_writer import VideoWriter

    video = VideoWriter(video_name="chap5_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)

# use compute_trim function to compute trim state and trim input
Va = 25.
gamma = 0. * np.pi / 180.
Radius = np.inf  # straight line: np.inf

trim_state, trim_input = compute_trim(mav, Va, gamma, Radius)
# print("trim states: \r\n", trim_state)
# print("trim input: a e r t \r\n", trim_input)

mav.external_set_state(trim_state)  # set the initial state of the mav to the trim state
mav._update_velocity_data()  # update Va, alpha, beta to match the new state above.
delta = trim_input  # set input to constant constant trim input

# # compute the state space model linearized about trim
compute_model(mav, trim_state, trim_input)

# this signal will be used to excite modes
input_signal = Signals(amplitude=.05,
                       duration=0.01,
                       start_time=2.0)

# initialize the simulation time
sim_time = SIM.start_time

# main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:

    # -------physical system-------------
    # current_wind = wind.update()  # get the new wind vector
    current_wind = np.zeros((6, 1))
    # this input excites the phugoid mode by adding an impulse at t=5.0
    # delta.elevator += input_signal.impulse(sim_time)
    # delta.rudder += input_signal.doublet(sim_time)
    mav.update(delta, current_wind)  # propagate the MAV dynamics

    # -------update viewer-------------
    mav_view.update(mav.true_state)  # plot body of MAV
    data_view.update(mav.true_state,  # true states
                     mav.true_state,  # estimated states
                     mav.true_state,  # commanded states
                     delta,  # input to aircraft
                     SIM.ts_simulation)
    if VIDEO is True:
        video.update(sim_time)

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if VIDEO is True:
    video.close()
