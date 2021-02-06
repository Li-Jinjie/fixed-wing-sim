"""
mavSimPy 
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        1/1/2019 - RWB
        1/29/2019 - RWB
        2/2/2019 - RWB
"""
import sys

sys.path.append('..')
import numpy as np
import parameters.simulation_parameters as SIM

from Coordinate_Frames.spacecraft_viewer import SpacecraftViewer
from Coordinate_Frames.video_writer import VideoWriter
from Kinematics_and_Dynamics.data_viewer import DataViewer
from Forces_and_Moments.mav_dynamics import MavDynamics
from Forces_and_Moments.wind_simulation import WindSimulation
from Linear_Design_Models.trim import compute_trim
from Linear_Design_Models.compute_models import compute_ss_model, compute_tf_model

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
mav_view = SpacecraftViewer()  # initialize the mav viewer
data_view = DataViewer()  # initialize view of data plots
if VIDEO == True:
    video = VideoWriter(video_name="chap5_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)

# use compute_trim function to compute trim state and trim input
Va = 30.
gamma = 5. * np.pi / 180.
Radius = -150  # straight line: np.inf

trim_state, trim_input = compute_trim(mav, Va, gamma, Radius)
# print("trim states: \r\n", trim_state)
# print("trim input: a e r t \r\n", trim_input)

mav._state = trim_state  # set the initial state of the mav to the trim state
mav._update_velocity_data()  # update Va, alpha, beta to match the new state above.
delta = trim_input  # set input to constant constant trim input

# # compute the state space model linearized about trim
A_lon, B_lon, A_lat, B_lat = compute_ss_model(mav, trim_state, trim_input)

T_phi_delta_a, T_chi_phi, T_beta_delta_r, T_theta_delta_e, T_h_theta, T_h_Va, T_Va_delta_t, T_Va_theta \
    = compute_tf_model(mav, trim_state, trim_input)

# initialize the simulation time
sim_time = SIM.start_time

# main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:

    # -------physical system-------------
    # current_wind = wind.update()  # get the new wind vector
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
