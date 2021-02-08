'''
some coefficients refer to the code on uavbook-gitlab
'''

import sys

sys.path.append('..')
import numpy as np
import parameters.aerosonde_parameters as MAV
from Linear_Design_Models.transfer_function_coef import TfCoef
from Linear_Design_Models.trim import compute_trim
from Forces_and_Moments.mav_dynamics import MavDynamics
import parameters.simulation_parameters as SIM

mav = MavDynamics(SIM.ts_simulation)

gravity = MAV.gravity
sigma = 0.05
Va0 = 35  # m/s according to the Comments on https://uavbook.byu.edu/doku.php

trim_state, trim_input = compute_trim(mav, Va0, 0, np.inf)
tf = TfCoef(Va0, trim_state, trim_input)

delta_a_max = np.radians(45)  # aileron
phi_max = np.radians(30)  # roll
delta_r_max = np.radians(45)  # rudder
delta_e_max = limit = np.radians(45)  # elevator
theta_max = np.radians(30)  # pitch
delta_t_max = 1.0

# ----------roll loop-------------
wn_roll = 10
zeta_roll = 0.707
roll_kp = (wn_roll ** 2) / tf.a_phi_2
roll_kd = (2. * zeta_roll * wn_roll - tf.a_phi_1) / tf.a_phi_2

# ----------course loop-------------
wn_course = wn_roll / 10.0
zeta_course = 1.0
course_kp = 2. * zeta_course * wn_course * Va0 / gravity
course_ki = (wn_course ** 2) * Va0 / gravity

# ----------sideslip loop-------------
wn_sideslip = 20
zeta_sideslip = 0.707
sideslip_kp = (2. * zeta_sideslip * wn_sideslip - tf.a_beta_1) / tf.a_beta_2
sideslip_ki = (wn_sideslip ** 2) / tf.a_beta_2

# ----------yaw damper-------------
yaw_damper_p_wo = 0.45
yaw_damper_kr = 0.2

# ----------pitch loop-------------
wn_pitch = 24.0
zeta_pitch = 0.707
pitch_kp = (wn_pitch ** 2 - tf.a_theta_2) / tf.a_theta_3
pitch_kd = (2. * zeta_pitch * wn_pitch - tf.a_theta_1) / tf.a_theta_3
K_theta_DC = pitch_kp * tf.a_theta_3 / (tf.a_theta_2 + pitch_kp * tf.a_theta_3)

# ----------altitude loop-------------
wn_altitude = wn_pitch / 15  # 5-15 refer to the book
zeta_altitude = 1.0
altitude_kp = 2. * zeta_altitude * wn_altitude / (K_theta_DC * Va0)
altitude_ki = wn_altitude ** 2 / (K_theta_DC * Va0)
altitude_zone = 10.0  # moving saturation limit around current altitude

# ---------airspeed hold using throttle---------------
wn_airspeed_throttle = 3.0
zeta_airspeed_throttle = 2  # 0.707

airspeed_throttle_kp = (2 * zeta_airspeed_throttle * wn_airspeed_throttle - tf.a_V_1) / tf.a_V_2
airspeed_throttle_ki = wn_airspeed_throttle ** 2 / tf.a_V_2

del tf
