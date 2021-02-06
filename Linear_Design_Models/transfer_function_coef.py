#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
Author: LI Jinjie
File: transfer_function_coef.py
Date: 2021/2/6 9:53
LastEditors: LI Jinjie
LastEditTime: 2021/2/6 9:53
Description: Calculate the coefficients of transfer functions.
'''
from tools.tools import Euler2Quaternion, Quaternion2Euler, Euler2Rotation
import parameters.aerosonde_parameters as MAV
import numpy as np


class tf_coef:
    def __init__(self, Va, trim_state, trim_input):
        # Lateral Transfer Functions:
        # 1. Roll Angle
        rV2Sb_2 = (1 / 2) * MAV.rho * (Va ** 2) * MAV.S_wing * MAV.b
        self.a_phi_1 = -rV2Sb_2 * MAV.C_p_p * MAV.b / (2 * Va)
        self.a_phi_2 = rV2Sb_2 * MAV.C_p_delta_a

        # 2. Course and Heading

        # 3. Sideslip
        rVS_2m = (MAV.rho * Va * MAV.S_wing / (2 * MAV.mass))
        self.a_beta_1 = - rVS_2m * MAV.C_Y_beta
        self.a_beta_2 = rVS_2m * MAV.C_Y_delta_r

        # Longitudinal Transfer Functions:
        # 1. Pitch Angle
        rV2cS_2Jy = (MAV.rho * (Va ** 2) * MAV.c * MAV.S_wing / (2 * MAV.Jy))
        self.a_theta_1 = - rV2cS_2Jy * MAV.C_m_q * (MAV.c / (2 * Va))
        self.a_theta_2 = - rV2cS_2Jy * MAV.C_m_alpha
        self.a_theta_3 = rV2cS_2Jy * MAV.C_m_delta_e

        # 2. Altitude

        # 3. Airspeed
        # Pay Attention: Use trim condition!!!
        u_trim, v_trim, w_trim = trim_state.item(3), trim_state.item(4), trim_state.item(5)
        Va_trim = u_trim  # u in body frame, no wind
        delta_e_trim = trim_input[1].item()
        delta_t_trim = trim_input[3].item()
        phi_trim, theta_trim, psi_trim = Quaternion2Euler(trim_state[6:10])
        alpha_trim = np.arctan2(w_trim, u_trim)
        R_trim = Euler2Rotation(phi_trim, theta_trim, psi_trim)  # R: body to inertial
        [[_], [v_i_trim], [w_i_trim]] = R_trim @ trim_state[3:6]  # in the inertial frame
        chi_trim = np.arctan2(v_i_trim, w_i_trim)  # -pi to pi

        self.a_V_1 = (MAV.rho * Va_trim * MAV.S_wing / MAV.mass) * \
                (MAV.C_D_0 + MAV.C_D_alpha * alpha_trim + MAV.C_D_delta_e * delta_e_trim) + \
                (MAV.rho * MAV.S_prop / MAV.mass) * MAV.C_prop * Va_trim
        self.a_V_2 = (MAV.rho * MAV.S_prop / MAV.mass) * MAV.C_prop * (MAV.k_motor ** 2) * delta_t_trim
        self.a_V_3 = MAV.gravity * np.cos(theta_trim - chi_trim)
