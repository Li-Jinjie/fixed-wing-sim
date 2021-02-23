#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
Author: LI Jinjie
File: transfer_function_coef.py
Date: 2021/2/1 9:53
LastEditors: LI Jinjie
LastEditTime: 2021/2/1 9:53
Description: The wind parameters from page 56.
'''

import sys

sys.path.append('..')
import numpy as np

# --------steady state wind defined in the inertial frame--------
# n, e, d; m/s
w_ns = 3
w_es = -3
w_ds = 0

# ---------Dryden gust model parameters-------
mode = 0
# Choose one condition:
if mode == 0:
    # low altitude (50m), light turbulence
    L_u = 200  # m
    L_v = L_u  # m
    L_w = 50  # m
    sigma_u = 1.06  # m/s
    sigma_v = sigma_u  # m/s
    sigma_w = 0.7  # m/s
elif mode == 1:
    # low altitude (50m), moderate turbulence
    L_u = 200  # m
    L_v = L_u  # m
    L_w = 50  # m
    sigma_u = 2.12  # m/s
    sigma_v = sigma_u  # m/s
    sigma_w = 1.4  # m/s
elif mode == 2:
    # medium altitude (600m), light turbulence
    L_u = 533  # m
    L_v = L_u  # m
    L_w = 533  # m
    sigma_u = 1.5  # m/s
    sigma_v = sigma_u  # m/s
    sigma_w = 1.5  # m/s
elif mode == 3:
    # medium altitude (600m), moderate turbulence
    L_u = 533  # m
    L_v = L_u  # m
    L_w = 533  # m
    sigma_u = 3.0  # m/s
    sigma_v = sigma_u  # m/s
    sigma_w = 3.0  # m/s
