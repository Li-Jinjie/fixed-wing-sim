"""
Class to determine wind velocity at any given moment,
calculates a steady wind speed and uses a stochastic
process to represent wind gusts. (Follows section 4.4 in uav book)
"""
import sys

sys.path.append('..')
from tools.transfer_function import TransferFunction
import numpy as np
import parameters.wind_parameters as w_para

# TODO: According to the book, Va here should be the same as Va in mav_dynamics.py. But I give it a constant here.
#  However, if Va is time-variant, the wind system becomes a time-variant system, which should be considered carefully.
Va = 20  # m/s^2
sqrt_3 = np.sqrt(3)


class WindSimulation:
    def __init__(self, Ts):
        # steady state wind defined in the inertial frame
        # change the transfer model according to the page 14 of addendum by McLain
        self._steady_state = np.array([[w_para.w_ns],
                                       [w_para.w_es],
                                       [w_para.w_ds]])  # n, e, d; m/s
        a1 = w_para.sigma_u * np.sqrt(2. * Va / (np.pi * w_para.L_u))
        b1 = Va / w_para.L_u
        self.u_w = TransferFunction(num=np.array([[a1]]),
                                    den=np.array([[1., b1]]),
                                    Ts=Ts)
        a2 = w_para.sigma_v * np.sqrt(3. * Va / (np.pi * w_para.L_v))
        a3 = a2 * Va / (sqrt_3 * w_para.L_v)
        b2 = Va / w_para.L_v
        self.v_w = TransferFunction(num=np.array([[a2, a3]]),
                                    den=np.array([[1., 2. * b2, b2 ** 2.]]),
                                    Ts=Ts)
        a4 = w_para.sigma_w * np.sqrt(3. * Va / (np.pi * w_para.L_w))
        a5 = a4 * Va / (sqrt_3 * w_para.L_w)
        b3 = Va / w_para.L_w
        self.w_w = TransferFunction(num=np.array([[a4, a5]]),
                                    den=np.array([[1., 2. * b3, b3 ** 2.]]),
                                    Ts=Ts)
        self._Ts = Ts

    def update(self):
        # returns a six vector.
        #   The first three elements are the steady state wind in the inertial frame
        #   The second three elements are the gust in the body frame
        gust = np.array([[self.u_w.update(np.random.randn())],
                         [self.v_w.update(np.random.randn())],
                         [self.w_w.update(np.random.randn())]])
        # gust = np.array([[0.],[0.],[0.]])
        return np.concatenate((self._steady_state, gust))
