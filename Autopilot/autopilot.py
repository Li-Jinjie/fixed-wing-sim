"""
autopilot block for mavsim_python
    - Beard & McLain, PUP, 2012
    - Last Update:
        2/6/2019 - RWB
"""
import sys
import numpy as np

sys.path.append('..')
import parameters.control_parameters as AP
from tools.transfer_function import TransferFunction
from tools.wrap import wrap
from Autopilot.pi_control import PIControl
from Autopilot.pd_control_with_rate import PDControlWithRate
from message_types.msg_state import MsgState
from message_types.msg_delta import MsgDelta


class Autopilot:
    def __init__(self, ts_control):
        # instantiate lateral controllers
        self.roll_from_aileron = PDControlWithRate(
            kp=AP.roll_kp,
            kd=AP.roll_kd,
            limit=np.radians(45))
        self.course_from_roll = PIControl(
            kp=AP.course_kp,
            ki=AP.course_ki,
            Ts=ts_control,
            limit=np.radians(30))
        self.yaw_damper = TransferFunction(
            num=np.array([[AP.yaw_damper_kr, 0]]),
            den=np.array([[1, AP.yaw_damper_p_wo]]),
            Ts=ts_control)

        # instantiate lateral controllers
        self.pitch_from_elevator = PDControlWithRate(
            kp=AP.pitch_kp,
            kd=AP.pitch_kd,
            limit=np.radians(45))
        self.altitude_from_pitch = PIControl(
            kp=AP.altitude_kp,
            ki=AP.altitude_ki,
            Ts=ts_control,
            limit=np.radians(30))
        self.airspeed_from_throttle = PIControl(
            kp=AP.airspeed_throttle_kp,
            ki=AP.airspeed_throttle_ki,
            Ts=ts_control,
            limit=1.0)
        self.commanded_state = MsgState()

    def update(self, cmd, state):

        # lateral autopilot
        chi_c = wrap(cmd.course_command, state.chi)  # pay attention to wrap() here
        phi_c = self.saturate(cmd.phi_feedforward + self.course_from_roll.update(chi_c, state.chi), -np.radians(30),
                              np.radians(30))
        delta_a = self.roll_from_aileron.update(phi_c, state.phi, state.p)  # phi_dot = p
        delta_r = self.yaw_damper.update(state.r)  # uncomment this line if the side-slip angle is unmeasurable.

        # longitudinal autopilot
        # saturate the altitude command
        altitude_c = self.saturate(cmd.altitude_command, state.altitude - AP.altitude_zone,
                                   state.altitude + AP.altitude_zone)
        theta_c = self.altitude_from_pitch.update(altitude_c, state.altitude)
        delta_e = self.pitch_from_elevator.update(theta_c, state.theta, state.q)  # theta_dot is approximately ~ q
        delta_t = self.airspeed_from_throttle.update(cmd.airspeed_command, state.Va)  # delta_t_* is not known
        delta_t = self.saturate(delta_t, 0.0, 1.0)

        # construct output and commanded states
        delta = MsgDelta(elevator=delta_e,
                         aileron=delta_a,
                         rudder=delta_r,
                         throttle=delta_t)
        self.commanded_state.altitude = cmd.altitude_command
        self.commanded_state.Va = cmd.airspeed_command
        self.commanded_state.phi = phi_c
        self.commanded_state.theta = theta_c
        self.commanded_state.chi = cmd.course_command
        return delta, self.commanded_state

    def saturate(self, input, low_limit, up_limit):
        if input <= low_limit:
            output = low_limit
        elif input >= up_limit:
            output = up_limit
        else:
            output = input
        return output
