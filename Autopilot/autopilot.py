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
from Autopilot.pid_control import PidControl, PiControl, PdControlWithRate
from message_types.msg_state import MsgState


class AutoPilot:
    def __init__(self, ts_control):
        # instantiate lateral controllers
        self.roll_from_aileron = PdControlWithRate(
            kp=AP.roll_kp,
            kd=AP.roll_kd,
            limit=np.radians(45))
        self.course_from_roll = PiControl(
            kp=AP.course_kp,
            ki=AP.course_ki,
            Ts=ts_control,
            limit=np.radians(30))
        self.sideslip_from_rudder = PiControl(
            kp=AP.sideslip_kp,
            ki=AP.sideslip_ki,
            Ts=ts_control,
            limit=np.radians(45))

        # Please read the page 38 in uavbook_supplement.pdf
        self.yaw_damper = TransferFunction(
            num=np.array([[AP.yaw_damper_kr, 0]]),
            den=np.array([[1, 1 / AP.yaw_damper_p_wo]]),
            Ts=ts_control)

        # instantiate lateral controllers
        self.pitch_from_elevator = PdControlWithRate(
            kp=AP.pitch_kp,
            kd=AP.pitch_kd,
            limit=np.radians(45))
        self.altitude_from_pitch = PiControl(
            kp=AP.altitude_kp,
            ki=AP.altitude_ki,
            Ts=ts_control,
            limit=np.radians(30))
        self.airspeed_from_throttle = PiControl(
            kp=AP.airspeed_throttle_kp,
            ki=AP.airspeed_throttle_ki,
            Ts=ts_control,
            limit=1.0)
        self.commanded_state = MsgState()

    def update(self, cmd, state):
        # state is MsgState(), store a lot of information

        # lateral autopilot
        chi_c = wrap(cmd.course_command, state.chi)  # attention to wrap() here
        phi_c = self.saturate(cmd.phi_feedforward + self.course_from_roll.update(chi_c, state.chi), -np.radians(30),
                              np.radians(30))
        delta_a = self.roll_from_aileron.update(phi_c, state.phi, state.p)  # phi_dot = p
        # delta_r = self.sideslip_from_rudder.update(0, state.beta)
        delta_r = self.yaw_damper.update(state.r)  # uncomment this line if the side-slip angle is unmeasurable.

        # longitudinal autopilot
        h_c = self.saturate(cmd.altitude_command, state.h - AP.altitude_zone, state.h + AP.altitude_zone)
        theta_c = self.altitude_from_pitch.update(h_c, state.h)
        delta_e = self.pitch_from_elevator.update(theta_c, state.theta, state.q)  # theta_dot is approximately ~ q
        delta_t = self.airspeed_from_throttle.update(cmd.airspeed_command, state.Va)  # delta_t_* is not known
        delta_t = self.saturate(delta_t, 0.0, 1.0)

        # construct output and commanded states
        delta = np.array([[delta_e], [delta_a], [delta_r], [delta_t]])
        self.commanded_state.h = cmd.altitude_command
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
