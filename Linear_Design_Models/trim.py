"""
compute_trim 
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        12/29/2018 - RWB
"""
import sys

sys.path.append('..')
import numpy as np
from scipy.optimize import minimize
from tools.rotations import euler_2_quaternion, quaternion_2_euler
from message_types.msg_delta import MsgDelta


def compute_trim(mav, Va, gamma, Radius):
    # define initial state and input
    state0 = mav._state
    # PAY ATTENTION: The initial value is VERY CRITICAL to the final solution. [0,0,0,1] in p278.
    delta0 = MsgDelta(throttle=1)

    x0 = np.concatenate((state0, delta0.to_array()), axis=0)
    # define equality constraints
    # An explanation about lambda: https://www.w3schools.com/python/python_lambda.asp
    # An explanation about scipy.optimize.minimize: https://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html
    cons = ({'type': 'eq',
             'fun': lambda x: np.array([
                 x[3] ** 2 + x[4] ** 2 + x[5] ** 2 - Va ** 2,  # magnitude of velocity vector is Va
                 x[4],  # v=0, force side velocity to be zero
                 x[6] ** 2 + x[7] ** 2 + x[8] ** 2 + x[9] ** 2 - 1.,  # force quaternion to be unit length
                 # x[7],  # e1=0  - forcing e1=e3=0 ensures zero roll and zero yaw in trim
                 # x[9],  # e3=0
                 # x[10],  # p=0  - angular rates should all be zero
                 # x[11],  # q=0
                 # x[12],  # r=0
             ]),
             'jac': lambda x: np.array([
                 [0., 0., 0., 2 * x[3], 2 * x[4], 2 * x[5], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0., 2 * x[6], 2 * x[7], 2 * x[8], 2 * x[9], 0., 0., 0., 0., 0., 0., 0.],
                 # [0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                 # [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.],
                 # [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],
                 # [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
                 # [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.],
             ])
             })
    # 'jac' means the Jacobin matrix of constrain functions to x
    # solve the minimization problem to find the trim states and inputs
    res = minimize(trim_objective_fun, x0, method='SLSQP', args=(mav, Va, gamma, Radius),
                   constraints=cons, options={'ftol': 1e-10, 'disp': True})
    # extract trim state and input and return
    trim_state = np.array([res.x[0:13]]).T
    trim_input = MsgDelta(elevator=res.x.item(13), aileron=res.x.item(14),
                          rudder=res.x.item(15), throttle=res.x.item(16))
    # trim_input.print()
    # print('trim_state=', trim_state.T)
    return trim_state, trim_input


# objective function to be minimized
def trim_objective_fun(x, mav, Va, gamma, Radius):
    # Read Algorithm 14 in page 283.
    # Watch the solution about chapter 5 on Youtube.
    # Use euler angle only. Do not use quaternion here since quaternion_trim_dot cannot be calculated without e0_trim,
    # e1_trim, e2_trim and e3_trim, which are impossible to get based on Va, gamma and Radius.

    # Step 1: compute x_euler_trim_dot using, state_*, Va, gamma and Radius

    x_euler_trim_dot = np.array([[0],  # (0) pn_dot, don't care
                                 [0],  # (1) pe_dot, don't care
                                 [-Va * np.sin(gamma)],  # (2) pd_dot = -h.
                                 [0],  # (3) u_dot
                                 [0],  # (4) v_dot
                                 [0],  # (5) w_dot
                                 [0],  # (6) phi_dot
                                 [0],  # (7) theta_dot
                                 [Va / Radius * np.cos(gamma)],  # (8) psi_dot
                                 [0],  # (9) p_dot
                                 [0],  # (10) q_dot
                                 [0]])  # (11) r_dot

    # Step 2: compute x_quat_dot = f(x,u) using x and {mav._forces_moments(), mav._derivatives()}
    state = np.array([x[0:13]]).T
    delta = MsgDelta(elevator=x.item(13), aileron=x.item(14), rudder=x.item(15), throttle=x.item(16))

    mav._state = state
    mav._update_velocity_data()  # update Va, alpha, beta to match new states
    forces_moments = mav._forces_moments(delta)
    x_quat_dot = mav._derivatives(state, forces_moments)

    # Step 3: # convert x_quat_dot to x_euler_dot using x, refer to (3.3) in page 31
    phi, theta, psi = quaternion_2_euler(x[6:10])
    p, q, r = x.item(10), x.item(11), x.item(12)

    rot_tran = np.array([[1, np.sin(phi) * np.tan(theta), np.cos(phi) * np.tan(theta)],
                         [0, np.cos(phi), -np.sin(phi)],
                         [0, np.sin(phi) / np.cos(theta), np.cos(phi) / np.cos(theta)]])
    ptp_dot = rot_tran @ np.array([[p, q, r]]).T
    phi_dot = ptp_dot.item(0)
    theta_dot = ptp_dot.item(1)
    psi_dot = ptp_dot.item(2)
    x_euler_dot = np.zeros([12, 1])
    x_euler_dot[0:6] = x_quat_dot[0:6]
    x_euler_dot[6:9] = np.array([[phi_dot, theta_dot, psi_dot]]).T
    x_euler_dot[9:] = x_quat_dot[10:]

    # Step 4: get objective function J
    # PAY ATTENTION: don't care about pn and pe

    J = np.sum(np.square((x_euler_trim_dot - x_euler_dot)[2:]))

    return J


if __name__ == "__main__":
    pass
