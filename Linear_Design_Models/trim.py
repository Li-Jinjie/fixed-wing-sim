"""
compute_trim 
    - Chapter 5 assignment for Beard & McLain, PUP, 2012
    - Update history:  
        2/5/2019 - RWB
"""
import sys

sys.path.append('..')
import numpy as np
from Forces_and_Moments.mav_dynamics import mav_dynamics
from scipy.optimize import minimize
from tools.tools import Euler2Quaternion, Quaternion2Euler


def compute_trim(mav, Va, gamma, Radius):
    # define initial state and input
    state0 = mav._state
    delta0 = np.zeros([4, 1])  # delta_a, delta_e, delta_r, delta_t
    # PAY ATTENTION: initial value is VERY critical to the final solution
    delta0[1] = -0.2
    delta0[3] = 1.8  # delta_t

    x0 = np.concatenate((state0, delta0), axis=0)
    # define equality constraints
    # An explanation about lambda: https://www.w3schools.com/python/python_lambda.asp
    # An explanation about scipy.optimize.minimize: https://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html
    cons = ({'type': 'eq',
             'fun': lambda x: np.array([
                 x[3] ** 2 + x[4] ** 2 + x[5] ** 2 - Va ** 2,  # magnitude of velocity vector is Va
                 x[4],  # v=0, force side velocity to be zero
                 x[6] ** 2 + x[7] ** 2 + x[8] ** 2 + x[9] ** 2 - 1.,  # force quaternion to be unit length
                 x[7],  # e1=0  - forcing e1=e3=0 ensures zero roll and zero yaw in trim
                 x[9],  # e3=0
                 x[10],  # p=0  - angular rates should all be zero
                 x[11],  # q=0
                 x[12],  # r=0
             ]),
             'jac': lambda x: np.array([
                 [0., 0., 0., 2 * x[3], 2 * x[4], 2 * x[5], 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0., 2 * x[6], 2 * x[7], 2 * x[8], 2 * x[9], 0., 0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.],
             ])
             })
    # 'jac' means the Jacobin matrix of constrain functions to x
    # solve the minimization problem to find the trim states and inputs
    res = minimize(trim_objective, x0, method='SLSQP', args=(mav, Va, gamma, Radius),
                   constraints=cons, options={'ftol': 1e-10, 'disp': True})
    # extract trim state and input and return
    trim_state = np.array([res.x[0:13]]).T
    trim_input = np.array([res.x[13:17]]).T
    return trim_state, trim_input


# objective function to be minimized
def trim_objective(x, mav, Va, gamma, Radius):
    # Read Algorithm 14 in page 283.
    # Watch the solution about chapter 5 on Youtube.
    # Use quaternion instead of euler angle.

    # Step 1: compute x_trim_dot using Va, gamma and Radius
    #       (1) Compute p, q, r according to (3.3) in page 31.
    phi_dot = 0
    theta_dot = 0
    psi_dot = Va / Radius * np.cos(gamma)
    phi, theta, psi = Quaternion2Euler(x[6:10])

    rot_tran = np.array([[1, 0, -np.sin(theta)],
                         [0, np.cos(phi), np.sin(phi) * np.cos(theta)],
                         [0, -np.sin(phi), np.cos(phi) * np.cos(theta)]])
    pqr = rot_tran @ np.array([[phi_dot, theta_dot, psi_dot]]).T
    p = pqr.item(0)
    q = pqr.item(1)
    r = pqr.item(2)

    #       (2) Compute e0_dot, e1_dot, e2_dot, e3_dot according to B.3 in page 256
    e0, e1, e2, e3 = x.item(6), x.item(7), x.item(8), x.item(9)
    e0_dot = -1 / 2 * (p * e1 + q * e2 + r * e3)
    e1_dot = 1 / 2 * (p * e0 + r * e2 - q * e3)
    e2_dot = 1 / 2 * (q * e0 - r * e1 + p * e3)
    e3_dot = 1 / 2 * (r * e0 + q * e1 - p * e2)

    #       (3) Compute x_trim_dot
    x_trim_dot = np.array([[0],  # (0) pn_dot, don't care
                           [0],  # (1) pe_dot, don't care
                           [-Va * np.sin(gamma)],  # (2) pd_dot = -h.
                           [0],  # (3) u_dot
                           [0],  # (4) v_dot
                           [0],  # (5) w_dot
                           [e0_dot],  # (6) e0_dot
                           [e1_dot],  # (7) e1_dot
                           [e2_dot],  # (8) e2_dot
                           [e3_dot],  # (9) e3_dot
                           [0],  # (10) p_dot
                           [0],  # (11) q_dot
                           [0]])  # (12) r_dot

    # Step 2: compute x_dot = f(x,u) using x and {mav._forces_moments(), mav._derivatives()}
    state = np.array([x[0:13]]).T
    delta = np.array([x[13:17]]).T
    mav._state = state
    mav._update_velocity_data()  # update Va, alpha, beta to match new states
    forces_moments = mav._forces_moments(delta)
    x_dot = mav._derivatives(state, forces_moments)

    # Step 3: get objective function J
    # PAY ATTENTION: don't care about pn and pe

    J = np.sum(np.square((x_trim_dot - x_dot)[2:]))

    return J


if __name__ == "__main__":
    pass
