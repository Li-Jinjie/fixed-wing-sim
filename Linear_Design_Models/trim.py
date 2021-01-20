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
from tools.tools import Euler2Quaternion

def compute_trim(mav, Va, gamma):
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
    res = minimize(trim_objective, x0, method='SLSQP', args=(mav, Va, gamma),
                   constraints=cons, options={'ftol': 1e-10, 'disp': True})
    # extract trim state and input and return
    trim_state = np.array([res.x[0:13]]).T
    trim_input = np.array([res.x[13:17]]).T
    return trim_state, trim_input


# objective function to be minimized
# wings-level, constant-altitude flight, R* = infinite
# TODO: add R*
def trim_objective(x, mav, Va, gamma):
    # read Algorithm 14 in page 283.
    # Use quaternion instead of euler angle. Because the p,q,r = 0, according to the equations B.11-14 in page 258，
    # I set e0_dot, e1_dot, e2_dot, e3_dot = 0.

    # step 1: compute x_trim_dot using Va and gamma
    x_trim_dot = np.array([[0],  # (0) pn., don't care
                           [0],  # (1) pe., don't care
                           [-Va * np.sin(gamma)],  # (2) pd. = -h.
                           [0],  # (3) u.
                           [0],  # (4) v.
                           [0],  # (5) w.
                           [0],  # (6) e0.
                           [0],  # (7) e1.
                           [0],  # (8) e2.
                           [0],  # (9) e3.
                           [0],  # (10) p.
                           [0],  # (11) q.
                           [0]])  # (12) r.

    # step 2: compute x_dot = f(x,u) using x and {mav._forces_moments(), mav._derivatives()}
    state = np.array([x[0:13]]).T
    delta = np.array([x[13:17]]).T
    forces_moments = mav._forces_moments(delta)
    x_dot = mav._derivatives(state, forces_moments)

    # step 3: get objective function J
    # PAY ATTENTION: don't care about pn and pe

    J = np.sum(np.square((x_trim_dot - x_dot)[2:]))

    return J


if __name__ == "__main__":
    mav = mav_dynamics(0.01)

    pass
