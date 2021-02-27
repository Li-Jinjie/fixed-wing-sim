# dubins_parameters
#   - Dubins parameters that define path between two configurations
#
# mavsim_matlab 
#     - Beard & McLain, PUP, 2012
#     - Update history:  
#         3/26/2019 - RWB
#         4/2/2020 - RWB

import numpy as np
import sys

sys.path.append('..')


class DubinsParameters:
    def __init__(self, ps=9999 * np.ones((3, 1)), chis=9999,
                 pe=9999 * np.ones((3, 1)), chie=9999, R=9999):
        if R == 9999:  # fly in a circle
            L = R
            cs = ps
            lams = R
            ce = ps
            lame = R
            w1 = ps
            q1 = ps
            w2 = ps
            w3 = ps
            q3 = ps
        else:
            L, cs, lams, ce, lame, w1, q1, w2, w3, q3 \
                = compute_parameters(ps, chis, pe, chie, R)
        self.p_s = ps
        self.chi_s = chis
        self.p_e = pe
        self.chi_e = chie
        self.radius = R
        self.length = L
        self.center_s = cs
        self.dir_s = lams
        self.center_e = ce
        self.dir_e = lame
        self.r1 = w1
        self.n1 = q1
        self.r2 = w2
        self.r3 = w3
        self.n3 = q3

    def update(self, ps, chis, pe, chie, R):
        L, cs, lams, ce, lame, w1, q1, w2, w3, q3 \
            = compute_parameters(ps, chis, pe, chie, R)
        self.p_s = ps
        self.chi_s = chis
        self.p_e = pe
        self.chi_e = chie
        self.radius = R
        self.length = L
        self.center_s = cs
        self.dir_s = lams
        self.center_e = ce
        self.dir_e = lame
        self.r1 = w1
        self.n1 = q1
        self.r2 = w2
        self.r3 = w3
        self.n3 = q3


def compute_parameters(ps, chis, pe, chie, R):
    ell = np.linalg.norm(ps - pe)
    pi = np.pi
    if ell < 2 * R:  # TODO: figure out why, not 3R on the book.
        print('Error in Dubins Parameters: The distance between nodes must be larger than 2R.')
    else:
        # compute start and end circles
        crs = ps + R * rotz(pi / 2.) @ np.array([[np.cos(chis), np.sin(chis), 0]]).T
        cls = ps + R * rotz(-pi / 2.) @ np.array([[np.cos(chis), np.sin(chis), 0]]).T
        cre = pe + R * rotz(pi / 2.) @ np.array([[np.cos(chie), np.sin(chie), 0]]).T
        cle = pe + R * rotz(-pi / 2.) @ np.array([[np.cos(chie), np.sin(chie), 0]]).T

        # compute L1
        line = cre - crs
        theta_var = np.arctan2(line.item(1), line.item(0))

        angle_1 = mod(2. * pi + mod(theta_var - pi / 2.) - mod(chis - pi / 2.))
        angle_2 = mod(2. * pi + mod(chie - pi / 2.) - mod(theta_var - pi / 2.))

        L1 = np.linalg.norm(crs - cre) + R * angle_1 + R * angle_2

        # compute L2
        line = cle - crs
        ell = np.linalg.norm(line)
        theta_var = np.arctan2(line.item(1), line.item(0))
        theta_var_2 = theta_var - pi / 2. + np.arcsin(2. * R / ell)
        angle_1 = mod(2. * pi + mod(theta_var_2) - mod(chis - pi / 2.))
        angle_2 = mod(2. * pi + mod(theta_var_2 + pi) - mod(chie + pi / 2.))
        if not np.isreal(angle_2):  # TODO: figure out why. I guess it should be np.inf
            L2 = np.inf
            print("angle_2 is not a real numble!")
        else:
            L2 = np.sqrt(ell ** 2 - 4. * (R ** 2)) + R * angle_1 + R * angle_2

        # compute L3
        line = cre - cls
        ell = np.linalg.norm(line)
        theta_var = np.arctan2(line.item(1), line.item(0))
        theta_var_2 = np.arccos(2. * R / ell)
        angle_1 = mod(2. * pi + mod(chis + pi / 2.) - mod(theta_var + theta_var_2))
        angle_2 = mod(2. * pi + mod(chie - pi / 2.) - mod(theta_var + theta_var_2 - pi))
        if not np.isreal(angle_2):  # TODO: figure out why
            L3 = np.inf
        else:
            L3 = np.sqrt(ell ** 2 - 4. * (R ** 2)) + R * angle_1 + R * angle_2

        # compute L4
        line = cle - cls
        ell = np.linalg.norm(line)
        theta_var = np.arctan2(line.item(1), line.item(0))
        angle_1 = mod(2. * pi + mod(chis + pi / 2.) - mod(theta_var + pi / 2.))
        angle_2 = mod(2. * pi + mod(theta_var + pi / 2.) - mod(chie + pi / 2.))
        L4 = np.linalg.norm(cls - cle) + R * angle_1 + R * angle_2
        # L is the minimum distance
        L = np.min([L1, L2, L3, L4])
        idx = np.argmin([L1, L2, L3, L4])
        e1 = np.array([[1., 0., 0.]]).T
        if idx == 0:
            cs = crs
            lams = +1
            ce = cre
            lame = +1
            q1 = (ce - cs) / np.linalg.norm(ce - cs)
            w1 = cs + R * rotz(-pi / 2.) @ q1
            w2 = ce + R * rotz(-pi / 2.) @ q1
        elif idx == 1:
            cs = crs
            lams = +1
            ce = cle
            lame = -1

            line = ce - cs
            ell = np.linalg.norm(line)
            theta_var = np.arctan2(line.item(1), line.item(0))
            theta_var_2 = theta_var - pi / 2. + np.arcsin(2. * R / ell)

            q1 = rotz(theta_var_2 + pi / 2.) @ e1
            w1 = cs + R * rotz(theta_var_2) @ e1
            w2 = ce + R * rotz(theta_var_2 + pi) @ e1
        elif idx == 2:
            cs = cls
            lams = -1
            ce = cre
            lame = +1
            q1 = rotz(theta_var + theta_var_2 - pi / 2.) @ e1
            w1 = cs + R * rotz(theta_var + theta_var_2) @ e1
            w2 = ce + R * rotz(theta_var + theta_var_2 - pi) @ e1
        elif idx == 3:
            cs = cls
            lams = -1
            ce = cle
            lame = -1
            q1 = (ce - cs) / np.linalg.norm(ce - cs)
            w1 = cs + R * rotz(pi / 2.) @ q1
            w2 = ce + R * rotz(pi / 2.) @ q1
        w3 = pe
        q3 = rotz(chie) @ e1

        return L, cs, lams, ce, lame, w1, q1, w2, w3, q3


def rotz(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta), np.cos(theta), 0],
                     [0, 0, 1]])


def mod(x):
    while x < 0:
        x += 2 * np.pi
    while x > 2 * np.pi:
        x -= 2 * np.pi
    return x
