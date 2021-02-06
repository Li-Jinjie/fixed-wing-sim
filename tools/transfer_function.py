"""
transfer function block (SISO)
num = np.array([[b1, b0]])   b1*s + b0
den = np.array([[a3, a2, a1, a0]])   a3*s^3 + a2*s^2 + a1*s^1 + a0
"""
import numpy as np
import matplotlib.pyplot as plt


class transfer_function:
    def __init__(self, num, den, Ts):
        # expects num and den to be numpy arrays of shape (1,m) and (1,n)
        m = num.shape[1]
        n = den.shape[1]
        # set initial conditions
        self._state = np.zeros((n - 1, 1))
        # make the leading coef of den == 1
        if den[0][0] != 1:
            num = num / den[0][0]
            den = den / den[0][0]  # This line will change den, so it should be executed after the change of num.
        self.num = num
        self.den = den
        # set up state space equations in control canonical form
        self._A = np.eye(n - 1)
        self._B = np.zeros((n - 1, 1))
        self._C = np.zeros((1, n - 1))
        self._B[0][0] = Ts

        for i in range(0, n - 1):
            self._A[0][i] += - Ts * den[0][i + 1]
        for i in range(1, n - 1):
            self._A[i][i - 1] += Ts

        if m == n:
            self._D = num[0][0]
            for i in range(0, m - 1):
                self._C[0][n - i - 2] = num[0][m - i - 1] - num[0][0] * den[0][n - i - 1]

        else:
            self._D = 0.0
            for i in range(0, m):
                self._C[0][n - i - 2] = num[0][m - i - 1]

    def update(self, u):
        '''Update state space model'''
        self._state = self._A @ self._state + self._B * u
        y = self._C @ self._state + self._D * u
        return y[0][0]


if __name__ == "__main__":
    # instantiate the system
    Ts = 0.01  # simulation step size
    num = np.array([[2, 3, 1]])  # numerator polynomial
    den = np.array([[5, 7, 5, 6]])  # denominator polynomial (no leading 1: s^3+4s^2+5s+6)
    system = transfer_function(num, den, Ts)
    print("A:", system._A)
    print("B:", system._B)
    print("C:", system._C)
    print("D:", system._D)

    # main simulation loop
    sim_time = 0.0
    time = [sim_time]
    output = [system.update(0.)]
    while sim_time < 10.0:
        # u = np.random.randn()  # white noise
        u = 1
        y = system.update(u)  # update based on current input
        sim_time += Ts  # increment the simulation time

        # update date for plotting
        time.append(sim_time)
        output.append(y)

    # plot output vs time
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xticks(np.arange(0, 10, 0.5))
    plt.plot(time, output)
    plt.grid()
    plt.show()
