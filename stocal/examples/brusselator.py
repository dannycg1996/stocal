"""Brusselator

A stochastic realization of the famous Brusselator system, first
proposed in I. Prigogine and R. Lefever, Symmetry Breaking Instabilities
in Dissipative Systems, J. Chem. Phys. 48, 1695 (1968).

This is a simple example of a process with only static (non-infered)
reactions. The deterministic system exhibits ascillations when b>a+1.
"""

import stocal

a = 2.
b = 10.

# Two reactions: Two molecules of A forming a dimer A2 and reverse
#r1 = MassAction({'A': 2}, {'A2': 1}, 1.0)
#r2 = MassAction({'A2': 1}, {'A': 2}, 10.0)  # Last value is the stochastic rate at which the reaction occurs

process = stocal.Process([
    stocal.MassAction({}, {"x": 1}, a),
    stocal.MassAction({"x": 2, "y": 1}, {"x": 3}, 1.),
    stocal.MassAction({"x": 1}, {"y": 1, "c": 1}, b),
    stocal.MassAction({"x": 1}, {"d": 1}, 1.),
])

if __name__ == '__main__':
    traj = process.sample({}, tmax=50)
    print("# time\tx\ty\tc\td")
    for dt, transitions in traj:
        print(traj.time, '\t'.join(str(traj.state[s]) for s in "xycd"))
