# quick simulation to see if there's a bias factor

# ask: is this an ok way to simulate?
# note: the bias factor changes with delta... does delta depend on wavelength?
# (because 

import numpy as np

# first assume sigma_lambda = 1 and apply eqn 14

xis = np.random.normal(size=1000000)
deltas = np.random.normal(scale=2, size=1000000)
xs = xis + deltas

numer = np.sum(xis**2)
denom = np.sum(xs**2)

print("bias:")
print(numer / denom)

# now try to incorporate tau:
constants = [0.1, 0.2, 0.5, 1, 2, 5, 10]

for c in constants:

    tau = c + np.random.normal(size=1000000) / 4
    tau_factor = (1 - np.exp(-xs)) / xs
    xis_star = xis * tau_factor
    xs_star = xis_star + deltas

    numer = np.sum(xis_star**2)
    denom = np.sum(xs_star**2)
    print("bias with tau = " + str(c) + " + noise")
    print(numer / denom)