# toy model to see the effect of switching to nonlinear model
# hopefully we'll see it's possible for alphas to decrease when switching to nonlinear model.

# results so far:
# if y = x (perfect correlation) and beta = x, alpha decreases from 2 to 1.5.
# if y = x and beta is random: similar result
# if y = beta*x (i.e. nonlinear model is correct), alpha decreases from 2.5 to 2.
# BUT once I calculate beta using the exponential formula, seems like not really a way to get alpha to drop.

import numpy as np

# generate some values for i100:
N = 1000
i100 = np.random.rand(N)

b_factor = 0.01
beta = np.divide(1-np.exp(-b_factor*i100), b_factor*i100)
y = np.multiply(i100, beta)

alpha_num = np.dot(y, i100 - np.mean(i100))
alpha_denom = np.dot(y, np.power(i100 - np.mean(i100), 2))

alpha2_num = np.dot(y, np.multiply(i100, beta) - np.mean(np.multiply(i100, beta)))
alpha2_denom = np.dot(y, np.power(np.multiply(i100, beta) - np.mean(np.multiply(i100, beta)), 2))

print(alpha_num)
print(alpha_denom)
print(alpha_num/alpha_denom)

print(alpha2_num)
print(alpha2_denom)
print(alpha2_num/alpha2_denom)
