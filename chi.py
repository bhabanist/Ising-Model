import numpy as np
h = 1.0
# load data from dat file
data = np.loadtxt('bethe_mag_1.0.txt')
beta = data[:, 0]
magnetization = data[:, 1]
hp = data[:, 2]

# calculate lambd
lam = 3*np.sinh(2*beta)/(np.cosh(2*(hp+h)) + np.cosh(2*beta))

# calculate dh'/dh
dhdh = lam/(1-lam)

# calculate chi
chi = 2*(1 + dhdh)*(1+np.exp(-2*beta)*np.cosh(2*(hp+h)))/np.power(np.cosh(2*(hp+h)) + np.exp(-2*beta), 2)

# save results to a new dat file
results = np.column_stack((beta, magnetization, hp, chi))
np.savetxt('bethe_1.0.txt', results)

