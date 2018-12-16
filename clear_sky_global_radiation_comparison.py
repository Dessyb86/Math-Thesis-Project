import numpy as np
import math
import pylab
import matplotlib.pyplot as plt

Q_on_max = 1374

Q_on_min = 1360

B_10 = 0.1310428

B_11 = 0.7342016

k_t = 0.1368551

x = np.linspace(0,1.57,500)

Phi_1 = 990 * np.cos(x) - 30

Phi_2_max = (B_10 + B_11 * np.exp(-k_t / np.cos(x))) * Q_on_max + Q_on_max * np.cos(x) * (0.271 - 0.294 * (B_10 + B_11 * np.exp(-k_t / np.cos(x))))

Phi_2_min = (B_10 + B_11 * np.exp(-k_t / np.cos(x))) * Q_on_min + Q_on_min * np.cos(x) * (0.271 - 0.294 * (B_10 + B_11 * np.exp(-k_t / np.cos(x))))

# compose plot
pylab.plot(x,Phi_1,color='red',linestyle = 'solid')

pylab.plot(x,Phi_2_max,color='green',linestyle = 'solid')

pylab.plot(x,Phi_2_min,color='blue',linestyle = 'solid')

plt.xlabel('Solar Zanith Angle (radians)', fontsize=12)
plt.ylabel('Clear Sky Global Radiation (W.m**(-2))', fontsize=12)
plt.savefig('SC_plot')
plt.clf()
plt.cla()
plt.close()
