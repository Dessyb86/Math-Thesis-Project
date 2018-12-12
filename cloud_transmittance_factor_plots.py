import numpy as np
import pylab
import matplotlib.pyplot as plt

x = np.linspace(0,1,500)
y = 1-0.75*x**3.4

plt.plot(x,y,color='black',linestyle='dashed')

xlist=np.linspace(0.0001,1,100)
ylist=np.linspace(0,1,100)
X,Y = np.meshgrid(xlist,ylist)
Z = 1 - ((1 - Y)/X**2)
CS = plt.contour(X,Y,Z,[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],linestyles='solid')
plt.clabel(CS, inline=1, fontsize=10)
plt.title('The Effect of Total Cloud Cover on Global Solar Radiation')
plt.xlabel('Fractional Total Cloud Cover (dimensionless)', fontsize=12)
plt.ylabel('Cloud Transmittance Factor (dimensionless)', fontsize=12)
plt.savefig('a_plot.png')
pylab.show() # show the plot
plt.clf()
plt.cla()
plt.close()
