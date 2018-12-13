import numpy as np
import pylab
import matplotlib.pyplot as plt

#creating line plot using phi_sc_1 function

x = np.linspace(0,1,500) # creates x - an array of 500 evenly spaced values of fractional total cloud cover on in interval [0,1]

y = 1-0.75*x**3.4 #creates y - an array of 500 values of the cloud transmittance factor on the interval [0,1] that corresponds to the array x

plt.plot(x,y,color='black',linestyle='dashed') #creates line plot of y against x

xlist=np.linspace(0.0001,1,100) #creates an array of 100 evenly spaced values of fractional total cloud cover on in interval [0,1]

ylist=np.linspace(0,1,100) #creates an array of 100 evenly spaced values of the cloud transmittance factor on in interval [0,1]

X,Y = np.meshgrid(xlist,ylist) #creates a multi-dimensional mesh grid

Z = 1 - ((1 - Y)/X**2) #creates Z - an array of the cloud attenuation factor for each corresponding values of X and Y

CS = plt.contour(X,Y,Z,[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],linestyles='solid') #creates level curve plots for different values of Z on a plot of Y against X

plt.clabel(CS, inline=1, fontsize=10) #labels each level curve

plt.title('The Effect of Total Cloud Cover on Global Solar Radiation') #creates title for the plot

plt.xlabel('Fractional Total Cloud Cover (dimensionless)', fontsize=12) #labels the horizontal axis

plt.ylabel('Cloud Transmittance Factor (dimensionless)', fontsize=12) #labels the y-axis

plt.savefig('a_plot.png') #creates "a_plot.png" - a png image file of the aforementioned plots on the same graph

plt.clf() #clears the current figure

plt.cla() #clears the current axis

plt.close()#closes the image file
