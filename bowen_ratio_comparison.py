#########################################################################################################
# File Name:        bowen_ratio_comparison                                                              #
# Editor:           Patrick Barrett                                                                     #
# Last update:      September 10, 2018                                                                  #
# Descrition :      This script reads a csv file of weather variables and outputs statistical data      #                                               #
#                   comparing the estimation of the bowen ratio as given by aermod vs Lin et al., 2016  #
#                                                                                                       #
# input file:       The input file must be set to the variable "weather_data_file". The input file must #
#                   must contain columns of the following variables: (1) year (2) month (3) day (4) hour#
#                   (5) wind speed (kts) at the reference height (6) wind direction (tens of degrees) at#
#                   the reference height for wind (7) station level pressure (Hpa) (8) total cloud cover#
#                   (okta) (9) relative humidity (%)                                                    #
#                                                                                                       #
# output file:      The output file is entitled "bowen_comparison.txt"                                  #
#########################################################################################################

############################# Mandatory User Defined Parameters ##############################

weather_data_file = 'new_data.csv' #name of the input file with the weather dat

year_col = 1 #number of the column holding the year of each instance

month_col = 2 #number of the column holding the month of each instance

day_col = 3 #number of the column holding the day of each instance

hour_col = 4 #number of the column holding the hour of each instance

tot_cloud_col = 6 #number of the column holding the total cloud cover of each instance

wind_dir_col = 7 #number of the column holding the wind direction of each instance

wind_spd_col = 8 #number of the column holding the wind speed of each instance

ref_temp_col = 9 #number of the column holding the reference temperature of each instance

rel_hum_col = 11 #number of the column holding the relative humidity of each instance

stn_pres_col = 12 #number of the column holding the station level pressure of each instance


############################# User Defined Parameters ##############################

rainy_season = [5,6,7,8,9,10,11] #array of the number of the months in the rainy season

bowen_rainy = 0.6 #bowen ratio value for rainy season

bowen_dry = 2.5 #bowen ratio value for dry season

############################# DO NOT TAMPER BEYOND THIS POINT #####################

import numpy as np
import math
import csv
import calendar
import scipy
import scipy.integrate
import matplotlib.pyplot as plt

import scipy.stats

from math import cos, sin, acos as arccos, asin as arcsin, tan as tg, degrees, radians

from scipy.special import lambertw
import os
import pandas as pd
from datetime import datetime
from datetime import timedelta
from astral import Astral

import sys

############################## Reading input data file #############################

csv1 = np.genfromtxt (weather_data_file, delimiter=",")

year = csv1[:,year_col] #the year the measurement was taken

month = csv1[:,month_col] #the month the measurement was taken

day = csv1[:,day_col] #the day the measurement was taken

hour = csv1[:,hour_col] #the hour (in UTC) the measurement was taken

tot_cloud = csv1[:,tot_cloud_col] #total cloud (oktas)

wind_dir = csv1[:,wind_dir_col] #wind direction (tens of degrees)

wind_spd = csv1[:,wind_spd_col] #wind speed (knots)

ref_temp = csv1[:,ref_temp_col] #reference air temperature (Celcius)

rel_hum = csv1[:,rel_hum_col] #relative humidity (percent)

pres = csv1[:,stn_pres_col] #atmospheric pressure measured at the height of the station (millibars)

########################### compilation variables #####################################

date_stack = [] #stack of date values

dry_temp_stack = [] #stack of temperature values during dry season

dry_rh_stack = [] #stack of relative humidity values during dry season

wet_temp_stack = [] #stack of temperature values during wet season

wet_rh_stack = [] #stack of relative humidity values during wet season

B_AERMOD_stack = [] #stack of bowen ratio values computed by aermod formulation

B_ALT_stack = [] #stack of bowen ratio values computed by alternate formulation (Lin et. al, 2016)

dry_stack = [] #stack of the calculated bowen ratio values for the dry season

rainy_stack = [] #stack of the calculated bowen ratio values for the rainy season

jan_stack = [] #stack of the calculated bowen ratio for the month of january

feb_stack = [] #stack of the calculated bowen ratio for the month of february

mar_stack = [] #stack of the calculated bowen ratio for the month of march

apr_stack = [] #stack of the calculated bowen ratio for the month of april

may_stack = [] #stack of the calculated bowen ratio for the month of may

jun_stack = [] #stack of the calculated bowen ratio for the month of june

jul_stack = [] #stack of the calculated bowen ratio for the month of july

aug_stack = [] #stack of the calculated bowen ratio for the month of august

sep_stack = [] #stack of the calculated bowen ratio for the month of september

oct_stack = [] #stack of the calculated bowen ratio for the month of october

nov_stack = [] #stack of the calculated bowen ratio for the month of november

dec_stack = [] #stack of the calculated bowen ratio for the month of december

overall_stack = [] #stack of all calculated values of bowen ratio

n_factor_dry = [] #stack of the coefficient of net radiation in calculating sensible surface heat flux for dry season

n_factor_wet = [] #stack of the coefficient of net radiation in calculating sensible surface heat flux for wet season

n_factor_overall = [] #stack of the coefficient of net radiation in calculating sensible surface heat flux overall

#################################### Output File ###############################

outfile = "bowen_comparison.txt"

dryfile = "bowen_dry_data.txt"

wetfile = "bowen_wet_data.txt"

#################################### Data Conversions ############################

for i in range(len(year)):

    #if-else statement calculates the angle between the wind vector and positive x-axis
    if wind_dir[i] >= 1 and wind_dir[i] <= 27:
        
        wind_dir[i] = 270 - 10 * wind_dir[i]

    elif wind_dir[i] >= 28 and wind_dir[i] <= 36:
        
        wind_dir[i] = 630 - 10 * wind_dir[i]

    else:
        
        wind_dir[i] = float('nan')

    wind_spd[i] *= 0.54444444444 #converts wind from knots to metres per second

    ref_temp[i] += 273.16 #converts reference temperature from Celcius to Kelvin
        
    rel_hum[i] /= 100 # converts relative humidty from percentage to decimal

    tot_cloud[i] /= 8 #converts total cloud from oktas to fractional form

#################################### Defining Global Constants ############################

date_GMT = datetime.strptime(str(int(year[i]))+'-'+str(int(month[i]))+'-'+str(int(day[i]))+':'+str(int(hour[i])),'%Y-%m-%d:%H') #datetime for the ith instance

date_LST = date_GMT - timedelta(hours = 4) #changing time from UTC to LST

################################## Calculations #########################################

for i in range(len(year)):

    B_0 = 1.46 * (1.0/rel_hum[i]) * ((ref_temp[i] / 273.0)**2) * (np.exp(-19.83 * (1.0 - (273.0 / ref_temp[i])))) #bowen ratio

    overall_stack.append(B_0)

    n_factor_overall.append(0.9 / (1 + (1 / B_0)))
    
    indicator = 1 #indicates dry season

    for j in range(len(rainy_season)): #checking if the instance occurs within the rainy season

        if month[i] == rainy_season[j]:

            indicator = 0 #indicates rainy season

    if indicator == 1: #checking if the instance occurs in the dry season

        dry_stack.append(B_0) #appending to stack for the dry season

        n_factor_dry.append(0.9 / (1 + (1 / B_0)))

        dry_temp_stack.append(ref_temp[i])

        dry_rh_stack.append(rel_hum[i])

        B_AERMOD_stack.append(2.0)

        B_ALT_stack.append(B_0)

    else:

        rainy_stack.append(B_0) #appending to stack for the rainy season

        n_factor_wet.append(0.9 / (1 + (1 / B_0)))

        wet_temp_stack.append(ref_temp[i])

        wet_rh_stack.append(rel_hum[i])

        B_AERMOD_stack.append(0.6)

        B_ALT_stack.append(B_0)

    if month[i] == 1:

        jan_stack.append(B_0)

    elif month[i] == 2:

        feb_stack.append(B_0)

    elif month[i] == 3:

        mar_stack.append(B_0)

    elif month[i] == 4:

        apr_stack.append(B_0)

    elif month[i] == 5:

        may_stack.append(B_0)

    elif month[i] == 6:

        jun_stack.append(B_0)

    elif month[i] == 7:

        jul_stack.append(B_0)

    elif month[i] == 8:

        aug_stack.append(B_0)

    elif month[i] == 9:

        sep_stack.append(B_0)

    elif month[i] == 10:

        oct_stack.append(B_0)

    elif month[i] == 11:

        nov_stack.append(B_0)

    elif month[i] == 12:

        dec_stack.append(B_0)

    else:
        print("error")
        
output_file = open(outfile,"w")

output_file.write('The following are mean values for the bowen ratio obtained by the method described by Lin et al., 2016 \n\n')

output_file.write('Jan - {} - {}\n'.format(np.mean(jan_stack),np.std(jan_stack)))

output_file.write('Feb - {} - {}\n'.format(np.mean(feb_stack),np.std(feb_stack)))

output_file.write('Mar - {} - {}\n'.format(np.mean(mar_stack),np.std(mar_stack)))

output_file.write('Apr - {} - {}\n'.format(np.mean(apr_stack),np.std(apr_stack)))

output_file.write('May - {} - {}\n'.format(np.mean(may_stack),np.std(may_stack)))

output_file.write('Jun - {} - {}\n'.format(np.mean(jun_stack),np.std(jun_stack)))

output_file.write('Jul - {} - {}\n'.format(np.mean(jul_stack),np.std(jul_stack)))

output_file.write('Aug - {} - {}\n'.format(np.mean(aug_stack),np.std(aug_stack)))

output_file.write('Sep - {} - {}\n'.format(np.mean(sep_stack),np.std(sep_stack)))

output_file.write('Oct - {} - {}\n'.format(np.mean(oct_stack),np.std(oct_stack)))

output_file.write('Nov - {} - {}\n'.format(np.mean(nov_stack),np.std(nov_stack)))

output_file.write('Dec - {} - {}\n\n'.format(np.mean(dec_stack),np.std(dec_stack)))

output_file.write('Dry Season - {} - {}\n'.format(np.mean(dry_stack),np.std(dry_stack)))

output_file.write('Rainy Season - {} - {}\n'.format(np.mean(rainy_stack),np.std(rainy_stack)))

output_file.write('Annual Mean - {} - {}\n'.format(np.mean(overall_stack),np.std(overall_stack)))

output_file.write('N-Factor Dry - {} - {}\n'.format(np.mean(n_factor_dry),np.std(n_factor_dry)))
                  
output_file.write('N-Factor Wet - {} - {}\n'.format(np.mean(n_factor_wet),np.std(n_factor_wet)))

output_file.write('N-Factor Overall - {} - {}\n'.format(np.mean(n_factor_overall),np.std(n_factor_overall)))
                  
output_file.close()

print(len(dry_stack))

print(len(rainy_stack))

xlist = np.linspace(0.0001, 1.0, 100)
ylist = np.linspace(283,313, 100)
X,Y = np.meshgrid(xlist, ylist)
Z = 1.46*(1/X)*(Y/273)**2*np.exp(-19.83*(1-273/Y))

CS1 = plt.contour(X, Y, Z, [0.6], colors = ['blue'], linestyles = 'solid')
plt.clabel(CS1, inline=1, fontsize=10)
plt.scatter(wet_rh_stack, wet_temp_stack, c='blue')
CS2 = plt.contour(X, Y, Z, [0.3526], colors = ['black'], linestyles = 'dashed')
plt.clabel(CS2, inline=1, fontsize=10)
plt.xlabel('Relative Humidity (dimensionless)', fontsize=12)
plt.ylabel('Absolute Temperature (Kelvin)', fontsize=12)

plt.savefig('bowen_plot_wet.png')
plt.clf()
plt.cla()
plt.close()


CS3 = plt.contour(X, Y, Z, [2.0], colors = ['red'], linestyles = 'solid')
plt.clabel(CS3, inline=1, fontsize=10)
plt.scatter(dry_rh_stack, dry_temp_stack, c='red')
CS4 = plt.contour(X, Y, Z, [0.4007], colors = ['black'], linestyles = 'dotted')
plt.clabel(CS4, inline=1, fontsize=10)
plt.xlabel('Relative Humidity (dimensionless)', fontsize=12)
plt.ylabel('Absolute Temperature (Kelvin)', fontsize=12)

plt.savefig('bowen_plot_dry.png')
plt.clf()
plt.cla()
plt.close()

CS5 = plt.contour(X, Y, Z, [0.6], colors = ['blue'], linestyles = 'solid')
plt.clabel(CS5, inline=1, fontsize=10)
plt.scatter(wet_rh_stack, wet_temp_stack, c='blue')
CS6 = plt.contour(X, Y, Z, [2.0], colors = ['red'], linestyles = 'solid')
plt.clabel(CS6, inline=1, fontsize=10)
plt.scatter(dry_rh_stack, dry_temp_stack, c='red')
CS7 = plt.contour(X, Y, Z, [0.3526], colors = ['black'], linestyles = 'dashed')
plt.clabel(CS7, inline=1, fontsize=10)
CS8 = plt.contour(X, Y, Z, [0.4007], colors = ['black'], linestyles = 'dotted')
plt.clabel(CS8, inline=1, fontsize=10)
plt.xlabel('Relative Humidity (dimensionless)', fontsize=12)
plt.ylabel('Absolute Temperature (Kelvin)', fontsize=12)


plt.savefig('bowen_plot_both.png')
plt.clf()
plt.cla()
plt.close()

xlist = np.linspace(0.0001, 1.0, 100)
ylist = np.linspace(260,350, 100)
X,Y = np.meshgrid(xlist, ylist)
Z = 1.46*(1/X)*(Y/273)**2*np.exp(-19.83*(1-273/Y))

CS = plt.contour(X, Y, Z, [0.6,2.0], colors = ['blue','red'], linestyles = 'solid')
plt.clabel(CS, inline=1, fontsize=10)
plt.savefig('bowen_plot_big.png')
plt.clf()
plt.cla()
plt.close()

print(scipy.stats.kruskal(dry_stack,rainy_stack))
H,p = scipy.stats.kruskal(dry_stack,rainy_stack)

print(H)
print(p)

