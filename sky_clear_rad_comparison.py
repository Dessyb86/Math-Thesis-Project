#########################################################################################################
# File Name:        sky_clear_rad_comparison                                                            #
# Editor:           Patrick Barrett                                                                     #
# Last update:      September 11, 2018                                                                  #
# Descrition :      This script reads a csv file of weather variables and outputs the values of the     #
#                   and variance of sky clear global radiation obtain by the two estimation methods     #
#                   specified                                                                           #
#                                                                                                       #
# input file:       The input file must be set to the variable "weather_data_file". The input file must #
#                   must contain columns of the following variables: (1) year (2) month (3) day (4) hour#
#                   (5) wind speed (kts) at the reference height (6) wind direction (tens of degrees) at#
#                   the reference height for wind (7) station level pressure (Hpa) (8) total cloud cover#
#                   (okta) (9) relative humidity (%)                                                    #
#                                                                                                       #
# output file:      The output file is entitled "sky_clear_rad.txt"                                      #
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

z_0 = 0.04 #roughness length for the area in metres

lat = 13.1919 #latitude of the location

albedo_0 = 0.165 #mean albedo at solar noon

cloud1_col = 17 #number of the columns with cloud information

cloud2_col = 18 #number of the columns with cloud information

cloud3_col = 19 #number of the columns with cloud information

cloud4_col = 20 #number of the columns with cloud information

B_10 = 0.1310428 #constant coefficient in calculating atmospheric transmittance (Hotel, 1976)

B_11 = 0.7342016 #coefficient of the exponential in calculating atmospheric transmittance (Hotel, 1976)

k_t = 0.1368551 #transmittance coefficient in calculating atmospheric transmittance (Hotel, 1976)
############################# DO NOT TAMPER BEYOND THIS POINT #####################

import numpy as np
import math
import csv
import calendar
import scipy
import scipy.integrate
from math import cos, sin, acos as arccos, asin as arcsin, tan as tg, degrees, radians

from scipy.special import lambertw
import os
import pandas as pd
from datetime import datetime
from datetime import timedelta
from astral import Astral

import sys

import patcon
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

########################## Defining Storage Variables #############################

skrad_diff_stack = [] #stack of the running sum of the sensible surface heatflux in the CBL

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
        
    tot_cloud[i] /= 8 #converts total cloud from oktas to fractional form

#################################### Defining Global Constants ############################

phi_lat = lat * np.pi / 180 #latitudinal angle of the location

omega = 0.0000729 #angular velocity of the earth
                
f = 2 * omega * np.sin(phi_lat) #coriolis parameter for the tropics (Blackadar and Tennekes, 1968) 

rho = 1.225 #density of the standard atmosphere

cp = 1006.4 #specific heat capacity of dry air at constant pressure

kv = 0.4 #vonKarman constant

g_0 = 9.80665 #standard gravity

heat_indicator = 0 #indicator that controls the reset of the heatflux_sum_stack and heatflux_time_stack variables

############################ Calculations#######################################           
for i in range(len(year)):

    quality_indicator = 0 #indicator if the input data for a particular instance is in good standing

    if ref_temp[i] >= 0:

        quality_indicator += 1

    if tot_cloud[i] >= 0 and tot_cloud[i] <= 1:

        quality_indicator += 1

    if quality_indicator == 2: #setting weather parameters to nan as there are errors in the key parameters

        ##########determining the phase of the PBL (pbl_indocator) ##########

        #calculating sunrise and sunset time
        date_GMT = datetime.strptime(str(int(year[i]))+'-'+str(int(month[i]))+'-'+str(int(day[i]))+':'+str(int(hour[i])),'%Y-%m-%d:%H') #datetime for the ith instance

        date_LST = date_GMT - timedelta(hours = 4) #changing time from UTC to LST


        city_name = 'Bridgetown' 
        a = Astral()
        a.solar_depression = 'civil'
        city = a[city_name]
        timezone = city.timezone

        date2 = datetime.strptime(str(int(date_LST.year))+'-'+str(int(date_LST.month))+'-'+str(int(date_LST.day)),'%Y-%m-%d')
        
        sun = city.sun(date2, local=True)

        lat = city.latitude * (np.pi / 180) #changing latitudinal angle from degrees to radians

        sunrise = sun['sunrise'].hour + (sun['sunrise'].minute / 60) + (sun['sunrise'].second / 3600)
        
        sunset = sun['sunset'].hour + (sun['sunset'].minute / 60) + (sun['sunset'].second / 3600)
        
        if date_LST.hour >= sunrise and date_LST.hour <= sunset: #current instance is within the daytime

            ############### calculating sensible surface heat flux (Phi_S) #################

            n_j = float(date_LST.strftime('%j')) #day number of the year for the current instance

            phi_D = ((-23.45*np.pi)/180.0)* np.cos(2*np.pi*((n_j + 10.0)/365.25)) #declination angle

            phi_H_current = np.pi*(float(date_LST.hour - sunrise)/float(sunset - sunrise))-((np.pi)/2) #hour angle of the current instance

            phi_Z_current = float(np.arccos(np.sin(phi_lat) * np.sin(phi_D) + np.cos(phi_lat) * np.cos(phi_D) * np.cos(phi_H_current))) #zenith angle of the current instance

##            albedo1 = albedo_0 + (1 - albedo_0) * np.exp(((-18 * phi_Z_previous) / np.pi) - 0.5 * (1 - albedo_0)**2) #albedo at current instance
##
##            albedo2 = albedo_0 + (1 - albedo_0) * np.exp(((-18 * phi_Z_current) / np.pi) - 0.5 * (1 - albedo_0)**2) #albedo oe hour prior to current instance 
##
##            mean_albedo = (albedo1 + albedo2)/2 #mean albedo

            Phi_SC_1 = 990.0 * (np.cos(phi_Z_current)) - 30.0 #sky clear global radiation
               
            if phi_Z_current == np.pi / 2:#calculating tb when zenith angle is pi/2

                tb = B_10

            else: #calculating sky clear radiation

                tb = B_10 + B_11 * np.exp(-k_t / np.cos(phi_Z_current)) #atmospheric transmittance for direct beam

            td = 0.271 - 0.294 * tb #atmospheric transmittance for diffused beam

            Q_on = 1367.0 * (1.0+(7.0/1367.0)*np.cos((2.0*np.pi*(n_j-3))/365.25)) #value of solar irradiance at the top of the atmosphere

            Qb = Q_on * tb #direct beam portion of the global solar radiation

            Qd = Q_on * td * np.cos(phi_Z_current) #diffused beam portion of the global solar radiation
            
            Phi_SC_2 = Qb + Qd #sky clear global solar radiation

            if Phi_SC_1 > 0:

                skrad_diff_stack.append(Phi_SC_2/Phi_SC_1)

            
output_file = open('sky_clear_rad.txt',"w")

output_file.write('The alternative method estimates the sky clear global radiation over the aermod method by a factor of {} +- {}'.format(np.mean(skrad_diff_stack),np.std(skrad_diff_stack)))

output_file.close()
  
