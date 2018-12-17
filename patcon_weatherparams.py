#########################################################################################################
# File Name:        patcon_weatherparams                                                                #
# Editor:           Patrick Barrett                                                                     #
# Last update:      September 10, 2018                                                                  #
# Descrition :      This script reads a csv file of weather variables and outputs a csv file containing #
#                   the weather parameters required to estimate model parameters of the Gaussian plume  #
#                   dispersion model. Additionally, the user must specify a group of indicators and     #
#                   variables to guide the calculations.                                                #
#                                                                                                       #
# input file:       The input file must be set to the variable "weather_data_file". The input file must #
#                   must contain columns of the following variables: (1) year (2) month (3) day (4) hour#
#                   (5) wind speed (kts) at the reference height (6) wind direction (tens of degrees) at#
#                   the reference height for wind (7) station level pressure (Hpa) (8) total cloud cover#
#                   (okta) (9) relative humidity (%)                                                    #
#                                                                                                       #
# output file:      The output file is entitled "weatherparams_i_j_k_l, where i,j,k and l are the user  #
#                   defined indicators to calculate convective mixing layer height, bowen ratio, sky    #
#                   clear global radiation and the attenuation factor, respectfully. Another output file#
#                   entitled "heatflux_comparison_j_k_l" with values of the bowen ratio, attenuation    #
#                   factor, clear sky global solar radiation and associated heatflux will be printed    #
#########################################################################################################

############################# Mandatory User Defined Parameters ##############################

weather_data_file = 'dataset_1.csv' #name of the input file with the weather dat

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

z_uref = 10 #height of the wind speed measurements in metres

z_tref = 1.5 #height of the temperature measurements in metres

z_station = 77 #height of the measure station in metres

lat = 13.1919 #latitude of the location

albedo_0 = 0.165 #mean albedo at solar noon

############################# Conditional User Defined Parameters ##############################

############## CBL Parameters (If CBL calculations are involved)###############

#The following conditions define the formula with which the heatflux for the CBL should be calculated

B_0_indicator = 1 #indicates what formula should be used to calculate the Bowen ratio. 

#B_0_indicator = 1 - Bowen ratio will be calculated by aermod formulation (using tables of mean seasonal values)

#B_0_indicator = 2 - Bowen ratio will be calculated from reference temperature and relative humidity (Lin et al., 2016)

R_SC_indicator = 1 # indicates what formula should be used to calculate the sky clear global radiation. 

#R_SC_indicator = 1 - sky clear radiation will be calculated by aermod formulation (Kasten and Czeplak, 1980)

#R_SC_indicator = 2 - sky clear global radiation will be calculated from julian day and zenith angle (Khan and Ahmad, 2012)

A_indicator = 1 # indicates what formula should be used to calculate the attenuation factor. 

#A_indicator = 1 - attenuation factor will be calculated by aermod formulation (Kasten and Czeplak, 1980)

#A_indicator = 2 - attenuation factor will be calculated by weight formula using cloud cover and transmittance for each layer from Katen and Czeplak, 1980

#A_indicator = 3 - attenuation factor will be calculated by weight formula using cloud cover and transmittance for each layer from Walled and Yussra, 2006

#The following must be specified if A_indicator was set to either 2 or 3

cloud1_col = 17 #number of the columns with cloud information

cloud2_col = 18 #number of the columns with cloud information

cloud3_col = 19 #number of the columns with cloud information

cloud4_col = 20 #number of the columns with cloud information

B_10 = 0.1310428 #constant coefficient in calculating atmospheric transmittance (Hotel, 1976)

B_11 = 0.7342016 #coefficient of the exponential in calculating atmospheric transmittance (Hotel, 1976)

k_t = 0.1368551 #transmittance coefficient in calculating atmospheric transmittance (Hotel, 1976)

#The following must be specified if B_0_indicator  was set to 1

rainy_season = [5,6,7,8,9,10,11] #array of the number of the months in the rainy season

bowen_rainy = 0.6 #bowen ratio value for rainy season

bowen_dry = 2.5 #bowen ratio value for dry season

zic_indicator = 2 #indicates what formula will be used to calculate the convective mixing layer height in the CBL

#zic_indicator = 1 - convective mixing height will be calculated by the temporal evolution (Carson, 1973; Weil and Brower, 1983)

#zic_indicator = 2 - convective mixing height will be calculated as the exact expression of the lifting condesation level (Romps, 2017)

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

rel_hum = csv1[:,rel_hum_col] #relative humidity (percent)

pres = csv1[:,stn_pres_col] #atmospheric pressure measured at the height of the station (millibars)

if A_indicator == 2 or A_indicator == 3:
    
    cloud1 = csv1[:,cloud1_col] #first cloud group

    cloud2 = csv1[:,cloud2_col] #second cloud group

    cloud3 = csv1[:,cloud3_col] #third cloud group

    cloud4 = csv1[:,cloud4_col] #fourth cloud group

########################## Defining Storage Variables #############################

a_stack =[] #stack of cloud transmittance values

B0_stack =[] #stack of bowen ratio values

SC_stack =[] #stack of clear sky global radiation values

z_im_stack =[] #stack of mechanical mixing height values

z_ic_stack = [] #stack of convective mixing height values

zi_stack = [] #stack of mixing layer height values

heatflux_stack =[] #stack of heatflux values

L_stack =[] #stack of obukhov length values

sfv_stack =[] #stack of shear velocity values

pts_stack =[] #stack of potential temperature scale values

w_stack = [] #stack of Deardorff values

date_stack = [] #stack of date values

heatflux_sum_stack = [] #stack of the running sum of the sensible surface heatflux in the CBL

heatflux_time_stack = [] #stack of the time elapsed in the CBL for the particular day
#################################### Output File ###############################

outfile = "weatherparams_"+str(int(zic_indicator))+"_"+str(int(B_0_indicator))+"_"+str(int(R_SC_indicator))+"_"+str(A_indicator)

outfile2 = "heatflux_comparison_"+str(int(B_0_indicator))+"_"+str(int(R_SC_indicator))+"_"+str(A_indicator)
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

date1 = datetime.strptime(str(2000)+'-'+str(9)+'-'+str(17)+':'+str(11),'%Y-%m-%d:%H') #arbitrary datetime variable (2000-11-17 11:00:00)
date2 = datetime.strptime(str(2000)+'-'+str(9)+'-'+str(17)+':'+str(12),'%Y-%m-%d:%H') #arbitrary datetime variable (2000-11-17 12:00:00)

time_check = date2-date1 #stores the value of one hour to be used in comparisons

phi_lat = lat * np.pi / 180 #latitudinal angle of the location

omega = 0.0000729 #angular velocity of the earth
                
f = 2 * omega * np.sin(phi_lat) #coriolis parameter for the tropics (Blackadar and Tennekes, 1968) 

rho = 1.225 #density of the standard atmosphere

cp = 1006.4 #specific heat capacity of dry air at constant pressure

kv = 0.4 #vonKarman constant

g_0 = 9.80665 #standard gravity

heat_indicator = 0 #indicator that controls the reset of the heatflux_sum_stack and heatflux_time_stack variables

################################## checking rainy_season array #################

if B_0_indicator == 1:
    
    if len(rainy_season) == 0:

        print("rainy season array is empty")

        sys.exit()

    for i in range(len(rainy_season)):

        if rainy_season[i] < 1 or rainy_season[i] > 12 or (rainy_season[i] % 1 != 0):

            print("month values in rainy_season array must be whole numbers between 1 and 12 inclusive")

            sys.exit()
            
############################ weather parameter functions ###################

def transmittance1(cloud_info):#calculates the transmittance of the cloud layer for A_indicator = 2

    x = int((cloud_info % 10000) / 1000)
    
    if x == 0:#tramsmissivity for cirrus
        
        x = 0.61
        
    elif x == 1:#tramsmissivity for cirrocumulus
        
        x = 0.61
        
    elif x== 2:#tramsmissivity for cirrostratus
        
        x = 0.61
        
    elif x == 3:#tramsmissivity for altocumulus
        
        x = 0.27
        
    elif x == 4:#tramsmissivity for altostratus

        x = 0.27

    elif cloud_info == 5:#tramsmissivity for nimbostratus

        x = 0.26

    elif x == 6:#tramsmissivity for stratocumulus

        x = 0.25

    elif x == 7:#tramsmissivity for stratus

        x=0.18

    elif x == 8:#tramsmissivity for cumulus

        x=0.25

    elif x == 9:#tramsmissivity for cumulonimbus

        x=0.25

    else:#unknown cloud type

        x= float('nan')
    
    return x

def transmittance2(cloud_info):

    x = int((cloud_info % 10000) / 1000)
    
    if x == 0:#tramsmissivity for cirrus
        
        x = 0.82
        
    elif x == 1:#tramsmissivity for cirrocumulus
        
        x = 0.84
        
    elif x== 2:#tramsmissivity for cirrostratus
        
        x = 0.82
        
    elif x == 3:#tramsmissivity for altocumulus
        
        x = 0.49
        
    elif x == 4:#tramsmissivity for altostratus

        x = 0.50

    elif cloud_info == 5:#tramsmissivity for nimbostratus

        x = 0.55

    elif x == 6:#tramsmissivity for stratocumulus

        x = 0.21

    elif x == 7:#tramsmissivity for stratus

        x = float('nan')

    elif x == 8:#tramsmissivity for cumulus

        x=0.19

    elif x == 9:#tramsmissivity for cumulonimbus

        x=0.19

    else:#unknown cloud type

        x= float('nan')
    
    return x

def windspd(zr,z_0,z_i,L,SFV,pbl_indicator,kv):

    if pbl_indicator == 1:
        
        m_0 = lambda z: (1 - 16 * (z_0 / L))**0.25

        m_u = lambda z: (1 - 16 * (z / L))**0.25

        w_u = lambda z: 2 * np.log((1 + m_u(z)) / 2) + np.log((1 + (m_u(z))**2) / 2) - 2 * np.arctan(m_u(z)) + (np.pi / 2) #universal stability function evaluated at the reference height for wind speed measurements

        w_0 = lambda z: 2 * np.log((1 + m_0(z)) / 2) + np.log((1 + (m_0(z))**2) / 2) - 2 * np.arctan(m_0(z)) + (np.pi / 2) #universal stability function evaluated at the roughness length

        u = lambda z: (SFV / kv) * (np.log(z / z_0) - w_u(z) + w_0(z)) # general formulation for wind speed according to MOS theory

        u_7 = u(7 * z_0) # wind speed at 0.7 * z_0

        u_zi = u(z_i)

        meanspd = lambda z: 0 if z < z_0 else (z / (7 * z_0)) * u_7 if z >= z_0 and z <= 7 * z_0 else u(z) if z > 7 * z_0 and z <= z_i else u_zi  #wind speed with adjustment according to Cenedese et al., 1998)

    else:
        u = lambda z: (SFV / kv) * (np.log(z / z_0) + (5 * z / L))

        u_7 = u(7 * z_0) # wind speed at 0.7 * z_0

        u_zi = u(z_i)

        meanspd = lambda z: 0 if z < z_0 else (z / (7 * z_0)) * u_7 if z >= z_0 and z <= 7 * z_0 else u(z) if z > 7 * z_0 and z <= z_i else u_zi  #wind speed with adjustment according to Cenedese et al., 1998)

    return meanspd(zr)

def lcl(p,T,rh=None,rhl=None,rhs=None,return_ldl=False,return_min_lcl_ldl=False):

   import math
   import scipy.special

   # Parameters
   Ttrip = 273.16     # K
   ptrip = 611.65     # Pa
   E0v   = 2.3740e6   # J/kg
   E0s   = 0.3337e6   # J/kg
   ggr   = 9.81       # m/s^2
   rgasa = 287.04     # J/kg/K 
   rgasv = 461        # J/kg/K 
   cva   = 719        # J/kg/K
   cvv   = 1418       # J/kg/K 
   cvl   = 4119       # J/kg/K 
   cvs   = 1861       # J/kg/K 
   cpa   = cva + rgasa
   cpv   = cvv + rgasv

   # The saturation vapor pressure over liquid water
   def pvstarl(T):
      return ptrip * (T/Ttrip)**((cpv-cvl)/rgasv) * \
         np.exp( (E0v - (cvv-cvl)*Ttrip) / rgasv * (1/Ttrip - 1/T) )
   
   # The saturation vapor pressure over solid ice
   def pvstars(T):
      return ptrip * (T/Ttrip)**((cpv-cvs)/rgasv) * \
         np.exp( (E0v + E0s - (cvv-cvs)*Ttrip) / rgasv * (1/Ttrip - 1/T) )

   # Calculate pv from rh, rhl, or rhs
   rh_counter = 0
   if rh  is not None:
      rh_counter = rh_counter + 1
   if rhl is not None:
      rh_counter = rh_counter + 1
   if rhs is not None:
      rh_counter = rh_counter + 1
   if rh_counter != 1:
      print(rh_counter)
      exit('Error in lcl: Exactly one of rh, rhl, and rhs must be specified')
   if rh is not None:
      # The variable rh is assumed to be 
      # with respect to liquid if T > Ttrip and 
      # with respect to solid if T < Ttrip
      if T > Ttrip:
         pv = rh * pvstarl(T)
      else:
         pv = rh * pvstars(T)
      rhl = pv / pvstarl(T)
      rhs = pv / pvstars(T)
   elif rhl is not None:
      pv = rhl * pvstarl(T)
      rhs = pv / pvstars(T)
      if T > Ttrip:
         rh = rhl
      else:
         rh = rhs
   elif rhs is not None:
      pv = rhs * pvstars(T)
      rhl = pv / pvstarl(T)
      if T > Ttrip:
         rh = rhl
      else:
         rh = rhs
   if pv > p:
      return float('nan')

   # Calculate lcl_liquid and lcl_solid
   qv = rgasa*pv / (rgasv*p + (rgasa-rgasv)*pv)
   rgasm = (1-qv)*rgasa + qv*rgasv
   cpm = (1-qv)*cpa + qv*cpv
   if rh == 0:
      return cpm*T/ggr
   aL = -(cpv-cvl)/rgasv + cpm/rgasm
   bL = -(E0v-(cvv-cvl)*Ttrip)/(rgasv*T)
   cL = pv/pvstarl(T)*np.exp(-(E0v-(cvv-cvl)*Ttrip)/(rgasv*T))
   aS = -(cpv-cvs)/rgasv + cpm/rgasm
   bS = -(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*T)
   cS = pv/pvstars(T)*np.exp(-(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*T))
   lcl = cpm*T/ggr*( 1 - \
      bL/(aL*scipy.special.lambertw(bL/aL*cL**(1/aL),-1).real) )
   ldl = cpm*T/ggr*( 1 - \
      bS/(aS*scipy.special.lambertw(bS/aS*cS**(1/aS),-1).real) )

   # Return either lcl or ldl
   if return_ldl and return_min_lcl_ldl:
      exit('return_ldl and return_min_lcl_ldl cannot both be true')
   elif return_ldl:
      return ldl
   elif return_min_lcl_ldl:
      return min(lcl,ldl)
   else:
      return lcl

def integrand(z):
    
    return theta_0 + (PTS / kv) * (np.log(z / z_0) - (2 * np.log((1 + (1 - 16 * z / L)**0.5) / 2)+np.log((1 + (1 - 16 * z / L)) / 2) - np.arctan((1 - 16 * z / L)**0.5) + np.pi / 2) + (2 * np.log((1 + (1 - 16 * z_0 / L)**0.5) / 2)+np.log((1 + (1 - 16 * z_0 / L)) / 2) - np.arctan((1 - 16 * z_0 / L)**0.5) + np.pi / 2))

def func(z):
   
    integral,err = scipy.integrate.quad(integrand, z_0, z)
    
    return RHS - (theta_0 + (PTS / kv) * (np.log(z / z_0) - (2 * np.log((1 + (1 - 16 * z / L)**0.5) / 2)+np.log((1 + (1 - 16 * z / L)) / 2) - np.arctan((1 - 16 * z / L)**0.5) + np.pi / 2) + (2 * np.log((1 + (1 - 16 * z_0 / L)**0.5) / 2)+np.log((1 + (1 - 16 * z_0 / L)) / 2) - np.arctan((1 - 16 * z_0 / L)**0.5) + np.pi / 2))) * z + integral

vfunc = np.vectorize(func)

z_guess = 850

############################ Calculations#######################################           
for i in range(len(year)):

    quality_indicator = 0 #indicator if the input data for a particular instance is in good standing
    
    if wind_spd[i] > 0:

        quality_indicator += 1

    if ref_temp[i] >= 0:

        quality_indicator += 1

    if rel_hum[i] >= 0 and rel_hum[i] <= 1:

        quality_indicator += 1

    if tot_cloud[i] >= 0 and tot_cloud[i] <= 1:

        quality_indicator += 1

    if A_indicator ==2 or A_indicator == 3:

        if cloud1[i] >=0 and cloud1[i] <= 99999:

            quality_indicator +=1

        if cloud2[i] >=0 and cloud2[i] <= 99999:

            quality_indicator +=1

        if cloud3[i] >=0 and cloud3[i] <= 99999:

            quality_indicator +=1

        if cloud4[i] >=0 and cloud4[i] <= 99999:

            quality_indicator +=1

        quality_indicator /= 2

    if quality_indicator != 4: #setting weather parameters to nan as there are errors in the key parameters
  
        z_im = float('nan')

        z_ic = float('nan')

        z_i = float('nan')

        Phi_S = float('nan')

        L = float('nan')

        SFV = float('nan')

        PTS = float('nan')

        w = float('nan')

        date_LST = float('nan')

        heat_indicator = 0

    else: #calculating model parameters as there are no errors in the key parameters
 
        ###########creating datetime variables (date_LST and time_indicator) #############

        date_GMT = datetime.strptime(str(int(year[i]))+'-'+str(int(month[i]))+'-'+str(int(day[i]))+':'+str(int(hour[i])),'%Y-%m-%d:%H') #datetime for the ith instance

        date_LST = date_GMT - timedelta(hours = 4) #changing time from UTC to LST


        if i > 0: #calculating date of previous row of data
        
            prev_date = datetime.strptime(str(int(year[i-1]))+'-'+str(int(month[i-1]))+'-'+str(int(day[i-1]))+':'+str(int(hour[i-1])),'%Y-%m-%d:%H') #datetime for the previous instance

        else: #setting the datetime variable for the previous instance

            prev_date = date_GMT #sets datetime variable for previous instace to present instance

        if (date_GMT - prev_date) == time_check:
            
            time_indicator = 1 #indicates the mechanical mixing layer height will be calculated by time evolution of the equivalent mixing layer height

        else:

            time_indicator = 0 #indicates the mechanical mixing layer height will be given as the equivalent mixing height

        ##########determining the phase of the PBL (pbl_indocator) ##########

        pbl_indicator = 0 #indicates the SBL phase

        #calculating sunrise and sunset time
    
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

            phi_H_previous = np.pi*(float((date_LST.hour - 1) - sunrise)/float(sunset - sunrise))-((np.pi)/2) #hour angle of the hour before the current instance

            phi_Z_current = float(np.arccos(np.sin(phi_lat) * np.sin(phi_D) + np.cos(phi_lat) * np.cos(phi_D) * np.cos(phi_H_current))) #zenith angle of the current instance
            
            phi_Z_previous = float(np.arccos(np.sin(phi_lat) * np.sin(phi_D) + np.cos(phi_lat) * np.cos(phi_D) * np.cos(phi_H_previous))) #zenith angle of the hour before the current instance

            albedo1 = albedo_0 + (1 - albedo_0) * np.exp(((-18 * phi_Z_previous) / np.pi) - 0.5 * (1 - albedo_0)**2) #albedo at current instance

            albedo2 = albedo_0 + (1 - albedo_0) * np.exp(((-18 * phi_Z_current) / np.pi) - 0.5 * (1 - albedo_0)**2) #albedo oe hour prior to current instance 

            mean_albedo = (albedo1 + albedo2)/2 #mean albedo

            if R_SC_indicator == 1: #calculating sky clear global radiation

                Phi_SC = 990.0 * (np.cos(phi_Z_current)) - 30.0 #sky clear global radiation

                SC_stack.append(Phi_SC)

            elif R_SC_indicator == 2: #calculating sky clear global radiation
                
                if phi_Z_current == np.pi / 2:#calculating tb when zenith angle is pi/2

                    tb = B_10

                else: #calculating sky clear radiation

                    tb = B_10 + B_11 * np.exp(-k_t / np.cos(phi_Z_current)) #atmospheric transmittance for direct beam

                td = 0.271 - 0.294 * tb #atmospheric transmittance for diffused beam

                Q_on = 1367.0 * (1.0+(7.0/1367.0)*np.cos((2.0*np.pi*(n_j-3))/365.25)) #value of solar irradiance at the top of the atmosphere

                Qb = Q_on * tb #direct beam portion of the global solar radiation

                Qd = Q_on * td * np.cos(phi_Z_current) #diffused beam portion of the global solar radiation
                
                Phi_SC = Qb + Qd #sky clear global solar radiation

                SC_stack.append(Phi_SC)

            else: #R_SC_indicator value erroneous
 
                print("Invalid value given for R_SC_indicator")

                sys.exit()

                SC_stack.append(Phi_SC)

            if B_0_indicator == 1:#calculating bowen ratio

                wet_month_indicator = 1 #indicates dry season

                for j in range(len(rainy_season)): #checking whether the month is in the rainy season or not

                    if date_LST.month == rainy_season[j]:

                        wet_month_indicator = 0 #indicates rainy season
                        
                if wet_month_indicator == 0:#current instance is in the wet season

                    B_0 = bowen_rainy #bowen ratio for wet season

                else:#current instance is in the dry season

                    B_0 = bowen_dry #bowen ratio for dry season

                B0_stack.append(B_0)
                    
            else: #calculating bowen ratio                

                B_0 = 1.46 * (1.0/rel_hum[i]) * ((ref_temp[i] / 273.0)**2) * (np.exp(-19.83 * (1.0 - (273.0 / ref_temp[i])))) #bowen ratio

                B0_stack.append(B_0)
            
            if A_indicator == 1:#calculating cloud attenuation factor

                A = 1 - 0.75 * (tot_cloud[i])**3.4

                a_stack.append(A)

            elif A_indicator == 2:

                c1 = (int(cloud1[i]/10000))/8.0 #fractional cloud cover for 1st cloud layer

                c2 = (int(cloud2[i]/10000))/8.0 #fractional cloud cover for 2nd cloud layer

                c3 = (int(cloud3[i]/10000))/8.0 #fractional cloud cover for 3rd cloud layer

                c4 = (int(cloud4[i]/10000))/8.0 #fractional cloud cover for 4th cloud layer

                k1 = transmittance1(cloud1[i]) #cloud transmittance for the 1st cloud layer

                k2 = transmittance1(cloud2[i]) #cloud transmittance for the 2nd cloud layer

                k3 = transmittance1(cloud3[i]) #cloud transmittance for the 3rd cloud layer

                k4 = transmittance1(cloud4[i]) #cloud transmittance for the 4th cloud layer

                A = (1 - tot_cloud[i]) + tot_cloud[i] * (1 - (1 - c1)* k1) * (1 - (1 - c2) * k2) * (1 - (1 - c3) * k3) * (1 - (1 - c4) * k4) #attenuation factor

                a_stack.append(A)
                
            else:
                
                c1 = (int(cloud1[i]/10000))/8.0 #fractional cloud cover for 1st cloud layer

                c2 = (int(cloud2[i]/10000))/8.0 #fractional cloud cover for 2nd cloud layer

                c3 = (int(cloud3[i]/10000))/8.0 #fractional cloud cover for 3rd cloud layer

                c4 = (int(cloud4[i]/10000))/8.0 #fractional cloud cover for 4th cloud layer

                k1 = transmittance2(cloud1[i]) #cloud transmittance for the 1st cloud layer

                k2 = transmittance2(cloud2[i]) #cloud transmittance for the 2nd cloud layer

                k3 = transmittance2(cloud3[i]) #cloud transmittance for the 3rd cloud layer

                k4 = transmittance2(cloud4[i]) #cloud transmittance for the 4th cloud layer

                A = (1 - tot_cloud[i]) + tot_cloud[i] * (1 - (1 - c1)* k1) * (1 - (1 - c2) * k2) * (1 - (1 - c3) * k3) * (1 - (1 - c4) * k4) #attenuation factor

                a_stack.append(A)
 
            Phi_SO = Phi_SC * A #global solar radiation

            Phi_N = ((1 - mean_albedo) * Phi_SO + (5.31 * 10**(-13)) * (ref_temp[i])**6 - (5.67 * 10**(-8)) * ref_temp[i]**4 + 60 * tot_cloud[i]) / (1 + 0.12) #net radiation

            Phi_S = (0.9 * Phi_N) / (1 + (1 / B_0))

            if Phi_S > 0:#setting pbl_indicator to indicate CBL phase

                pbl_indicator = 1 #PBL set to CBL phase

        if pbl_indicator == 1: #calculating weather parameters for CBL phase

            #calculating Obukhov length

            L_1 = (- rho * cp * ref_temp[i] * (kv**2) * (wind_spd[i])**3) / (g_0 * Phi_S * (np.log(z_uref / z_0))**3) #first iteration of Obukhov length, L_0 is assumed to be infinity

            m_u = (1 - 16 * (z_uref / L_1))**0.25 #stability function evaluated at the reference height for wind measurements
            
            m_0 = (1 - 16 * (z_0 / L_1))**0.25 #stability function evaluated at the roughness length
            
            w_u = 2 * np.log((1 + m_u) / 2) + np.log((1 + m_u**2) / 2) - 2 * np.arctan(m_u) + (np.pi / 2) #universal stability function evaluated at the reference height for wind speed measurements

            w_0 = 2 * np.log((1 + m_0) / 2) + np.log((1 + m_0**2) / 2) - 2 * np.arctan(m_0) + (np.pi / 2) #universal stability function evaluated at the roughness length

            L_2 = (- rho * cp * ref_temp[i] * (kv**2) * (wind_spd[i])**3) / (g_0 * Phi_S * (np.log(z_uref / z_0) - w_u + w_0)**3) #second iteration of Obukhov length, L_0 is assumed to be infinity
            
            stoppage_criteria = abs((L_2 - L_1) / L_1) #stoppage criteria for iteration

            L = L_2
     
            while stoppage_criteria > 0.00000001:

                m_u = (1 - 16 * (z_uref / L))**0.25 #stability function evaluated at the reference height for wind measurements
                
                m_0 = (1 - 16 * (z_0 / L))**0.25 #stability function evaluated at the roughness length
                
                w_u = 2 * np.log((1 + m_u) / 2) + np.log((1 + m_u**2) / 2) - 2 * np.arctan(m_u) + (np.pi / 2) #universal stability function evaluated at the reference height for wind speed measurements

                w_0 = 2 * np.log((1 + m_0) / 2) + np.log((1 + m_0**2) / 2) - 2 * np.arctan(m_0) + (np.pi / 2) #universal stability function evaluated at the roughness length
                   
                L_new = (- rho * cp * ref_temp[i] * (kv**2) * (wind_spd[i])**3) / (g_0 * Phi_S * (np.log(z_uref / z_0) - w_u + w_0)**3) #second iteration of Obukhov length, L_0 is assumed to be infinity

                stoppage_criteria = abs((L_new - L) / L) #stoppage criteria for iteration

                L = L_new #updating value of Obukhov length

            SFV = (kv * wind_spd[i]) / (np.log(z_uref / z_0) - w_u + w_0) #shear velocity

            PTS = -Phi_S / (rho * cp * SFV) #potential temperature scale

            if abs(SFV / (f * L)) <= 4:

                z_ie = 0.2 * SFV / f

            else:
                z_ie = 0.4 * abs((SFV * L) / f)

            #z_ie2 = -f * L + np.sqrt((f**2) * (L**2) + 0.57 * SFV * f * L)            
    
            if time_indicator == 0: #mechinal mixing height as the equivalent mixing height

                z_im = z_ie#mechanical mixing layer height set to equivalent mixing height

            elif np.isnan(z_im_stack[-1]):#no estimate of the equivalent mixing height for the previous instance

               z_im = z_ie #mechanical mixing layer height set to equivalent mixing height

            elif heatflux_stack[-1] <= 0 :#previous instance was in the SBL

                z_im = z_ie #mechanical mixing layer height set to equivalent mixing height
                
            else: #mechanical mixing height by time evolution equation

                tau = z_im_stack[-1] / (2 * SFV) #time scale

                z_im = z_im_stack[-1] * np.exp(-3600 / tau) + z_ie * (1 - np.exp(-3600 / tau)) #mechanical mixing height by time evolution (Venkatram 1982)

            if zic_indicator == 1:

                if date_LST.hour - sunrise > 0 and date_LST.hour - sunrise <=1: #first hour after sunrise

                    stoppage_criteria = len(heatflux_sum_stack) 

                    while stoppage_criteria > 0:

                        heatflux_sum_stack.pop(0)

                        heatflux_time_stack.pop(0)

                        stoppage_criteria = len(heatflux_sum_stack) 

                    heatflux_sum_stack.append(0) #appending the heatflux at the time of sunrise

                    heatflux_sum_stack.append(Phi_S) #appending the heatflux at the first hour after sunrise

                    heatflux_time_stack.append(0) #appending the time sunrise which corresponds to a heatflux of zero

                    heatflux_time_stack.append(date_LST.hour - sunrise)#appending the time of the first hour after sunrise

                    heat_indicator = 1 #indicates that sunrise has already passed and there has been no missing data since sunrise

                    A = 0.2

                    RHS_integral = scipy.integrate.simps(heatflux_sum_stack,heatflux_time_stack)

                    RHS = ((1 + 2 * A) / (rho * cp)) * RHS_integral

                    theta_0 = ref_temp[i] * (1000 / pres [i])**0.286

                    z_ic, = scipy.optimize.fsolve(vfunc, z_guess)

                elif heat_indicator == 1 and Phi_S > 0 : #heatflux after sunrise

                    heatflux_sum_stack.append(heatflux_sum_stack[-1] + Phi_S)#appending the heatflux for the current instance

                    heatflux_time_stack.append(date_LST.hour - sunrise) #appending the hour after sunrise for the current instance

                    A = 0.2#empirical constant (Deardorff, 1980)

                    RHS_integral = scipy.integrate.simps(heatflux_sum_stack,heatflux_time_stack) #calculating the integral of the heatflux w.r.t. time from sunrise to the time of the current instance (hours)

                    RHS = ((1 + 2 * A) / (rho * cp)) * RHS_integral

                    theta_0 = ref_temp[i] * (1000 / pres [i])**0.286

                    z_ic, = scipy.optimize.fsolve(vfunc, z_guess)

                else:

                    stoppage_criteria = len(heatflux_sum_stack) 

                    while stoppage_criteria > 0:

                        heatflux_sum_stack.pop(0)

                        heatflux_time_stack.pop(0)

                        stoppage_criteria = len(heatflux_sum_stack) 
                    
                    heat_indicator = 0
                    z_ic = float('nan')

            elif zic_indicator == 2:#calculating convective mixing layer height

                z_ic = lcl(pres[i]*100,ref_temp[i],rel_hum[i]) #convective mixing layer height

            else: #error in z_ic_indicator

                print("erroneous value given for z_ic_indicator")

                sys.exit()
                
            z_i = z_ic #setting mixing layer height as the convective mixing layer height

            w = ((g_0 * Phi_S * z_i) / (rho * cp * ref_temp[i]))**(1/3) #Deardorff velocity


        else: #calculating weather parameters for the SBL phase

            PTS = 0.09 * (1 - 0.5 * tot_cloud[i]**2) #potential temperature using fractional total cloud cover (Van Ulden and Holtslag, 1985):

            CD = kv / np.log(z_uref / z_0) #drag coefficient

            u_0 = np.sqrt((5.0 * z_uref * g_0 * PTS) / ref_temp[i])

            u_crit = 4.0 * u_0 / CD #critical wind

            SFV =(CD * wind_spd[i] / 2.0) * (-1 + (1 + ((2 * u_0) / (wind_spd[i] * CD**0.5))**2)**0.5)

            L = (ref_temp[i] * SFV**2) / (kv * g_0 * PTS)

            if abs(SFV / (f * L)) <= 4:

                z_ie = 0.2 * SFV / f

            else:
                z_ie = 0.4 * abs((SFV * L) / f)

            #z_ie2 = -f * L + np.sqrt((f**2) * (L**2) + 0.57 * SFV * f * L)            
    
            if time_indicator == 0: #mechinal mixing height as the equivalent mixing height

                z_im = z_ie#mechanical mixing layer height set to equivalent mixing height

            elif np.isnan(z_im_stack[-1]):#no estimate of the equivalent mixing height for the previous instance

               z_im = z_ie #mechanical mixing layer height set to equivalent mixing height

            elif heatflux_stack[-1] >= 0 :#previous instance was in the CBL

                z_im = z_ie #mechanical mixing layer height set to equivalent mixing height
                
            else: #mechanical mixing height by time evolution equation

                tau = z_im_stack[-1] / (2 * SFV) #time scale

                z_im = z_im_stack[-1] * np.exp(-3600 / tau) + z_ie * (1 - np.exp(-3600 / tau)) #mechanical mixing height by time evolution (Venkatram 1982)

            z_ic = float('nan') #convective mixing layer does not exist in the SBL

            z_i = z_im #mixing layer height set to mechanical mixing layer height

            Phi_S = -rho * cp * SFV * PTS #sensible surface heat flux

            w = float('nan') #no Deardorff velocity (convective scale) in the SBL
            
            heat_indicator = 0

        z_im_stack.append(z_im) #adding mechanical mixing layer height of current instance to stack

        z_ic_stack.append(z_ic) #adding convective mixing layer height of current instance to stack

        zi_stack.append(z_i) #adding mixing layer height of current instance to stack

        heatflux_stack.append(Phi_S) #adding sensible surface heatflux of current instance to stack

        L_stack.append(L) #adding Obukhov length of current instance to stack

        sfv_stack.append(SFV) #adding shear velocity of current instance to stack

        pts_stack.append(PTS) #adding potential temperature scale of current instance to stack

        w_stack.append(w) #adding Deardorff of current instance to stack

        date_stack.append(date_LST) #adding date of current instance to stack


file = open(outfile,"w")

for i in range(len(date_stack)):
    file.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(date_stack[i],round(float(wind_spd[i]),3),float(wind_dir[i]),float(pres[i]),float(ref_temp[i]),float(rel_hum[i]),round(float(heatflux_stack[i]),6),round(float(L_stack[i]),6),round(float(sfv_stack[i]),6),round(float(pts_stack[i]),6),round(float(w_stack[i]),6),round(float(zi_stack[i]),6),round(float(z_ic_stack[i]),6),round(float(z_im_stack[i]),6)))

file.close()

file = open(outfile2,"w")

for i in range(len(date_stack)):

    if heatflux_stack[i] >= 0: #indicates cbl values

        file.write("{},{},{},{},{}\n".format(date_stack[i],round(float(B0_stack[i]),6),round(float(SC_stack[i]),6),round(float(a_stack[i]),6),round(float(heatflux_stack[i]),6)))
  
