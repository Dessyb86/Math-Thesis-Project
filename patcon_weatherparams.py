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
#                   entitled "heatflux_comparison_i_j_k_l" with values of the bowen ratio, attenuation  #
#                   factor, clear sky global solar radiation, net radiation, sensible surface           #
#                   heatflux and convective mixing height will be printed                               #
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

############## CBL Parametecrs (If CBL calculations are involved)###############

#The following conditions define the formula with which the heatflux for the CBL should be calculated

B_10 = 0.1310428 #constant coefficient in calculating atmospheric transmittance (Hotel, 1976)

B_11 = 0.7342016 #coefficient of the exponential in calculating atmospheric transmittance (Hotel, 1976)

k_t = 0.1368551 #transmittance coefficient in calculating atmospheric transmittance (Hotel, 1976)

#The following must be specified if B_0_indicator  was set to 1

rainy_season = [5,6,7,8,9,10,11] #array of the number of the months in the rainy season

bowen_rainy = 0.6 #bowen ratio value for rainy season

bowen_dry = 2.0 #bowen ratio value for dry season

#The following must be specified if A_indicator was set to either 2 or 3

cloud1_col = 17 #number of the columns with cloud information

cloud2_col = 18 #number of the columns with cloud information

cloud3_col = 19 #number of the columns with cloud information

cloud4_col = 20 #number of the columns with cloud information

##################### DO NOT TAMPER BEYOND THIS POINT #####################

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

#################################### Defining Global Constants ############################

date1 = datetime.strptime(str(2000)+'-'+str(9)+'-'+str(17)+':'+str(11),'%Y-%m-%d:%H') #arbitrary datetime variable (2000-11-17 11:00:00)
date2 = datetime.strptime(str(2000)+'-'+str(9)+'-'+str(17)+':'+str(12),'%Y-%m-%d:%H') #arbitrary datetime variable (2000-11-17 12:00:00)

time_check = date2 - date1 #stores the value of one hour to be used in comparisons

null_time = date2 - date2

phi_lat = lat * np.pi / 180 #latitudinal angle of the location

omega = 0.0000729 #angular velocity of the earth
                
f = 2 * omega * np.sin(phi_lat) #coriolis parameter for the tropics (Blackadar and Tennekes, 1968) 

cpa = 1003.5 #specific heat capacity of dry air at constant pressure

cva = 1930 #specific heat capacity of water vapour at constant pressure

Ra = 287.05 #gas constant for dry air

Rv = 461.52 #gas constant for water vapour

kv = 0.4 #vonKarman constant

g_0 = 9.80665 #standard gravity

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

############################## EXTRACTING REQUISITE COLUMNS FROM INPUT DATA FILE #############################

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

#################################### CONVERTING UNITS OF INPUT DATA ############################
                
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

################################## PROCESSING INPUT DATA TO PRODUCE REQUISITE OUTPUT DATA ##################################

for zic_ind in range(3):
    
    for B_0_ind in range(2):

        for R_SC_ind in range(2):

            for A_ind in range(3):

                B_0_indicator = B_0_ind + 1 

                #B_0_indicator = 1 - Bowen ratio will be calculated by aermod formulation (using tables of mean seasonal values)

                #B_0_indicator = 2 - Bowen ratio will be calculated from reference temperature and relative humidity (Lin et al., 2016)

                R_SC_indicator = R_SC_ind + 1 # indicates what formula should be used to calculate the sky clear global radiation. 

                #R_SC_indicator = 1 - sky clear radiation will be calculated by aermod formulation (Kasten and Czeplak, 1980)

                #R_SC_indicator = 2 - sky clear global radiation will be calculated from julian day and zenith angle (Khan and Ahmad, 2012)

                A_indicator = A_ind + 1 # indicates what formula should be used to calculate the attenuation factor. 

                #A_indicator = 1 - attenuation factor will be calculated by aermod formulation (Kasten and Czeplak, 1980)

                #A_indicator = 2 - attenuation factor will be calculated by weight formula using cloud cover and transmittance for each layer from Katen and Czeplak, 1980

                #A_indicator = 3 - attenuation factor will be calculated by weight formula using cloud cover and transmittance for each layer from Walled and Yussra, 2006

                zic_indicator =  zic_ind + 1 #indicates what formula will be used to calculate the convective mixing layer height in the CBL

                #zic_indicator = 1 - convective mixing height will be calculated by the temporal evolution (Batchvarova and Gryning (1991))

                #zic_indicator = 2 - convective mixing height will be calculated as the exact expression of the lifting condesation level (Romps, 2017)

                #zic_indicator = 3 - convective mixing height will be calculated by the temporal evolution (Thomson (1992), Thomson (2000))


                #################################### Output Files ###############################

                outfile = "weatherparams_"+str(int(zic_indicator))+"_"+str(int(B_0_indicator))+"_"+str(int(R_SC_indicator))+"_"+str(A_indicator)

                outfile2 = "heatflux_comparison_"+str(int(zic_indicator))+"_"+str(int(B_0_indicator))+"_"+str(int(R_SC_indicator))+"_"+str(A_indicator)

                ##################################### Setting Individual Cloud Layer Data #############################################
                
                if A_indicator == 2 or A_indicator == 3:
                    
                    cloud1 = csv1[:,cloud1_col] #first cloud group

                    cloud2 = csv1[:,cloud2_col] #second cloud group

                    cloud3 = csv1[:,cloud3_col] #third cloud group

                    cloud4 = csv1[:,cloud4_col] #fourth cloud group

                ################################## Checking Rainy Season Array #########################################

                if B_0_indicator == 1:
                    
                    if len(rainy_season) == 0:

                        print("rainy season array is empty")

                        sys.exit()

                    for i in range(len(rainy_season)):

                        if rainy_season[i] < 1 or rainy_season[i] > 12 or (rainy_season[i] % 1 != 0):

                            print("month values in rainy_season array must be whole numbers between 1 and 12 inclusive")

                            sys.exit()

                ########################## Defining Storage Variables #############################

                a_stack =[] #stack of cloud transmittance values for special print out

                B0_stack =[] #stack of bowen ratio values for special print out

                SC_stack =[] #stack of clear sky global radiation values for special print out

                z_im_stack =[] #stack of mechanical mixing height values

                z_ic_stack = [] #stack of convective mixing height values

                z_ic2_stack = [] #stack of convective mixing height values for CBL

                zi_stack = [] #stack of mixing layer height values

                heatflux_stack =[] #stack of heatflux values

                L_stack =[] #stack of obukhov length values

                sfv_stack =[] #stack of shear velocity values

                pts_stack =[] #stack of potential temperature scale values

                w_stack = [] #stack of Deardorff values

                date_stack = [] #stack of date values

                date2_stack = [] #stack of date values for special print out

                R_N_stack = [] #stack of net radiation values for special print out

                R_S_stack = [] #stack of sensible surface heat flux values for special print out

                z_stack = [] #stack of convective mixing height values for special print out

                heatflux_sum_stack = [] #stack of the running sum of the sensible surface heatflux in the CBL

                heatflux_time_stack = [] #stack of the time elapsed in the CBL for the particular day

                ############################ Calculation of Weather Parameters #######################################
                
                for i in range(len(year)):

                    es = 610.78 * np.exp(17.2694 * ((ref_temp[i] - 273.16)/(ref_temp[i] - 273.16 + 238.3))) #saturation vapour pressure

                    e = rel_hum[i] * es # vapour pressure 

                    w = (Ra/Rv) * (e / (pres[i] * 100 - e)) #mixing ratio

                    q = w / (w + 1) #specific humidity

                    Rm =  (1 - q) * Ra + q * Rv # gas constant for moist air

                    rho = (pres[i] * 100) / (Rm * ref_temp[i]) #calculating density of the mosit atmosphere using the ideal gas law

                    cp =  (1 - q) * cpa + q * cva # specific heat capacity at constant pressure for the atmosphere

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

                        ########## determining the phase of the PBL (pbl_indocator) ##########

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

                        sunrise = sun['sunrise'].hour + (sun['sunrise'].minute / 60) + (sun['sunrise'].second / 3600) #time of sunrise in hours
                        
                        sunset = sun['sunset'].hour + (sun['sunset'].minute / 60) + (sun['sunset'].second / 3600) #time of sunset in hours
                       
                        if date_LST.hour >= sunrise and date_LST.hour <= sunset: #current instance is within the daytime

                            ############### calculating sensible surface heat flux (Phi_S) #################

                            n_j = float(date_LST.strftime('%j')) #day number of the year for the current instance

                            phi_D = ((-23.45*np.pi)/180.0)* np.cos(2*np.pi*((n_j + 10.0)/365.25)) #declination angle

                            phi_H_current = np.pi*(float(date_LST.hour - sunrise)/float(sunset - sunrise))-((np.pi)/2) #hour angle of the current instance

                            phi_H_previous = np.pi*(float((date_LST.hour - 1) - sunrise)/float(sunset - sunrise))-((np.pi)/2) #hour angle of the hour before the current instance

                            phi_Z_current = float(np.arccos(np.sin(phi_lat) * np.sin(phi_D) + np.cos(phi_lat) * np.cos(phi_D) * np.cos(phi_H_current))) #zenith angle of the current instance

                            if date_LST.hour - 1 < sunrise: #adjusts the zenith angle of the previous hour to the zenith angle at sunrise if the previous hour was before sunrise

                                phi_Z_previous = np.pi / 2

                            else: # previous hour was between sunrise and sunset
                            
                                phi_Z_previous = float(np.arccos(np.sin(phi_lat) * np.sin(phi_D) + np.cos(phi_lat) * np.cos(phi_D) * np.cos(phi_H_previous))) #zenith angle of the hour before the current instance

                            albedo1 = albedo_0 + (1 - albedo_0) * np.exp(((-18 * phi_Z_previous) / np.pi) - 0.5 * (1 - albedo_0)**2) #albedo at current instance

                            albedo2 = albedo_0 + (1 - albedo_0) * np.exp(((-18 * phi_Z_current) / np.pi) - 0.5 * (1 - albedo_0)**2) #albedo oe hour prior to current instance 

                            mean_albedo = (albedo1 + albedo2)/2 #mean albedo                                                                                  
                                                                            
                            if R_SC_indicator == 1: #calculating sky clear global radiation

                                Phi_SC = 990.0 * (np.cos(phi_Z_current)) - 30.0 #sky clear global radiation

                            elif R_SC_indicator == 2: #calculating sky clear global radiation
                                
                                if phi_Z_current == np.pi / 2:#calculating tb when zenith angle is pi/2

                                    tb = B_10 #atmospheric transmittance for direct beam 

                                else: #calculating sky clear radiation

                                    tb = B_10 + B_11 * np.exp(-k_t / np.cos(phi_Z_current)) #atmospheric transmittance for direct beam

                                td = 0.271 - 0.294 * tb #atmospheric transmittance for diffused beam

                                Q_on = 1367.0 * (1.0+(7.0/1367.0)*np.cos((2.0*np.pi*(n_j-3))/365.25)) #value of solar irradiance at the top of the atmosphere

                                Qb = Q_on * tb #direct beam portion of the global solar radiation

                                Qd = Q_on * td * np.cos(phi_Z_current) #diffused beam portion of the global solar radiation
                                
                                Phi_SC = Qb + Qd #sky clear global solar radiation

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
                                    
                            else: #calculating bowen ratio                

                                B_0 = 1.46 * (1.0/rel_hum[i]) * ((ref_temp[i] / 273.0)**2) * (np.exp(-19.83 * (1.0 - (273.0 / ref_temp[i])))) #bowen ratio
                            
                            if A_indicator == 1:#calculating cloud attenuation factor

                                A = 1 - 0.75 * (tot_cloud[i])**3.4

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
                 
                            Phi_SO = Phi_SC * A #global solar radiation

                            Phi_N = ((1 - mean_albedo) * Phi_SO + (5.31 * 10**(-13)) * (ref_temp[i])**6 - (5.67 * 10**(-8)) * ref_temp[i]**4 + 60 * tot_cloud[i]) / (1 + 0.12) #net radiation

                            Phi_S = (0.9 * Phi_N) / (1 + (1 / B_0))

                            if Phi_S > 0:#setting pbl_indicator to indicate CBL phase

                                pbl_indicator = 1 #PBL set to CBL phase

                                hour_check += 1 #determines the nth hour of the CBL phase

                            else: #setting pbl_indicator to indicate the phase is not the CBL
                                
                                pbl_indicator = 0

                                hour_check = 0

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

                            if zic_indicator == 1:

                                if hour_check == 1: #first hour after sunrise

                                    gamma = (0.013**2 * ref_temp[i]) /g_0

                                    if R_SC_indicator == 1: #estimating number of seconds since the onset of the CBL

                                        hour_0 = sunrise + ((sunset - sunrise) / np.pi) * (-np.arccos(((1 / 33) - np.sin(phi_D) * np.sin(phi_lat)) / (np.cos(phi_D) * np.cos(phi_lat))) + np.pi/2)

                                        delta_t = min(3600,3600 * (date_LST.hour - hour_0))

                                    else:

                                        delta_t = 3600 * (date_LST.hour - sunrise)

                                    if SFV == 0 or np.isnan(SFV) : #calm wind conditions

                                        z_ic = ((((2 * 0.2 + 1) * Phi_S) / (rho * cp * gamma)) * delta_t)**0.5

                                    else: #positive wind speeds

                                        z_ic = min(((2 * 5 * ref_temp[i] * delta_t) / (g_0 * gamma))**(1 / 3) * SFV, (5 * SFV**3 * ref_temp[i] * rho * cp) / (g_0 * Phi_S * abs(0.2 - 1)))

                                elif hour_check > 1: #more than first hour after sunrise

                                    gamma = (0.013**2 * ref_temp) /g_0

                                    prev_gamma = (0.013**2 * prev_ref_temp) /g_0

                                    delta_t = 3600

                                    if np.isnan(prev_SFV) or prev_SFV == 0: #calm wind conditions

                                        z_ic = ((((2 * 0.2 + 1) * Phi_S) / (rho * cp * gamma)) * delta_t)**0.5

                                    else: #positive wind speeds
                                        
                                        S = ((0.2 * prev_Phi_S) / (prev_rho * prev_cp)) + ((5 * prev_SFV**3 * prev_ref_temp) / (g_0 * z_ic))
                                    
                                        z_ic = (((2*S + (prev_Phi_S / (prev_rho * prev_cp))) / (prev_gamma)) * delta_t)**0.5

                                else: #emptying the heatflux_sum_stack variable

                                    z_ic = float('nan')

                            elif zic_indicator == 2:#calculating convective mixing layer height

                                z_ic = lcl(pres[i]*100,ref_temp[i],rel_hum[i]) #convective mixing layer height

                            elif zic_indicator == 3:

                                if hour_check == 1: #first hour after sunrise

                                    gamma = (0.013**2 * ref_temp[i]) /g_0

                                    if R_SC_indicator == 1:

                                        hour_0 = sunrise + ((sunset - sunrise) / np.pi) * (-np.arccos(((1 / 33) - np.sin(phi_D) * np.sin(phi_lat)) / (np.cos(phi_D) * np.cos(phi_lat))) + np.pi/2)

                                        delta_t = min(3600,3600 * (date.LST.hour - hour_0))

                                    else:

                                        delta_t = 3600 * (date.LST - sunrise)

                                    if np.isnan(SFV) or SFV == 0: #calm wind conditions

                                        z_ic = ((((2 * 0.2 + 1) * Phi_S) / (rho * cp * gamma)) * delta_t)**0.5

                                    else: #positive wind speeds

                                        z_ic = min(((2 * 5 * ref_temp[i] * delta_t) / (g_0 * gamma))**(1 / 3) * SFV, (5 * SFV**3 * ref_temp[i] * rho * cp) / (g_0 * Phi_S * abs(0.2 - 1)))

                                elif hour_check > 1: #heatflux after sunrise

                                    gamma = (0.013**2 * ref_temp) /g_0

                                    prev_gamma = (0.013**2 * prev_ref_temp) /g_0

                                    delta_t = 3600

                                    if np.isnan(prev_SFV) or prev_SFV == 0:

                                        z_ic = ((((2 * 0.2 + 1) * Phi_S) / (rho * cp * gamma)) * delta_t)**0.5
                                    
                                    else:

                                        z_ic = z_ic + ((1.4) * (prev_Phi_S / (prev_rho * prev_cp * prev_gamma * z_ic)) + ((prev_SFV**3 * prev_ref_temp) / (g_0 * prev_gamma * z_ic**2))) * delta_t

                                else: #emptying the heatflux_sum_stack variable

                                    z_ic = float('nan')
                                    
                            else: #error in z_ic_indicator

                                print("erroneous value given for z_ic_indicator")

                                sys.exit()
                                
                            z_i = z_ic #setting mixing layer height as the convective mixing layer height

                            w = ((g_0 * Phi_S * z_i) / (rho * cp * ref_temp[i]))**(1/3) #Deardorff velocity

                            ######################## Appending Data to be Printed ##############################

                            date2_stack.append(date_LST)
                            a_stack.append(A)
                            SC_stack.append(Phi_SC)
                            B0_stack.append(B_0)
                            R_N_stack.append(Phi_N)
                            R_S_stack.append(Phi_S)
                            z_stack.append(z_ic)

                        else:
                            
                            A = float('nan')
                            B_0 = float('nan')
                            Phi_N = float('nan')
                            Phi_S = float('nan')
                            z_ic = float('nan')
                            L,SFV,PTS = patcon.sbl_scales_drag(wind_spd[i],ref_temp[i],tot_cloud[i]) #calculating obukhov length, surface frictional velocity and potential temperature scale
                            Phi_S = -PTS * rho * cp * SFV
                            w = float('nan')

                            hour_check = 0

                        prev_Phi_S = Phi_S

                        prev_ref_temp = ref_temp[i]

                        prev_rho = rho

                        prev_SFV = SFV

                        prev_cp = cp

                file = open(outfile2,"w")

                for i in range(len(date2_stack)):

                    if R_S_stack[i] > 0: #indicates cbl values

                        file.write("{},{},{},{},{},{},{}\n".format(date2_stack[i],np.round(float(B0_stack[i]),6),np.round(float(SC_stack[i]),6),np.round(float(a_stack[i]),6),np.round(R_N_stack[i],6),np.round(R_S_stack[i],6),np.round(z_stack[i],6)))
