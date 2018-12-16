#########################################################################################################
# File Name:        corrective_factors                                                                  #
# Editor:           Patrick Barrett                                                                     #
# Last update:      September 2, 2018                                                                   #
# Descrition :      This script reads the text file containing GIS heights across Barbados and          #
#                   calculates the mean values for the corrective factors                               #
#                                                                                                       #                                   
# input file:       dataset_3.txt                                                                       #
#                                                                                                       #
# output file:      mean values of the corrective factors                                               #
#########################################################################################################

import pandas as pd

data = pd.read_csv('dataset_3.txt', sep=" ", header=None) #reading text file to a pandas dataframe

a_10 = 0.95 #dimensionless constant as suggested by Hottel, 1976

a_11 = 0.98 #dimensionless constant as suggested by Hottel, 1976

a_k = 1.02 #dimensionless constant as suggested by Hottel, 1976

b_10_sum = 0 #sum of corrective factor values for b_10

b_11_sum = 0 #sum of corrective factor values for b_11

k_t_sum = 0 #sum of corrective factor values for k_t

count = 0 #counter

n_row, n_col = data.shape #getting the number of rows and columns in the dataframe

for i in range(n_row): #iterating through the ith row of the data frame

    for j in range(n_col-1): #iterating through the jth column of the dataframe

        if float(data[i][j]) != -9999: #checks whether the data point is within the geography of Barbados or is in the ocean

            print(i,j)

            count += 1 #increments the variable count

            b_10_sum += a_10 * (0.4237 - 0.00821 * (6 - float(data[i][j]))**2) #calculates an instance of b_10 and adds it to the sum

            b_11_sum += a_11 * (0.5055 + 0.00595 * (6.5 - float(data[i][j]))**2) #calculates an instance of b_11 and adds it to the sum

            k_t_sum += a_k * (0.2711 + 0.01858 * (2.5 - float(data[i][j]))**2) #calculates an instance of k_t and adds is to the sum

print(b_10_sum / count) #prints the mean value of B_10

print(b_11_sum / count) #prints the mean value of B_11

print(k_t_sum / count) #prints the mean value of k_t

        
