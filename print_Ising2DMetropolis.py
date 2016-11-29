from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# Open the file for reading
infile = open("testrun3_IsingMC.txt", "r")

# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)   # Tested, and this is done successfully

# Getting lists ready to store the data
betas = []
binhold_m = []
m_beta_average = []
all_ms = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

counter = 0
j = 0

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        if counter==(bins-1):
            beta = words[0]
            beta = float(beta)
            betas.append(beta)
        m = words[1]
        m = float(m)
        binhold_m.append(m)
        counter +=1
        # When we are done with each temperature:
        if(counter==bins):                              # After getting a result from each bin                          
            m_beta_average.append(0)                    # Creating a new element in the list
            for i in range(0,bins):                     # Loop over all the bins
                m_beta_average[j] += binhold_m[i]/bins  # Indexing said element
                # Do some error analysis here
            j += 1                                      # Updating index
            binhold_m = []                              # Reset binhold_m for the next temperature
            counter = 0                                 # Reset the counter so that we can do the same for the next temperature
        
        
# We prefer arrays
betas = array(betas)
m_beta_average = array(m_beta_average)
m_beta_average_divs = m_beta_average/(s*s)
mtotal_beta_average = N*N*m_beta_average
RS_mtotal_beta_average = sqrt(mtotal_beta_average)
RMS_m_beta_average = sqrt(m_beta_average)
RMS_m_beta_average_divs = sqrt(m_beta_average_divs)

# Remember to close the file
infile.close()

# Doing the plotting thing
figure()
plot(betas, m_beta_average, 'r')
title('Average squared magnetization vs temperature in the Ising model')
xlabel(r'$\beta$')
ylabel('m')
show()

figure()
plot(betas, m_beta_average_divs, 'r')
title('Average squared magnetization vs temperature in the Ising model, with s=1')
xlabel(r'$\beta$')
ylabel('m')
show()

figure()
plot(betas, mtotal_beta_average, 'r')
title('Total squared magnetization vs temperature in the Ising model, with s=1')
xlabel(r'$\beta$')
ylabel('m')
show()

figure()
plot(betas, RMS_m_beta_average, 'r')
title('RMS magnetization vs temperature in the Ising model')
xlabel(r'$\beta$')
ylabel('m')
show()

figure()
plot(betas, RMS_m_beta_average_divs, 'r')
title('RMS magnetization vs temperature in the Ising model, with s=1')
xlabel(r'$\beta$')
ylabel('m')
show()

figure()
plot(betas, RS_mtotal_beta_average, 'r')
title('Root squared magnetization vs temperature in the Ising model, with s=1')
xlabel(r'$\beta$')
ylabel('m')
show()








