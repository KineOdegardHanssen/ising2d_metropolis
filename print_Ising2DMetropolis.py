# Hey ho

from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# Set these as they were when making the file. Change these according to the data set you read


# Opening the file. Change the name as you wish
#infile = open("/home/ubu/Master/Tidlig_Arbeid/2particlePalHuse/pal_huse_2d-build-desktop-Qt_4_8_1_in_PATH__System__Release/randompot_h0to5_J1_av500_res301_filtered_cranked_coordinateFile.txt", "r")

#infile = open("randompot_h0to5_J1_av500_res301_filtered_cranked_coordinateFile.txt", "r")
infile = open("firstrun_IsingMC.txt", "r")

firstline = infile.readline() # Reads the first line (or so I hope)

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)


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
        all_ms.append(m)            # ?
        binhold_m.append(m)
        counter +=1
        # When we are done with each temperature:
        if(counter==bins):
            counter = 0
            m_beta_average.append(0)                  # Creating a new element
            for i in range(0,bins):
                m_beta_average[j] += binhold_m[i]     # Indexing said element
            j += 1                                    # Updating index
            binhold_m = []
        
        
# We prefer arrays
betas = array(betas)
m_beta_average = array(m_beta_average)

# Remember to close the file
infile.close()

# Doing the plotting thing
figure()
plot(betas, m_beta_average, 'r')
title('Average squared magnetization vs temperature in the Ising model') # This is probably too long
xlabel(r'$\beta$')
ylabel('m')
show()



