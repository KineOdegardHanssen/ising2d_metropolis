from numpy import *
from scipy.misc import *
from matplotlib.pyplot import *
import time
import sys

# Function to find the autocorrelation time ...
def corr(deltat, ms):
    lenms = len(ms)
    a = 0
    dt = int(deltat)
    for i in range(lenms):
        if (i+dt)>=lenms:
            # Are we cyclic in time? In that case, break.
            return a
            # Or periodic boundary conditions? Seems weird, though...
            #index = i+dt-lenms
            #a += ms[i]*ms[index]
        else:
            # Keep accumulating the sum
            a += ms[i]*ms[i+dt]
    return a

def autocorrelation_time(i, Nbeta, ms):
    lenms = len(ms)
    a = 0
    #upperlimitdt = int(floor(0.25*Nbeta)) + 1  # Experiment with this
    upperlimit = (i+1)*Nbeta 
    for j in range(i*Nbeta, upperlimit):
        if (j+dt)>=upperlimit:
            # Are we cyclic in time? In that case, break.
            return a
            # Or periodic boundary conditions? Seems weird, though...
            #index = i+dt-lenms
            #a += ms[i]*ms[index]
        else:
            # Keep accumulating the sum
            a += ms[j]*ms[j+dt]
    return a

## Reading in the data ##
# Set this manually to suit the data file
Nbeta = 100
# Open the file for reading
infile = open("beta0to5_Nbetas100__spin0p5_L15_J1_mcsteps1000_bins100_IsingMC_all.txt", "r")

ms = []
es = []

# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # ms
        m = float(words[0])
        ms.append(m)
        # es
        e = float(words[1])
        es.append(e)
      
        
# We prefer arrays
ms = array(ms)
es = array(es)

# Remember to close the file
infile.close()

## Finding the autocorrelation function ##
# For ms
lenms = len(ms)
#endpoint = int(floor(0.01*lenms))
endpoint = 100
deltats = linspace(0, endpoint, endpoint+1)     # Array of distance between spins, to be plotted against correlations
ams = zeros(endpoint+1)                         # Correlation array


for i in range(len(deltats)):
    print deltats[i]
    ams[i]= corr(deltats[i], ms)

# Sort after the value of beta
j = 0
betasort = linspace(0,5,Nbeta)
actime = zeros(Nbeta)
for i in range(lenms):
    # Loop over dt as well?
    actime[i] += autocorrelation_time(i, Nbeta, ms)
    

# Doing the plotting thing
# Not integrated
figure()
plot(deltats, ams, 'r')
title(r'Autocorrelation time $A(\Delta\tau)$ between magnetizations at temporal distance $\Delta\tau$')
xlabel(r'$\Delta\tau$')
ylabel(r'$A(\Delta\tau)$')
show()

"""
figure()
errorbar(betas, m_avs, yerr=blockvars)
title('$<m>^2$ in the Ising model with error bars')
xlabel(r'$\beta$')
ylabel(r'$<m>^2$')
show()

figure()
plot(betas, blockvars, 'r')
title(r'Variance $\sigma_B$ of $<m>^2$ in the Ising model')
xlabel(r'$\beta$')
ylabel(r'$\sigma_B$')
show()

figure()
plot(betas, cvs, 'r')
title(r'Heat capacity $C_v$ in the Ising model')
xlabel(r'$\beta$')
ylabel(r'$C_v$')
show()

figure()
plot(betas, eavs, 'r')
title(r'Energy $<E>$ in the Ising model')
xlabel(r'$\beta$')
ylabel(r'$<E>$')
show()

figure()
plot(betas, esqavs, 'r')
title(r'$<E>^2$ in the Ising model')
xlabel(r'$\beta$')
ylabel(r'$<E>$')
show()
"""

