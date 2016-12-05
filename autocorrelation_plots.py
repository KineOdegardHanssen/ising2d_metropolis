from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# Function to find the autocorrelation time
def corr(deltat, ms):
    a = 0
    for i in range(ms):
        if (i+deltat)>=len(ms):
            # Are we cyclic in time? In that case, break.
            return a
            # Or periodic boundary conditions? Seems weird, though...
            #index = (i+deltat)-len(ms)
            #a += ms(i)*ms(index)
        else:
            # Keep accumulating the sum
            a += ms(i)*ms(i+deltat)
    return a


## Reading in the data ##
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
endpoint = floor(0.25*len(ms))
deltats = linspace(0, endpoint, endpoint+1)     # Array of distance between spins, to be plotted against correlations
ams = zeros(endpoint+1)                         # Correlation array

for i in range(lenms):
    ams[i]= corr(deltats[i], ms)

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
