from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# Should divide into temperatures! Now we loose the information of autocorrelation time vs temp.
# Function to find the autocorrelation time

def autocorrelation_time(j, mcsteps, ms, plotbool='false'): # Want the option to plot
    lowerlimiti = j*mcsteps
    upperlimiti = (j+1)*mcsteps
    upperlimitdt = int(floor(0.1*mcsteps))    
    dts = range(0,upperlimitdt)
    A = zeros(len(dts))
    act = 0
    for i in range(lowerlimiti,upperlimiti):     # Choosings the spins in our bin  
        for dt in dts:                           # A(dt) for different dt 
            a = 0 
            counter = 0
            while (i+dt)<upperlimiti:     # Cutting off when we reach the end of the bin           
                A[dt] += ms[i]*ms[i+dt]   # Correlation between the spins in our bin
                counter +=1               # Number of terms
            A[dt] = A[dt]/counter         # Or drop the counter and divide by mcsteps-delta+1
            act +=A[dt]
    
    if plotbool=='true':                  # Just in case we want to plot.
        figure()
        plot(dts, A)
        title('Unintegrated autocorrelation time')
        xlabel(r'$\Delta\tau$')
        ylabel('A')
        show()
        
    return act
        

print "Opening file."
## Reading in the data ##
# Open the file for reading
infile = open("beta0to5_Nbetas100__spin0p5_L15_J1_mcsteps1000_bins100_IsingMC_all.txt", "r")

readbool = 'false'
if readbool == 'true':
    # The first line contains information about the system. We read that separately
    firstline = infile.readline() # Reads the first line only
    # Number of spins, spin size, number of bins
    N, s, bins, Nbeta, mcsteps, betamin, betamax  = firstline.split()
    N = int(N); s = float(s); bins = int(bins); Nbeta = int(Nbeta); 
    mcsteps = int(mcsteps); betamin = float(betamin); betamax = float(betamax)
    
    L = sqrt(N);
else:
    N = 225; L = 15; s = 0.5; bins = 100; Nbeta = 100; mcsteps = 1000;
    betamin = 0; betamax = 5;


ms = []
es = []

print "Lists initialized. Now splitting into lines."
# Read the lines
lines = infile.readlines()

print "Line splitting done. Reading values into list line for line."
time1 = time.time()
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
      
time2 = time.time()
timediff = time2 - time1
print "Done with reading values into list. Time taken to retrieve data from file: ", timediff
print "Now converting lists to arrays."        
# We prefer arrays
ms = array(ms)
es = array(es)

print "Lists successfully converted into arrays. Now closing the file."

# Remember to close the file
infile.close()

## Finding the autocorrelation function ##
# For ms
#endpoint = floor(0.1*len(ms))
endpoint = 100
betas = linspace(betamin, betamax, Nbeta)      # Beta array
macts = zeros(Nbeta)                           # Autocorrelation time array 
macts_std = zeros(Nbeta)                       
 

print "File closed and arrays for autocorrelation plots made. Now feeding in the correlation."

counter = 0
for i in range(Nbeta):
    # For each beta
    mact_av = 0                  # To find the average value
    macts_inbin = zeros(bins)    # For storing the standard deviation
    for j in range(bins):
        # Run over every bin
        act = autocorrelation_time(counter, mcsteps, ms)     # Finding the integrated autocorrelation time
        mact_av += act
        macts_inbin[j] = act   
    # Feeding the average value in
    mact_av = mact_av/bins
    macts[i] = mact_av
    # Calculating the standard deviation  
    for k in range(bins):
        macts_std[i] += (macts_inbin[k]-mact_av)*(macts_inbin[k]-mact_av)
    macts_std[i] = macts_std/(bins*(bins-1))
    
    

print "Done with finding the integrated autocorrelation times. Now plotting."
# Doing the plotting thing
# Not integrated
figure()
errorbar(betas, macts, yerr=macts_std)
title(r'Integrated autocorrelation time $\tau_{int}$ between magnetizations at temporal distance $\Delta\tau$')
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
