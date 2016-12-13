from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

def autocorrelation_time(j, mcsteps, bins, ms, plotbool='false'): # Want the option to plot
    tmax = mcsteps*bins
    lowerlimiti = j*tmax
    upperlimiti = (j+1)*tmax
    upperlimitdt = 10000 # Just like it was in the book #Or: int(floor(0.01*mcsteps*N))# Do I really want this large a limit?  
    dts = range(0,upperlimitdt)
    A = zeros(len(dts))
    # Finding A[0] for accumulating act
    term1 = 0
    term2 = 0
    for i in range(lowerlimiti, upperlimiti):
        term1 += ms[i]*ms[i]
        term2 += ms[i]
    A0 = term1/tmax - term2*term2/(tmax*tmax)
    A[0] = 1    # Normalized by A0
    act =  1
    counter = 0
    noofloops = upperlimiti*upperlimitdt
    for dt in range(1, upperlimitdt):                       # A(dt) for different dt 
        term1   = 0
        term2f1 = 0
        term2f2 = 0               
        for i in range(lowerlimiti,upperlimiti-dt):         # Choosing the spins in our bin
            term1 += ms[i]*ms[i+dt]
            term2f1 += ms[i]
            term2f2 += ms[i+dt]
        A[dt] = (term1/(tmax-dt) - term2f1*term2f2/(tmax-dt)**2)/A0
        act += A[dt] 
 
    if plotbool=='true':                  # Just in case we want to plot.
        figure()
        plot(dts, A)
        title('Time-displaced autocorrelation function')
        xlabel(r'$\Delta\tau$')
        ylabel('A')
        show()
        
    return act
        

print "Opening file."
## Reading in the data ##
# Open the file for reading
infile = open("beta0to5_Nbetas100_spin0p5_L15_J1_mcsteps1000_bins100_cverrorattempt_IsingMC_all.txt", "r")

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
betas = linspace(betamin, betamax, Nbeta)      # Beta array
macts = zeros(Nbeta)                           # Autocorrelation time array 
macts_std = zeros(Nbeta)                       
macts2 = zeros(Nbeta)                          # Autocorrelation time array 
macts_std2 = zeros(Nbeta)

"""
act = autocorrelation_time(0, mcsteps, bins, ms, plotbool='true')
print betas[35]
act = autocorrelation_time(35, mcsteps, bins, ms, plotbool='true')
"""           

print "File closed and arrays for autocorrelation plots made. Now feeding in the correlation."

totalnumber = Nbeta*bins
for i in range(Nbeta):                                # We want the integrated autocorrelation for each beta
    print "i = ", i
    act = autocorrelation_time(i, mcsteps, bins, ms)  # Finding the integrated autocorrelation time
    macts[i] = act

print "Done with finding the integrated autocorrelation times. Now plotting."
# Doing the plotting thing
figure()
plot(betas, macts)
title(r'Integrated autocorrelation time $\tau_{int}$ for magnetization')
xlabel(r'$\beta$')
ylabel(r'$\tau_{int}$')
show()

