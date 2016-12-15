from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

def autocorrelation_time(j, mcsteps, bins, ms, plotbool='false', comparebool='false'): # Want the option to plot
    tmax = mcsteps*bins
    lowerlimiti = j*tmax
    upperlimiti = (j+1)*tmax
    upperlimitdt = 100 # Just like it was in the book #Or: int(floor(0.01*mcsteps*N))# Do I really want this large a limit?  
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


s = 0.5; bins = 100; Nbeta = 100; mcsteps = 1000;    

## Reading in the data ##
# Open the file for reading
infile = open("beta1p762747174_Nbetas100_spin0p5_L40_J1_mcsteps1000_bins100_cverrorattempt_IsingMC_all.txt", "r")

ms40 = []
#es15 = []
# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # ms
        m = float(words[0])
        ms40.append(m)
        # es
        #e = float(words[1])
        #es15.append(e)

# Remember to close the file      
infile.close()

# Open the file for reading
infile = open("beta1p762747174_Nbetas100_spin0p5_L30_J1_mcsteps1000_bins100_cverrorattempt_IsingMC_all.txt", "r")

ms30 = []
#es15 = []
# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # ms
        m = float(words[0])
        ms30.append(m)
        # es
        #e = float(words[1])
        #es15.append(e)

# Remember to close the file      
infile.close()

# Open the file for reading
infile = open("beta1p762747174_Nbetas100_spin0p5_L20_J1_mcsteps1000_bins100_cverrorattempt_IsingMC_all.txt", "r")

ms20 = []
#es15 = []
# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # ms
        m = float(words[0])
        ms20.append(m)
        # es
        #e = float(words[1])
        #es15.append(e)

# Remember to close the file      
infile.close()
# Open the file for reading
infile = open("beta1p762747174_Nbetas100_spin0p5_L15_J1_mcsteps1000_bins100_cverrorattempt_IsingMC_all.txt", "r")

ms15 = []
#es15 = []
# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # ms
        m = float(words[0])
        ms15.append(m)
        # es
        #e = float(words[1])
        #es15.append(e)

# Remember to close the file      
infile.close()

# Open the file for reading
infile = open("beta1p762747174_Nbetas100_spin0p5_L10_J1_mcsteps1000_bins100_cverrorattempt_IsingMC_all.txt", "r")

ms10 = []
#es15 = []
# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # ms
        m = float(words[0])
        ms10.append(m)
        # es
        #e = float(words[1])
        #es15.append(e)

# Remember to close the file      
infile.close()

# Open the file for reading
infile = open("beta1p762747174_Nbetas100_spin0p5_L5_J1_mcsteps1000_bins100_cverrorattempt_IsingMC_all.txt", "r")

ms5 = []
#es15 = []
# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # ms
        m = float(words[0])
        ms5.append(m)
        # es
        #e = float(words[1])
        #es15.append(e)

# Remember to close the file      
infile.close()

# Open the file for reading
infile = open("beta1p762747174_Nbetas100_spin0p5_L4_J1_mcsteps1000_bins100_cverrorattempt_IsingMC_all.txt", "r")

ms4 = []
#es15 = []
# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # ms
        m = float(words[0])
        ms4.append(m)
        # es
        #e = float(words[1])
        #es15.append(e)

# Remember to close the file      
infile.close()

# Open the file for reading
infile = open("beta1p762747174_Nbetas100_spin0p5_L3_J1_mcsteps1000_bins100_cverrorattempt_IsingMC_all.txt", "r")

ms3 = []
#es15 = []
# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # ms
        m = float(words[0])
        ms3.append(m)
        # es
        #e = float(words[1])
        #es15.append(e)

# Remember to close the file      
infile.close()

# We prefer arrays
ms40 = array(ms40)
ms30 = array(ms30)
ms20 = array(ms20)
ms15 = array(ms15)
ms10 = array(ms10)
ms5 = array(ms5)
ms4 = array(ms4)
ms3 = array(ms3)

Ls = array([log(3),log(4),log(5),log(10),log(15),log(20),log(30),log(40)])
act = zeros(len(Ls))

## Finding the autocorrelation function ##       
act[0] = log(autocorrelation_time(0, mcsteps, bins, ms3))
act[1] = log(autocorrelation_time(0, mcsteps, bins, ms4))
act[2] = log(autocorrelation_time(0, mcsteps, bins, ms5))
act[3] = log(autocorrelation_time(0, mcsteps, bins, ms10))
act[4] = log(autocorrelation_time(0, mcsteps, bins, ms15))
act[5] = log(autocorrelation_time(0, mcsteps, bins, ms20))
act[6] = log(autocorrelation_time(0, mcsteps, bins, ms30))
act[7] = log(autocorrelation_time(0, mcsteps, bins, ms40))


# Doing the plotting thing
figure()
plot(Ls, act, label=r'$log(\tau_{int})$')
#hold('on')
#plot(Ls, 2*Ls, label='2log(L)')
title(r'Integrated autocorrelation time $\tau_{int}$ at $\beta_{c}$ vs L')
xlabel(r'$log(L)$', fontsize=20)
ylabel(r'$log(\tau_{int})$', fontsize=20)
#legend(loc='lower right')
show()

