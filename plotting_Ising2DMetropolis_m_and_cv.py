from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# Available data files
# beta0to10_Nbetas100_J1_mcsteps1000_bins100_IsingMC.txt
# beta0to0p01_Nbetas100_J1_mcsteps1000_bins100_IsingMC.txt
# beta0to0p1_Nbetas100_J1_mcsteps1000_bins100_IsingMC.txt

# Open the file for reading
infile = open("beta0to10_Nbetas100_J1_mcsteps1000_bins100_IsingMC.txt", "r")


# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)   # Tested, and this is done successfully
L = sqrt(N);

# Getting lists ready to store the data
betas = []      # List of beta values
m_avs = []      # Average m over all bins and mcsteps, listed according to beta value
blockvars = []  # Variance found by traversing over the bins
cvs = []        # Heat capacity
eavs = []       # Average energy of the total system
esqavs = []     # Average energy**2 of the total system

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # Betas
        beta = float(words[0])
        betas.append(beta)
        # m_avs
        m_av = float(words[1])
        m_avs.append(m_av)
        # blockvars
        blockvar = float(words[2])
        blockvars.append(blockvar)
        # cvs
        cv = float(words[3])
        cvs.append(cv)
        # eavs
        eav = float(words[4])
        eavs.append(eav)
        # esqavs
        esqav = float(words[5])
        esqavs.append(esqav)       
        
# We prefer arrays
betas = array(betas)
m_avs = array(m_avs)
blockvars = array(blockvars)
cvs = array(cvs)
eavs = array(eavs)
esqavs = array(esqavs)

# Remember to close the file
infile.close()

# Doing the plotting thing
figure()
plot(betas, m_avs, 'r')
title('Average squared magnetization $<m>^2$ in the Ising model, spin=', s, ' L=', L)
xlabel(r'$\beta$')
ylabel(r'$<m>^2$')
show()

figure()
errorbar(betas, m_avs, yerr=blockvars)
title('$<m>^2$ in the Ising model with error bars, s=', s, ' L=', L)
xlabel(r'$\beta$')
ylabel(r'$<m>^2$')
show()

figure()
plot(betas, blockvars, 'r')
title(r'Variance $\sigma_B$ of $<m>^2$ in the Ising model, s=', s, 'L=', L)
xlabel(r'$\beta$')
ylabel(r'$\sigma_B$')
show()

figure()
plot(betas, cvs, 'r')
title(r'Heat capacity $C_v$ in the Ising model, s=,'s, 'L=', L)
xlabel(r'$\beta$')
ylabel(r'$C_v$')
show()

figure()
plot(betas, eavs, 'r')
title(r'Energy $<E>$ in the Ising model, s=,'s, 'L=', L)
xlabel(r'$\beta$')
ylabel(r'$<E>$')
show()

figure()
plot(betas, esqavs, 'r')
title(r'$<E>^2$ in the Ising model, s=,'s, 'L=', L)
xlabel(r'$\beta$')
ylabel(r'$<E>$')
show()


