from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# Open the file for reading
#infile = open("beta0to5_Nbetas100__spin0p5_L15_J1_mcsteps1000_bins100_IsingMC.txt", "r")
infile = open("beta0to5_Nbetas100_spin0p5_L15_J1_mcsteps1000_bins100_cverrorattempt_IsingMC.txt", "r")


# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)   # Tested, and this is done successfully
L = sqrt(N);

# Getting lists ready to store the data
betas = []                  # List of beta values
m_avs = []                  # Average m over all bins and mcsteps, listed according to beta value
blockvars = []              # Variance found by traversing over the bins
cvs = []                    # Heat capacity
eavs = []                   # Average energy of the total system
esqavs = []                 # Average energy**2 of the total system
eavs_std = []               # Average energy of the total system
esqavs_std = []             # Average energy**2 of the total system
cv_avfromblocks_std = []    # Attempt to get the error in cv
cv_avfromblocks     = []    # To compare the value I get from this and the value I get when finding Cv the usual way

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
        # eavs_std
        eav_std = float(words[6])
        eavs_std.append(eav_std)
        # esqavs_std
        esqav_std = float(words[7])
        esqavs_std.append(esqav_std)
        # cv_avfromblocks_std
        cvblock_std = float(words[8])
        cv_avfromblocks_std.append(cvblock_std)
        # cv_avfromblocks
        cvblock = float(words[9])
        cv_avfromblocks.append(cvblock)  
        
# We prefer arrays
betas = array(betas)
m_avs = array(m_avs)
blockvars = array(blockvars)
cvs = array(cvs)
eavs = array(eavs)
esqavs = array(esqavs)
eavs_std = array(eavs_std)
esqavs_std = array(esqavs_std)
cv_avfromblocks_std = array(cv_avfromblocks_std)
cv_avfromblocks     = array(cv_avfromblocks)

# Remember to close the file
infile.close()

# Doing the plotting thing
figure()
plot(betas, m_avs, 'r')
title('Average squared magnetization $<m>^2$ in the Ising model, L =%s'%L)
xlabel(r'$\beta$')
ylabel(r'$<m>^2$')
show()

figure()
errorbar(betas, m_avs, yerr=blockvars)
title('Average squared magnetization $<m>^2$ in the Ising model, L =%s'%L)
xlabel(r'$\beta$')
ylabel(r'$<m>^2$')
show()

figure()
plot(betas, blockvars, 'r')
title(r'Standard deviation $\sigma_B$ of $<m>^2$ in the Ising model, L =%s'%L)
xlabel(r'$\beta$')
ylabel(r'$\sigma_B$')
show()

figure()
plot(betas, cvs, 'r')
title(r'Heat capacity $C_v$ in the Ising model, L =%s'%L)
xlabel(r'$\beta$')
ylabel(r'$C_v$')
show()

figure()
errorbar(betas, cvs, yerr=cv_avfromblocks_std)
title(r'Heat capacity $C_v$ in the Ising model with errobars, L =%s'%L)
xlabel(r'$\beta$')
ylabel(r'$C_v$')
show()

figure()
plot(betas, cvs, label='Usual')
hold('on')
plot(betas, cv_avfromblocks, label='Av from blocks')
title(r'Heat capacity $C_v$ in the Ising model, L =%s'%L)
xlabel(r'$\beta$')
ylabel(r'$C_v$')
legend(loc='upper right')
show()

figure()
errorbar(betas, eavs, yerr=eavs_std)
title(r'Energy $<E>$ in the Ising model, L =%s'%L)
xlabel(r'$\beta$')
ylabel(r'$<E>$')
show()

figure()
errorbar(betas, esqavs, yerr=esqavs_std)
title(r'$<E>^2$ in the Ising model, L =%s'%L)
xlabel(r'$\beta$')
ylabel(r'$<E>$')
show()


