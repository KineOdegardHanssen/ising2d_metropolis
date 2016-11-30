from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# Open the file for reading
infile = open("beta0to0p1_Nbetas50_mcsteps1000_bins100_IsingMC.txt", "r")

# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)   # Tested, and this is done successfully

# Getting lists ready to store the data
betas = []      # List of beta values
m_avs = []      # Average m over all bins and mcsteps, listed according to beta value
sigma_msqs = [] # sigma_m**2: variance/bins + covariance
variances = []  # Variances
covs = []       # Covariance terms
actimes = []    # Autocorrelation times

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
        # sigma_msq
        sigma_msq = float(words[2])
        sigma_msqs.append(sigma_msq)
        # variances
        var = float(words[3])
        variances.append(var)
        # covs
        cov = float(words[4])
        covs.append(cov)
        # actimes
        act = float(words[5])
        actimes.append(act)
           
        
# We prefer arrays
betas = array(betas)
m_avs = array(m_avs)
sigma_msqs = array(sigma_msqs)
variances = array(variances)
covs = array(covs)
actimes = array(actimes)

# Remember to close the file
infile.close()

# Doing the plotting thing
figure()
plot(betas, m_avs, 'r')
title('Average squared magnetization vs temperature in the Ising model')
xlabel(r'$\beta$')
ylabel(r'$<m>^2$')
show()

figure()
plot(betas, m_avs/(s*s), 'r')  # Normalizing the spin 
title('Average squared magnetization vs temperature in the Ising model, with s=1')
xlabel(r'$\beta$')
ylabel(r'$<m>^2$')
show()

figure()
plot(betas, sigma_msqs, 'r')
title(r'Total variance $\sigma^2_m$ of $<m>^2$ vs temperature in the Ising model')
xlabel(r'$\beta$')
ylabel(r'$\sigma_m$')
show()

figure()
plot(betas, variances, 'r')
title(r'Variance $\sigma^2$ of $<m>^2$ in the Ising model')
xlabel(r'$\beta$')
ylabel(r'$\sigma$')
show()

figure()
plot(betas, covs, 'r')
title(r'Covariance of $<m>^2$ in the Ising model')
xlabel(r'$\beta$')
ylabel('Covariance')
show()

figure()
plot(betas, actimes, 'r')
title(r'Autocorrelation time $\tau$ of $<m>^2$ in the Ising model')
xlabel(r'$\beta$')
ylabel(r'$\tau$')
show()


