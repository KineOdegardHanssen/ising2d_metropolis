from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# Open the file for reading
infile = open("beta0to5_Nbetas100_spin0p5_L15_J1_mcsteps1000_bins100_091216_1340_IsingMC.txt", "r")

# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)   # Tested, and this is done successfully
L = sqrt(N);

# Getting lists ready to store the data
betas           = []          # List of beta values
m_avs20         = []          # Average m over all bins and mcsteps, listed according to beta value
stdms20         = []          # Average m over all bins and mcsteps, listed according to beta value
mquadravs20     = []
mquadravs_std20 = []

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
        m_av20 = float(words[1])
        m_avs20.append(m_av20)
        # stdm20
        stdm20 = float(words[2])
        stdms20.append(stdm20)
        # mquadravs20
        mq = float(words[10])
        mquadravs20.append(mq)
        # mquadravs_std20
        mq_std = float(words[11])
        mquadravs_std20.append(mq_std)

        
 
# Remember to close the file
infile.close()

#infile = open("beta0to5_Nbetas100__spin0p5_L15_J1_mcsteps1000_bins100_IsingMC.txt", "r")
infile = open("beta0to5_Nbetas100_spin0p5_L10_J1_mcsteps1000_bins100_091216_1340_IsingMC.txt", "r")

# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)   # Tested, and this is done successfully
L = sqrt(N);

# Getting lists ready to store the data
m_avs10         = []          # Average m over all bins and mcsteps, listed according to beta value
stdms10         = []          # Average m over all bins and mcsteps, listed according to beta value
mquadravs10     = []
mquadravs_std10 = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av = float(words[1])
        m_avs10.append(m_av)
        # stdm20
        stdm = float(words[2])
        stdms10.append(stdm)
        # mquadravs20
        mq = float(words[10])
        mquadravs10.append(mq)
        # mquadravs_std20
        mq_std = float(words[11])
        mquadravs_std10.append(mq_std)  
 
# Remember to close the file
infile.close()

infile = open("beta0to5_Nbetas100_spin0p5_L15_J1_mcsteps1000_bins100_091216_1340_IsingMC.txt", "r")

# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)   # Tested, and this is done successfully
L = sqrt(N);

# Getting lists ready to store the data
m_avs15         = []          # Average m over all bins and mcsteps, listed according to beta value
stdms15         = []          # Average m over all bins and mcsteps, listed according to beta value
mquadravs15     = []
mquadravs_std15 = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av = float(words[1])
        m_avs15.append(m_av)
        # stdm15
        stdm = float(words[2])
        stdms15.append(stdm)
        # mquadravs20
        mq = float(words[10])
        mquadravs15.append(mq)
        # mquadravs_std20
        mq_std = float(words[11])
        mquadravs_std15.append(mq_std) 
        
# Remember to close the file
infile.close()
         
        

infile = open("beta0to5_Nbetas100_spin0p5_L5_J1_mcsteps1000_bins100_091216_1340_IsingMC.txt", "r")

# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)   # Tested, and this is done successfully
L = sqrt(N);

# Getting lists ready to store the data
m_avs5         = []          # Average m over all bins and mcsteps, listed according to beta value
stdms5         = []          # Average m over all bins and mcsteps, listed according to beta value
mquadravs5     = []
mquadravs_std5 = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av = float(words[1])
        m_avs5.append(m_av)
        # stdm5
        stdm = float(words[2])
        stdms5.append(stdm)
        # mquadravs5
        mq = float(words[10])
        mquadravs5.append(mq)
        # mquadravs_std5
        mq_std = float(words[11])
        mquadravs_std5.append(mq_std)  
     
# Remember to close the file
infile.close()    

infile = open("beta0to5_Nbetas100_spin0p5_L4_J1_mcsteps1000_bins100_091216_1340_IsingMC.txt", "r")

# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)   # Tested, and this is done successfully
L = sqrt(N);

# Getting lists ready to store the data
m_avs4         = []          # Average m over all bins and mcsteps, listed according to beta value
stdms4         = []          # Average m over all bins and mcsteps, listed according to beta value
mquadravs4     = []
mquadravs_std4 = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av = float(words[1])
        m_avs4.append(m_av)
        # stdm4
        stdm4 = float(words[2])
        stdms4.append(stdm)
        # mquadravs4
        mq = float(words[10])
        mquadravs4.append(mq)
        # mquadravs_std4
        mq_std = float(words[11])
        mquadravs_std4.append(mq_std)  
 
# Remember to close the file
infile.close()


infile = open("beta0to5_Nbetas100_spin0p5_L3_J1_mcsteps1000_bins100_091216_1340_IsingMC.txt", "r")

# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)   # Tested, and this is done successfully
L = sqrt(N);

# Getting lists ready to store the data
m_avs3         = []          # Average m over all bins and mcsteps, listed according to beta value
stdms3         = []          # Average m over all bins and mcsteps, listed according to beta value
mquadravs3     = []
mquadravs_std3 = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av = float(words[1])
        m_avs3.append(m_av)
        # stdm3
        stdm = float(words[2])
        stdms3.append(stdm)
        # mquadravs3
        mq = float(words[10])
        mquadravs3.append(mq)
        # mquadravs_std3
        mq_std = float(words[11])
        mquadravs_std3.append(mq_std)  
 
# Remember to close the file
infile.close()


infile = open("beta0to5_Nbetas100_spin0p5_L2_J1_mcsteps1000_bins100_091216_1340_IsingMC.txt", "r")

# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)   # Tested, and this is done successfully
L = sqrt(N);

# Getting lists ready to store the data
m_avs2         = []          # Average m over all bins and mcsteps, listed according to beta value
stdms2         = []          # Average m over all bins and mcsteps, listed according to beta value
mquadravs2     = []
mquadravs_std2 = []

# Read the rest of the lines
lines = infile.readlines()  # This works, I have checked it

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av = float(words[1])
        m_avs2.append(m_av)
        # stdm
        stdm = float(words[2])
        stdms2.append(stdm)
        # mquadravs2
        mq = float(words[10])
        mquadravs2.append(mq)
        # mquadravs_std2
        mq_std = float(words[11])
        mquadravs_std2.append(mq_std)  
 
# Remember to close the file
infile.close()

betas           = array(betas)
m_avs20         = array(m_avs20)
stdms20          = array(stdms20)
mquadravs20     = array(mquadravs20)
mquadravs_std20 = array(mquadravs_std20)
m_avs10         = array(m_avs10)
stdms10          = array(stdms10)
mquadravs10     = array(mquadravs10)
mquadravs_std10 = array(mquadravs_std10)
m_avs5          = array(m_avs5)
stdms5           = array(stdms5)
mquadravs5      = array(mquadravs5)
mquadravs_std5  = array(mquadravs_std5)
m_avs4          = array(m_avs4)
stdms4           = array(stdms4)
mquadravs4      = array(mquadravs4)
mquadravs_std4  = array(mquadravs_std4)
m_avs3          = array(m_avs3)
stdms3           = array(stdms3)
mquadravs3      = array(mquadravs3)
mquadravs_std3  = array(mquadravs_std3)
m_avs2          = array(m_avs2)
stdms2           = array(stdms2)
mquadravs2      = array(mquadravs2)
mquadravs_std2  = array(mquadravs_std2)

length = len(betas)
UL20 = zeros(length)
UL15 = zeros(length)  
UL10 = zeros(length) 
UL5 = zeros(length)
UL4 = zeros(length)
UL3 = zeros(length)
UL2 = zeros(length)
UL20delta = zeros(length)
UL15delta = zeros(length)  
UL10delta = zeros(length) 
UL5delta = zeros(length)
UL4delta = zeros(length)
UL3delta = zeros(length)
UL2delta = zeros(length)
   
for i in range(length):
    UL20[i] = 1 - mquadravs20[i]/(3*m_avs20[i]**2)
    UL15[i] = 1 - mquadravs15[i]/(3*m_avs15[i]**2)
    UL10[i] = 1 - mquadravs10[i]/(3*m_avs10[i]**2)
    UL5[i]  = 1 - mquadravs5[i]/(3*m_avs5[i]**2)
    UL4[i]  = 1 - mquadravs4[i]/(3*m_avs4[i]**2)
    UL3[i] = 1 - mquadravs3[i]/(3*m_avs3[i]**2)
    UL2[i] = 1 - mquadravs2[i]/(3*m_avs2[i]**2)
    UL20delta[i] = UL20[i]*sqrt((mquadravs_std20[i]/mquadravs20[i])**2+(stdms20[i]/m_avs20[i])**2)
    UL15delta[i] = UL15[i]*sqrt((mquadravs_std15[i]/mquadravs15[i])**2+(stdms15[i]/m_avs15[i])**2)
    UL10delta[i] = UL10[i]*sqrt((mquadravs_std10[i]/mquadravs10[i])**2+(stdms10[i]/m_avs10[i])**2)
    UL5delta[i] = UL5[i]*sqrt((mquadravs_std5[i]/mquadravs5[i])**2+(stdms5[i]/m_avs5[i])**2)
    UL4delta[i] = UL4[i]*sqrt((mquadravs_std4[i]/mquadravs4[i])**2+(stdms4[i]/m_avs4[i])**2)
    UL3delta[i] = UL3[i]*sqrt((mquadravs_std3[i]/mquadravs3[i])**2+(stdms3[i]/m_avs3[i])**2)
    UL2delta[i] = UL2[i]*sqrt((mquadravs_std2[i]/mquadravs2[i])**2+(stdms2[i]/m_avs2[i])**2)
    
    
    

#UL20 = 1 - mquadravs20/(3*m_av20**2)

figure()
plot(betas, UL3, label='L=3')
hold('on')
#errorbar(betas, UL2, yerr=UL2delta, label='L=2')
plot(betas, UL4, label='L=4')
plot(betas, UL5, label='L=5')
plot(betas, UL10, label='L=10')
plot(betas, UL15, label='L=15')
plot(betas, UL20, label='L=20')
title(r'The Binder cumulant vs temperature for different L')
xlabel(r'$\beta$')
ylabel(r'$U_L$')
legend(loc="upper right")
show()

figure()
errorbar(betas, UL3, yerr=UL3delta, label='L=3')
hold('on')
#errorbar(betas, UL2, yerr=UL2delta, label='L=2')
errorbar(betas, UL4, yerr=UL4delta, label='L=4')
errorbar(betas, UL5, yerr=UL5delta, label='L=5')
errorbar(betas, UL10, yerr=UL10delta, label='L=10')
errorbar(betas, UL15, yerr=UL15delta, label='L=15')
errorbar(betas, UL20, yerr=UL20delta, label='L=20')
title(r'The Binder cumulant vs temperature for different L')
xlabel(r'$\beta$')
ylabel(r'$U_L$')
legend(loc="upper right")
show()
