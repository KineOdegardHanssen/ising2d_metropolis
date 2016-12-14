from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# Open the file for reading
# 20
infile = open("beta0to5_Nbetas100_spin0p5_L20_J1_mcsteps1000_bins100_IsingMC.txt", "r")


# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only

# Number of spins, spin size, number of bins
N, s, bins  = firstline.split()
N = int(N); s = float(s); bins = int(bins)   # Tested, and this is done successfully
L = sqrt(N);

# Getting lists ready to store the data
betas     = []          # List of beta values
m_avs20   = []          # Average m over all bins and mcsteps, listed according to beta value
stdms20   = []          # Average m over all bins and mcsteps, listed according to beta value
eavs20     = []
eavs_std20 = []
cvs20      = []
#cv_avfromblocks_std20 = []

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
        # blockvars
        stdm20 = float(words[2])
        stdms20.append(stdm20)
        # cvs
        cv = float(words[3])
        cvs20.append(cv)
        # eavs
        eav = float(words[4])/400.0
        eavs20.append(eav)
        # esqavs
        #esqav = float(words[5])
        #esqavs.append(esqav)    
        # eavs_std
        eav_std = float(words[6])/400.0
        eavs_std20.append(eav_std)
        # esqavs_std
        #esqav_std = float(words[7])
        #esqavs_std.append(esqav_std) # Not sure I want these 
        # cv_avfromblocks_std
        #cvblock_std = float(words[8])
        #cv_avfromblocks_std20.append(cvblock_std)  
        
 
# Remember to close the file
infile.close()

# 10        
infile = open("beta0to5_Nbetas100_spin0p5_L10_J1_mcsteps1000_bins100_IsingMC.txt", "r")


# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only


# Getting lists ready to store the data
m_avs10 = []          # Average m over all bins and mcsteps, listed according to beta value
stdms10 = []          # Average m over all bins and mcsteps, listed according to beta value
eavs10     = []
eavs_std10 = []
cvs10      = []
#cv_avfromblocks_std10 = []

# Read the rest of the lines
lines = infile.readlines()  

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av10 = float(words[1])
        m_avs10.append(m_av10)
        # blockvars
        stdm10 = float(words[2])
        stdms10.append(stdm10)
        # cvs
        cv = float(words[3])
        cvs10.append(cv)
        # eavs
        eav = float(words[4])/100.0
        eavs10.append(eav)
        # esqavs
        #esqav = float(words[5])
        #esqavs.append(esqav)    
        # eavs_std
        eav_std = float(words[6])/100.0
        eavs_std10.append(eav_std)
        # esqavs_std
        #esqav_std = float(words[7])
        #esqavs_std.append(esqav_std) # Not sure I want these 
        #cvblock_std = float(words[8])
        #cv_avfromblocks_std10.append(cvblock_std)    
        
        
# Remember to close the file
infile.close()

# 15
infile = open("beta0to5_Nbetas100__spin0p5_L15_J1_mcsteps1000_bins100_IsingMC.txt", "r")


# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only


# Getting lists ready to store the data
m_avs15 = []          # Average m over all bins and mcsteps, listed according to beta value
stdms15 = []          # Average m over all bins and mcsteps, listed according to beta value
eavs15     = []
eavs_std15 = []
cvs15      = []
#cv_avfromblocks_std15 = []
# Read the rest of the lines
lines = infile.readlines()  

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av15 = float(words[1])
        m_avs15.append(m_av15)
        # blockvars
        stdm15 = float(words[2])
        stdms15.append(stdm15)
        # cvs
        cv = float(words[3])
        cvs15.append(cv)
        # eavs
        eav = float(words[4])/225.0
        eavs15.append(eav)
        # esqavs
        #esqav = float(words[5])
        #esqavs.append(esqav)    
        # eavs_std
        eav_std = float(words[6])/225.0
        eavs_std15.append(eav_std)
        # esqavs_std
        #esqav_std = float(words[7])
        #esqavs_std.append(esqav_std) # Not sure I want these   
        #cvblock_std = float(words[8])
        #cv_avfromblocks_std15.append(cvblock_std)  
        
        
# Remember to close the file
infile.close()
        
        
# 30
infile = open("beta0to5_Nbetas100_spin0p5_L30_J1_mcsteps1000_bins100_141216_1152_IsingMC.txt", "r")


# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only


# Getting lists ready to store the data
m_avs30 = []          # Average m over all bins and mcsteps, listed according to beta value
stdms30 = []          # Average m over all bins and mcsteps, listed according to beta value
eavs30     = []
eavs_std30 = []
cvs30      = []
#cv_avfromblocks_std4 = []

# Read the rest of the lines
lines = infile.readlines()  

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av = float(words[1])
        m_avs30.append(m_av)
        # blockvars
        stdm = float(words[2])
        stdms30.append(stdm)
        # cvs
        cv = float(words[3])
        cvs30.append(cv)
        # eavs
        eav = float(words[4])/16.0
        eavs30.append(eav)
        # esqavs
        #esqav = float(words[5])
        #esqavs.append(esqav)    
        # eavs_std
        eav_std = float(words[6])/16.0
        eavs_std30.append(eav_std)
        # esqavs_std
        #esqav_std = float(words[7])
        #esqavs_std.append(esqav_std) # Not sure I want these   
        #cvblock_std = float(words[8])
        #cv_avfromblocks_std4.append(cvblock_std)  
        
        
# Remember to close the file
infile.close()

# We prefer arrays
betas      = array(betas)
m_avs10    = array(m_avs10)
stdms10    = array(stdms10)
m_avs15    = array(m_avs15)
stdms15    = array(stdms15)
m_avs20    = array(m_avs20)
stdms20    = array(stdms20)
m_avs30    = array(m_avs30)
stdms30    = array(stdms30)
eavs10     = array(eavs10)
eavs_std10 = array(eavs_std10)
eavs15     = array(eavs15)
eavs_std15 = array(eavs_std15)
eavs20     = array(eavs20)
eavs_std20 = array(eavs_std20)
eavs30      = array(eavs30)
eavs_std30  = array(eavs_std30)
cvs30      = array(cvs30)
cvs20      = array(cvs20)
#cv_avfromblocks_std20 = array(cv_avfromblocks_std20)
cvs15      = array(cvs15)
#cv_avfromblocks_std15 = array(cv_avfromblocks_std15)
cvs10      = array(cvs10)


#Midvandering
eps10 = 1000
eps15 = 1000
eps20 = 1000
eps20 = 1000
eps30 = 1000
Ls = array([10,15,20, 30])
betacrit = zeros(len(Ls))
for i in range(len(m_avs2)):
    if(abs(m_avs10[i]-0.125)<eps10):
        betacrit[0] = betas[i]
        eps10 = abs(m_avs10[i]-0.125)
    if(abs(m_avs15[i]-0.125)<eps15):
        betacrit[1] = betas[i]
        eps15 = abs(m_avs15[i]-0.125)
    if(abs(m_avs20[i]-0.125)<eps20):
        betacrit[2] = betas[i]
        eps20 = abs(m_avs20[i]-0.125)
    if(abs(m_avs30[i]-0.125)<eps30):
        betacrit[3] = betas[i]
        eps30 = abs(m_avs30[i]-0.125)
         
exactbetacscalar = 2*log(1+sqrt(2))
exactbetac = zeros(len(Ls))+exactbetacscalar

def extrapolation(y1, y2, x1, x2, x):
    return y1+(x-x1)/(x2-x1)*(y2-y1)

extrapolated_betacrit = extrapolation(betacrit[2], betacrit[3], Ls[2], Ls[3], 1e7)
print "betacrit[5] = ", betacrit[1]
print "betacrit[6] = ", betacrit[2]
print "betacrit[7] = " , betacrit[3]
print "Exact value of beta_c: ", exactbetacscalar
print "Extrapolated value (when L=1e7): ", extrapolated_betacrit

figure()
plot(Ls, betacrit)
hold('on')
plot(Ls, exactbetac)
title(r'Value of $\beta$ where $<m^2>=0.125$ for different values of L')
xlabel('L')
ylabel(r'$\beta$')
show()
