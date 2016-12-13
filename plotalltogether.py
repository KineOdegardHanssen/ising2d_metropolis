from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# Open the file for reading
#infile = open("beta0to5_Nbetas100__spin0p5_L15_J1_mcsteps1000_bins100_IsingMC.txt", "r")
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

infile = open("beta0to5_Nbetas100_spin0p5_L2_J1_mcsteps1000_bins100_IsingMC.txt", "r")


# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only


# Getting lists ready to store the data
m_avs2 = []          # Average m over all bins and mcsteps, listed according to beta value
stdms2 = []          # Average m over all bins and mcsteps, listed according to beta value
eavs2     = []
eavs_std2 = []
cvs2      = []

# Read the rest of the lines
lines = infile.readlines()  

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av2 = float(words[1])
        m_avs2.append(m_av2)
        # blockvars
        stdm2 = float(words[2])
        stdms2.append(stdm2)
        # cvs
        cv = float(words[3])
        cvs2.append(cv)
        # eavs
        eav = float(words[4])/4.0
        eavs2.append(eav)
        # esqavs
        #esqav = float(words[5])
        #esqavs.append(esqav)    
        # eavs_std
        eav_std = float(words[6])/4.0
        eavs_std2.append(eav_std)
        # esqavs_std
        #esqav_std = float(words[7])
        #esqavs_std.append(esqav_std) # Not sure I want these   
        #cvblock_std = float(words[8])
        #cv_avfromblocks_std2.append(cvblock_std)  
        
        
# Remember to close the file
infile.close()

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
        
        
infile = open("beta0to5_Nbetas100_spin0p5_L5_J1_mcsteps1000_bins100_IsingMC.txt", "r")


# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only


# Getting lists ready to store the data
m_avs5 = []          # Average m over all bins and mcsteps, listed according to beta value
stdms5 = []          # Average m over all bins and mcsteps, listed according to beta value
eavs5     = []
eavs_std5 = []
cvs5      = []
cv_avfromblocks_std5 = []

# Read the rest of the lines
lines = infile.readlines()  

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av5 = float(words[1])
        m_avs5.append(m_av5)
        # blockvars
        stdm5 = float(words[2])
        stdms5.append(stdm5)
        # cvs
        cv = float(words[3])
        cvs5.append(cv)
        # eavs
        eav = float(words[4])/25.0
        eavs5.append(eav)
        # esqavs
        #esqav = float(words[5])
        #esqavs.append(esqav)    
        # eavs_std
        eav_std = float(words[6])/25.0
        eavs_std5.append(eav_std)
        # esqavs_std
        #esqav_std = float(words[7])
        #esqavs_std.append(esqav_std) # Not sure I want these 
        #cvblock_std = float(words[8])
        #cv_avfromblocks_std5.append(cvblock_std)  
          
        
        
# Remember to close the file
infile.close()

infile = open("beta0to5_Nbetas100_spin0p5_L3_J1_mcsteps1000_bins100_IsingMC.txt", "r")


# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only


# Getting lists ready to store the data
m_avs3 = []          # Average m over all bins and mcsteps, listed according to beta value
stdms3 = []          # Average m over all bins and mcsteps, listed according to beta value
eavs3     = []
eavs_std3 = []
cvs3      = []
#cv_avfromblocks_std3 = []

# Read the rest of the lines
lines = infile.readlines()  

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av3 = float(words[1])
        m_avs3.append(m_av3)
        # blockvars
        stdm3 = float(words[2])
        stdms3.append(stdm3)
        # cvs
        cv = float(words[3])
        cvs3.append(cv)
        # eavs
        eav = float(words[4])/9.0
        eavs3.append(eav)
        # esqavs
        #esqav = float(words[5])
        #esqavs.append(esqav)    
        # eavs_std
        eav_std = float(words[6])/9.0
        eavs_std3.append(eav_std)
        # esqavs_std
        #esqav_std = float(words[7])
        #esqavs_std.append(esqav_std) # Not sure I want these 
        #cvblock_std = float(words[8])
        #cv_avfromblocks_std3.append(cvblock_std)    
        
        
# Remember to close the file
infile.close()

infile = open("beta0to5_Nbetas100_spin0p5_L4_J1_mcsteps1000_bins100_IsingMC.txt", "r")


# The first line contains information about the system. We read that separately
firstline = infile.readline() # Reads the first line only


# Getting lists ready to store the data
m_avs4 = []          # Average m over all bins and mcsteps, listed according to beta value
stdms4 = []          # Average m over all bins and mcsteps, listed according to beta value
eavs4     = []
eavs_std4 = []
cvs4      = []
#cv_avfromblocks_std4 = []

# Read the rest of the lines
lines = infile.readlines()  

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        # m_avs
        m_av4 = float(words[1])
        m_avs4.append(m_av4)
        # blockvars
        stdm4 = float(words[2])
        stdms4.append(stdm4)
        # cvs
        cv = float(words[3])
        cvs4.append(cv)
        # eavs
        eav = float(words[4])/16.0
        eavs4.append(eav)
        # esqavs
        #esqav = float(words[5])
        #esqavs.append(esqav)    
        # eavs_std
        eav_std = float(words[6])/16.0
        eavs_std4.append(eav_std)
        # esqavs_std
        #esqav_std = float(words[7])
        #esqavs_std.append(esqav_std) # Not sure I want these   
        #cvblock_std = float(words[8])
        #cv_avfromblocks_std4.append(cvblock_std)  
        
        
# Remember to close the file
infile.close()


# We prefer arrays
betas      = array(betas)
m_avs2     = array(m_avs2)
stdms2     = array(stdms2)
m_avs3     = array(m_avs3)
stdms3     = array(stdms3)
m_avs4     = array(m_avs4)
stdms4     = array(stdms4)
m_avs5     = array(m_avs5)
stdms5     = array(stdms5)
m_avs10    = array(m_avs10)
stdms10    = array(stdms10)
m_avs15    = array(m_avs15)
stdms15    = array(stdms15)
m_avs20    = array(m_avs20)
stdms20    = array(stdms20)
eavs2      = array(eavs2)
eavs_std2  = array(eavs_std2)
eavs3      = array(eavs3)
eavs_std3  = array(eavs_std3)
eavs4      = array(eavs4)
eavs_std4  = array(eavs_std4)
eavs5      = array(eavs5)
eavs_std5  = array(eavs_std5)
eavs10     = array(eavs10)
eavs_std10 = array(eavs_std10)
eavs15     = array(eavs15)
eavs_std15 = array(eavs_std15)
eavs20     = array(eavs20)
eavs_std20 = array(eavs_std20)
cvs20      = array(cvs20)
#cv_avfromblocks_std20 = array(cv_avfromblocks_std20)
cvs15      = array(cvs15)
#cv_avfromblocks_std15 = array(cv_avfromblocks_std15)
cvs10      = array(cvs10)
#cv_avfromblocks_std10 = array(cv_avfromblocks_std10)
cvs5      = array(cvs5)
#cv_avfromblocks_std5 = array(cv_avfromblocks_std5)
cvs4      = array(cvs4)
#cv_avfromblocks_std4 = array(cv_avfromblocks_std4)
cvs3      = array(cvs3)
#cv_avfromblocks_std3 = array(cv_avfromblocks_std3)
cvs2      = array(cvs2)
#cv_avfromblocks_std2 = array(cv_avfromblocks_std2)
"""
cvs = array(cvs)
eavs = array(eavs)
esqavs = array(esqavs)
eavs_std = array(eavs_std)
esqavs_std = array(esqavs_std)"""


figure()
plot(betas, m_avs2, label='L=2')
hold('on')
plot(betas, m_avs3, label='L=3')
plot(betas, m_avs4, label='L=4')
plot(betas, m_avs5, label='L=5')
plot(betas, m_avs10, label='L=10')
plot(betas, m_avs15, label='L=15')
plot(betas, m_avs20, label='L=20')
title('Average squared magnetization $<m>^2$ in the Ising model for different L')
xlabel(r'$\beta$')
ylabel(r'$<m>^2$')
legend(loc="upper right")
show()

figure()
errorbar(betas, m_avs2, yerr=stdms2, label='L=2')
hold('on')
errorbar(betas, m_avs3, yerr=stdms3, label='L=3')
errorbar(betas, m_avs4, yerr=stdms4, label='L=4')
errorbar(betas, m_avs5, yerr=stdms5, label='L=5')
errorbar(betas, m_avs10, yerr=stdms10, label='L=10')
errorbar(betas, m_avs15, yerr=stdms15, label='L=15')
errorbar(betas, m_avs20, yerr=stdms20, label='L=20')
title('$<m>^2$ in the Ising model for different L, with errorbars')
xlabel(r'$\beta$')
ylabel(r'$<m>^2$')
legend(loc="upper right")
show()

figure()
plot(betas, eavs2, label='L=2')
hold('on')
plot(betas, eavs3, label='L=3')
plot(betas, eavs4, label='L=4')
plot(betas, eavs5, label='L=5')
plot(betas, eavs10, label='L=10')
plot(betas, eavs15, label='L=15')
plot(betas, eavs20, label='L=20')
title('Energy $<E>$ per particle in the Ising model for different L')
xlabel(r'$\beta$')
ylabel(r'$<E>$')
legend(loc="upper right")
show()

figure()
errorbar(betas, eavs2, yerr=eavs_std2, label='L=2')
hold('on')
errorbar(betas, eavs3, yerr=eavs_std3, label='L=3')
errorbar(betas, eavs4, yerr=eavs_std4, label='L=4')
errorbar(betas, eavs5, yerr=eavs_std5, label='L=5')
errorbar(betas, eavs10, yerr=eavs_std10, label='L=10')
errorbar(betas, eavs15, yerr=eavs_std15, label='L=15')
errorbar(betas, eavs20, yerr=eavs_std20, label='L=20')
title('Energy $<E>$ per particle in the Ising model for different L, with errorbars')
xlabel(r'$\beta$')
ylabel(r'$<E>$')
legend(loc="upper right")
show()

figure()
plot(betas, cvs2, label='L=2')
hold('on')
plot(betas, cvs3, label='L=3')
plot(betas, cvs4, label='L=4')
plot(betas, cvs5, label='L=5')
plot(betas, cvs10, label='L=10')
plot(betas, cvs15, label='L=15')
plot(betas, cvs20, label='L=20')
title(r'Total heat capacity $C_V$ in the Ising model for different L')
xlabel(r'$\beta$')
ylabel(r'$C_V$')
legend(loc="upper right")
show()

'''
figure()
plot(betas, cvs2)
hold('on')
title(r'Total heat capacity $C_V$ in the Ising model for L=2')
xlabel(r'$\beta$')
ylabel(r'$C_V$')
show()

figure()
plot(betas, cvs3)
title(r'Total heat capacity $C_V$ in the Ising model for L=3')
xlabel(r'$\beta$')
ylabel(r'$C_V$')
show()

figure()
plot(betas, cvs4)
title(r'Total heat capacity $C_V$ in the Ising model for L=4')
xlabel(r'$\beta$')
ylabel(r'$C_V$')
show()

figure()
plot(betas, cvs5)
title(r'Total heat capacity $C_V$ in the Ising model for L=5')
xlabel(r'$\beta$')
ylabel(r'$C_V$')
show()

figure()
plot(betas, cvs10)
title(r'Total heat capacity $C_V$ in the Ising model for L=10')
xlabel(r'$\beta$')
ylabel(r'$C_V$')
show()

figure()
plot(betas, cvs15)
title(r'Total heat capacity $C_V$ in the Ising model for L=15')
xlabel(r'$\beta$')
ylabel(r'$C_V$')
show()

figure()
plot(betas, cvs20)
title(r'Total heat capacity $C_V$ in the Ising model for L=20')
xlabel(r'$\beta$')
ylabel(r'$C_V$')
show()

'''


# Example from before:
"""
plot(t, qA_NA, label="In A")
xlabel("Time")
ylabel("Energy units")
hold("on")
plot(t, qB_NB, label="In B")
legend(loc="upper right")
title("Energy per oscillator qi/Ni")
show()
"""


