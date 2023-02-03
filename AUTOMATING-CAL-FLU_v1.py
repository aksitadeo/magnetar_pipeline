#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 14:38:40 2023

@author: deo010
"""

#%% Folder destination

print("""
      
      For section of the script, go to the ARCHIVEs/ folder destination containing the .ar files.
      
      """)

#%% Importing modules

import numpy as np
import math
import os
from matplotlib import pyplot as plt

#%% Plotting graphs

def plotgraphs(filename):
    
    file = np.genfromtxt("{}".format(filename), skip_header=1, delimiter=" ", dtype=float)[:,2:].T
    phase = file[0,:]
    flux = file[1,:]
    
    plt.figure(figsize=(20,12))
    plt.style.use('dark_background')
    plt.plot(phase,flux,'w')
    plt.xlabel('Time (ms)')
    plt.ylabel('Flux Density (mJy)')
    plt.title('{}'.format(filename))


#%% Find the noisy channels

def noise(filename):
    file = np.genfromtxt("{}".format(filename), skip_header=1, delimiter=" ", dtype=float)[:,2:].T
    phase = file[0,:]
    flux = file[1,:]
    
    length = len(phase)
    stds = []
    stdbins = []
    bins = 3
    
    for i in range(bins):                                                         
        fluxbins = flux[(i*int(length/bins)):((i+1)*int(length/bins))]         # create bins
        phasebins = phase[(i*int(length/bins)):((i+1)*int(length/bins))]       # pick out bins of phase    
        std = np.std(fluxbins)                                                 # take std of each bin
                                       
        stds += [std]
        stdbins += [phasebins]                  
        
    return stdbins[np.argmin(stds)]

#%% Calculate experimental rms

def rmsprac(filename,noise):
    
    ######### Finding RMS #########
    
    # the equation used is RMS = sqrt [ 1/n * sum (x_i) for i values ]  
    
    file = np.genfromtxt("{}".format(filename), skip_header=1, delimiter=" ", dtype=float)[:,3:].T
    I = file[0,:]
    
    start = int(min(noise))
    end = int(max(noise))
    obs = end-start
    
    rmsp = math.sqrt( (1/obs)*np.sum( (I[start:end])**2 ) )
    
    print("The practical rms is {}".format(rmsp))
    
    return rmsp, obs

#%% Calculate theoretical rms

def rmstheory(beta, SEFD, np, t_int, del_f):
    
    ######### Calculating RMS #########
    
    # Equation used is radiometer equation
    
    beta = float(beta)                                                         # correction factor
    SEFD = float(SEFD)                                                         # milliJansky, calculated using a fold-mode observation of J1550-5418 taken on 2015-03-12-22:14:02 UTC
    np = float(np)                                                             # polarisation number summed
    t_int = float(t_int)                                                       # integration time per phase bin
    del_f = float(del_f)                                                       # bandwidth
    
    rmst = (beta*SEFD)/(math.sqrt(np*t_int*del_f))
      
    print("The theoretical rms is {}".format(rmst))
    
    return rmst

#%% Find the correction factor

def correction(filename, rmst, rmsp):
    factor = rmst/rmsp
    print("The correction factor is {}.".format(factor))
    os.system("awk '{{if (NR!=1) {{{{ $4=$4*{} }}}} print }}' {} > cal_{}".format(factor,filename,filename))

#%% Find a window for the width of the fluence

def find_nearest(array,value):
     # define array and index of array elements
    array = np.asarray(array)
    idx = np.searchsorted(array, value, side="left")
    
    # if edge case
    if idx == len(array):                                                
        idx = idx-1
    
    # if left side of the peak
    if array[idx] < value:
       # return value of left and right side of the peak
        try:
            return array[idx], array[idx+1]
       # if no right side of peak, return edge value
        except IndexError:
            return array[idx], 1500
    # if not
    else:
        if idx == 0:
            return 0, array[idx]
        else:
            return array[idx-1],array[idx]

#%% Check if the value is less

def checkless(list1, val):
    index = []
    value = []
    
    for i,x in enumerate(list1):
        if x < val:
            index.append(i)
            value.append(x)
    
    return index, value
    

#%% Finding troughs
def troughfinds(filename,sig):
    file = np.genfromtxt("cal_{}".format(filename), skip_header=1, delimiter=" ", dtype=float)[:,2:].T
    flux = file[1,:]
    
    minind,minval = checkless(flux,-sig)
    maxi = max(flux)
    arr = np.argmax(flux)
    
    SN = maxi/(2*sig)
    print("The Signal-to-Noise ratio is {}".format(SN))
    minim, maxim = find_nearest(minind,arr)

    return minim, maxim, SN


#%% Find fluence

def window(filename,start,end):
    file = np.genfromtxt("cal_{}".format(filename), skip_header=1, delimiter=" ", dtype=float)[:,2:].T
    phase = file[0,:]
    flux = file[1,:]
    print(np.max(flux))
    
    # find burst window
    burst = end-start
    fluence = (sum(flux))/1000                                    # convert to Jy from mJy
    
    # plot the window
    plt.figure(figsize=(20,12))
    plt.style.use('dark_background')
    plt.plot(phase,flux,'w')
    plt.axvline(start,color='red')
    plt.axvline(end,color='red')
    plt.xlabel('Time (ms)')
    plt.ylabel('Flux Density (mJy)')
    plt.title('{}'.format(filename))
    
    print("Burst width is {} ms".format(burst))
    print("Fluence is {} Jy ms".format(fluence))
    
    return fluence

#%% Big function
def cal_fluence(filename):
    
    ##### CALIBRATION #####
    
    # Change these values when changing between sources
    
    beta = 1.1                  # correction factor
    SEFD = 25000                # miliJansky
    np = 2                      # polarisation number summed
    t = 1.5                     # integration time per phase bin
    del_f = 256e6               # bandwidth
    
    
    # Find which channels do not contain the peak
    noisetime = noise(filename)
    # Find the practical RMS 
    rmsp,obs = rmsprac(filename,noisetime)
    # Find the theoretical RMS
    rmst = rmstheory(beta, SEFD, np, t/obs, del_f)
    # Create text file with calibrated data
    correction(filename, rmst,rmsp)
    # Plot calibrated plot
    #plotgraphs('cal_{}'.format(filename))
    
    ##### FLUENCE #####
    
    # Find the peak window
    start, end, SN = troughfinds(filename,rmst)
    # Find the fluence of the burst
    if SN > 2.5:
        fluence = window(filename,start,end)
        return fluence
    
    else:
        print("Check {} for bursts".format(filename))

#%% Example code to run through all data

# creates a list of source that had errors to check for bursts manually
Errors = []
# creates a list of fluences
Fluences = []

# looping over all observations
for filename in os.listdir():
    if filename.endswith(".txt"):
        try: 
            print(" ")
            print("Now working on: {}".format(filename))
            fluence,flux = cal_fluence('{}'.format(filename))
            Fluences += [fluence]
        except Exception:
            print("Check file {}".format(filename))
            Errors += [filename]

#%% Save fluences files
# change this with each source and adjust directory as needed
source = 'J1818-1607'
np.savetxt('/DATA/MENSA_1/deo010/dir/FLUENCE_{}.txt'.format(source),Fluences,fmt='%s')

#%% Folder destination

print("""
      
      For section of the script, go your ORIGINAL directory containing the subdirectories holding the .sf files for each source.
      
      """)

#%% General notice

print("""
      
      TRANSIENT PHASE SPACE DIAGRAM
      
      The following code is specific to a project conducted and mentions files that that may not be applicable to what
      you are trying to achieve. If you are trying to plot a transient phase space diagram, change source variables and 
      file names to those relevent to your project. The FRBs and magnetar fluence values from literature are cited.
      
      """)
#%% Add our sources

#J1550-5418
J1550 = np.genfromtxt('FLUENCE_J1550-5418.txt',delimiter=',',filling_values=np.nan)
J1550 = J1550[~np.isnan(J1550)]
J1550dist = [4500]*len(J1550)

J1550N = np.genfromtxt('FLUENCE_J1550-5418NEW.txt',delimiter=',',filling_values=np.nan)
J1550N = J1550N[~np.isnan(J1550N)]
J1550Ndist = [4500]*len(J1550N)

#J1622-4950
J1622 = np.genfromtxt('FLUENCE_J1622-4950.txt',delimiter=',',filling_values=np.nan)
J1622 = J1622[~np.isnan(J1622)]
J1622dist = [5571]*len(J1622)

J1622N = np.genfromtxt('FLUENCE_J1622-4950NEW.txt',delimiter=',',filling_values=np.nan)
J1622N = J1622N[~np.isnan(J1622N)]
J1622Ndist = [5571]*len(J1622N)

#J1818-1607
J1818 = np.genfromtxt('FLUENCE_J1818-1607.txt',delimiter=',',filling_values=np.nan)
J1818 = J1818[~np.isnan(J1818)]
J1818dist = [7500]*len(J1818)


#%% Distance formula using Hubble's law

def z(z):
    c = 300000 #km/s
    v = z*c #km/s
    H0 = 67.3 #km/s/Mpc
    d = v/H0
    return [d*1e6]

#%% Literature sources

# =============================================================================
# ASKAP FRBs
# =============================================================================

# https://arxiv.org/pdf/2109.11535.pdf , https://iopscience.iop.org/article/10.3847/2041-8213/ac540f/pdf
FRB_20201124A = [125,640,51,21,58,51,27,41,35,29,28,18,13,16,15,39,0.94,3.68,0.61,1.49,1.53,6.5,0.94,4.48,1.4,1.31,5.92,0.63,3.39,6.56,1.98,1.95,1.78,0.91]
d_FRB_20201124A = [374e6]*len(FRB_20201124A)

# https://arxiv.org/pdf/2009.01214.pdf
FRB_20190711A = [1.4]
d_FRB_20190711A = z(0.522)*len(FRB_20190711A)

# https://arxiv.org/pdf/1810.04354.pdf
FRB_171020 = [200]
d_FRB_171020 = z(0.00867)*len(FRB_171020)

# https://arxiv.org/pdf/2005.13161.pdf
FRB_180924 = [16]
d_FRB_180924 = z(0.3214)*len(FRB_180924)

FRB_181112 = [23,26]
d_FRB_181112 = z(0.4755)*len(FRB_181112)

FRB_190102 = [14]
d_FRB_190102 = z(0.291)*len(FRB_190102)

FRB_190608 = [26]
d_FRB_190608 = z(0.1178)*len(FRB_190608)

FRB_190611 = [10]
d_FRB_190611 = z(0.378)*len(FRB_190611)

FRB_190711 = [34]
d_FRB_190711 = z(0.522)*len(FRB_190711)

# https://arxiv.org/pdf/2005.13160.pdf
FRB_121102 = [0.1]
d_FRB_121102 = z(0.19273)*len(FRB_121102)

FRB_190523 = [280]
d_FRB_190523 = z(0.66)*len(FRB_190523)

# https://iopscience.iop.org/article/10.3847/2041-8213/abb462/pdf
FRB_191001 = [143]
d_FRB_191001 = z(0.234)*len(FRB_191001)

# https://arxiv.org/pdf/2211.16790.pdf
FRB_20210117A = [36]
d_FRB_20210117A = z(0.214)*len(FRB_20210117A)

# https://arxiv.org/pdf/2108.01282.pdf
FRB_20180301A = [4.9]
d_FRB_20180301A = z(0.3304)*len(FRB_20180301A)

FRB_20191228A = [40]
d_FRB_20191228A = z(0.2432)*len(FRB_20191228A)

FRB_20200206A = [59]
d_FRB_20200206A = z(0.3688)*len(FRB_20200206A)

# =============================================================================
# EVN and DSA FRBs
# =============================================================================

# B 20200120E
FRB_20200120 = [0.13,0.63,0.53,0.71,0.09]
d_FRB_20200120 = [3.6e6]*len(FRB_20200120)

# https://arxiv.org/pdf/2001.02222.pdf
FRB_180916 = [0.72,0.2,0.62,2.53]
d_FRB_180916 = z(0.033)*len(FRB_180916)

# https://arxiv.org/pdf/2301.01000.pdf
FRB_20220319 = [11.2]
d_FRB_20220319 = [50e6]*len(FRB_20220319)

# https://arxiv.org/pdf/2211.09049.pdf
FRB_20220912A = [70,56,34,10,136,779,237,19,508,369,30]
d_FRB_20220912A = z(0.0771)*len(FRB_20220912A)

# =============================================================================
#  SGR & J1550 Bursts
# =============================================================================

# https://www.nature.com/articles/s41586-020-2863-y
SGR_1935 = [480000,220000]
d_SGR_1935 = [10000]*len(SGR_1935)

# https://arxiv.org/abs/2011.06607
J1550lit = [0.6e3]
d_J1550lit = [4500]*len(J1550lit)

crab = [583,321,210,4742]
d_crab = [1860]*len(crab)

#%% The formula

pc = np.logspace(3,10,1024)*3.086e18
e = [1e26,1e29,1e32,1e35,1e38,1e41,1e44]
Hz = 500e6

#%% Plotting the transient phase space diagram

plt.figure(figsize=(12,10))
plt.style.use('default')

# Lines of constant energy
for x in e:
    F = (1e26)*x*(1/(4*np.pi*(pc**2)*Hz))
    eq = "1e26 erg"
    plt.plot(pc/3e18,F,linestyle='dashed',alpha=0.3,zorder=0)

# Sources

# Our sources
plt.scatter(J1550dist,J1550,marker='s',label='IE 1547.0-5408',color='lightblue',zorder=5)
plt.scatter(J1550Ndist,J1550N,marker='s',color='lightblue',zorder=5)
plt.scatter(J1622dist,J1622,marker='s',label='J1622-4950',color='blue',zorder=5)
plt.scatter(J1622Ndist,J1622N,marker='s',color='blue',zorder=5)
plt.scatter(J1818dist,J1818,marker='s',label='J1818-1907',color='pink',zorder=5)

# Fast Radio Bursts
plt.scatter(d_FRB_20201124A,FRB_20201124A,marker='o',label='FRB 20201124A',zorder=5)
plt.scatter(d_FRB_20190711A,FRB_20190711A,marker='o',label='FRB 20190711A',zorder=5)
plt.scatter(d_FRB_171020,FRB_171020,marker='o',label='FRB 171020',zorder=5)
plt.scatter(d_FRB_180924,FRB_180924,marker='o',label='FRB 180924',zorder=5)
plt.scatter(d_FRB_181112,FRB_181112,marker='o',label='FRB 181112',zorder=5)
plt.scatter(d_FRB_190102,FRB_190102,marker='o',label='FRB 190102',zorder=5)
plt.scatter(d_FRB_190608,FRB_190608,marker='o',label='FRB 190608',zorder=5)
plt.scatter(d_FRB_190611,FRB_190611,marker='o',label='FRB 190611',zorder=5)
plt.scatter(d_FRB_190711,FRB_190711,marker='o',label='FRB 190711',zorder=5)
plt.scatter(d_FRB_121102,FRB_121102,marker='x',label='FRB 121102',zorder=5)
plt.scatter(d_FRB_190523,FRB_190523,marker='x',label='FRB 190523',zorder=5)
plt.scatter(d_FRB_20210117A,FRB_20210117A,marker='x',label='FRB 20210117A',zorder=5)
plt.scatter(d_FRB_20180301A,FRB_20180301A,marker='x',label='FRB 20180301A',zorder=5)
plt.scatter(d_FRB_20191228A,FRB_20191228A,marker='x',label='FRB 20191228A',zorder=5)
plt.scatter(d_FRB_20200206A,FRB_20200206A,marker='x',label='FRB 20200206A',zorder=5)
plt.scatter(d_FRB_20200120,FRB_20200120,marker='x',label='FRB 20200120',zorder=5)
plt.scatter(d_FRB_180916,FRB_180916,marker='x',label='FRB 180916',zorder=5)
plt.scatter(d_FRB_20220319,FRB_20220319,marker='x',label='FRB 20220319',zorder=5)
plt.scatter(d_FRB_20220912A,FRB_20220912A,marker='x',label='FRB 20220912A',zorder=5)
plt.scatter(d_SGR_1935,SGR_1935,marker='*',label='SGR_1935',zorder=5)
plt.scatter(d_J1550lit,J1550lit,marker='*',label='IE 1547.0-5408 Israel et al.',zorder=5)
plt.scatter(d_crab,crab,marker='*',label="Crab Pulsar",zorder=5)

plt.legend(bbox_to_anchor=(1, 1),loc='upper left',fancybox=True)
plt.title("Transient Phase Space",fontsize=15)
plt.xlabel('Distance (pc)',fontsize=12)
plt.ylabel('Flux Density (Jy ms)',fontsize=12)

# lines of constant energy plots
eq = ["1e26 erg","1e29 erg","1e32 erg","1e35 erg","1e38 erg","1e41 erg","1e44 erg"]
plt.text(10**3.5,10**-2, eq[0], fontsize=10,rotation=-45,color='lightblue')
plt.text(10**4.5,10**-1, eq[1], fontsize=10,rotation=-45,color='orange')
plt.text(10**5.5,10**0, eq[2], fontsize=10,rotation=-45,color='lightgreen')
plt.text(10**6.5,10**1, eq[3], fontsize=10,rotation=-45,color='red')
plt.text(10**7.2,10**2.6, eq[4], fontsize=10,rotation=-45,color='purple')
plt.text(10**8.1,10**3.8, eq[5], fontsize=10,rotation=-45,color='brown')
plt.text(10**9,10**4.97, eq[6], fontsize=10,rotation=-45,color='pink')

plt.xscale('log')
plt.yscale('log')
plt.xlim(1e3,1e10)
plt.ylim(1e-3,1e7)

