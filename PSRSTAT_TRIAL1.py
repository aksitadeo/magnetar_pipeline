#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 11:37:00 2022

@author: deo010
"""

#%% Folder destination

print("""
      
      For this script, go to the ARCHIVEs/ folder destination containing the .ar files.
      
      """)


#%% Importing modules

import os
import subprocess as sp
import pandas as pd

#%% BEFORE RUNNING THIS SCRIPT

print("""
      
      Before running this script, it is important to conduct source-specific RFI zapping. This can be done
in the terminal window of the system itself, using the following command as a template:
   
    paz -F "2670 2690" -F "3440 3465" -F "3545 3570" <filename>.ar -e .ar.zap
    
where the numbers* are the rows from which you want to delete RFI. Then, we should average in frequency
using the terminal command

    pam -D -F *ARCH.ar.zap -e .zap.F 
    
    *Tip: use the command 'pazi <filename>.ar' to output the row numbers you want, by left click selecting rows, right clicking, then pressing 'p'.
    
    """)

# J1550-5418 used paz -F "2670 2690" -F "3440 3465" -F "3545 3570" *.ar -e .ar.zap
# J1622-4950 used paz -F "2630 2670" -F "3444 3472" -F "3542 3570" *.ar -e .ar.zap


#%% Main Loop

x =[]

for filename in os.listdir():
    if filename.endswith(".ar.zap.F"):
        print("")
        print('Now working on: {}'.format(filename))
        output = sp.getoutput('psrstat -Q -c snr {}'.format(filename))
        xs = output.split()
        x += [xs]

#%% Convert to float

for thing in x:
    thing[1] = float(thing[1])
    print(thing)
print(x)       

#%% Create dataframe

df = pd.DataFrame(x)
df.columns = ['filename','S/N']
final_df = df.sort_values(by=['S/N'], ascending=False)
final_df

#%% Show the top 10 highest S/N candidates

top10 = final_df.head(20)

#%% View the sources

source = 'J1818-1607'

for file in top10['filename']:
    files = file[0:17]
    print("Now working on: {}".format(files))
    # /DATA/MENSA_1/deo010/dir/J1550-5418/F_s130416_212545.sf
    os.system('gv /DATA/MENSA_1/deo010/dir/{}/F_{}/{}TEST_rfifind.ps'.format(source,files,files))

#%% Averaging in polarisation

for filename in os.listdir():
    if filename.endswith(".ar.zap.F"):
        print("")
        print('Now working on: {}'.format(filename))
        os.system('pam -p {} -e Fp'.format(filename))

#%% Change to text file

for filename in os.listdir():
    if filename.endswith(".ar.zap.Fp"):
        print("")
        print('Now working on: {}'.format(filename))
        os.system('pdv -t {} >> {}.txt'.format(filename,filename))

#%% Code testing

print('test')

os.system('psrplot ARCHIVEs/t160114_192200.sfARCH.ar.zap.F')






