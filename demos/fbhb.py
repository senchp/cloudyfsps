#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
import sys
import numpy as np
import fsps
from cloudyfsps.ASCIItools import (writeASCII, compileASCII, checkCompiled, compiledExists)
from cloudyfsps.cloudyInputTools import (cloudyInput, printParFile)
from cloudyfsps.generalTools import calcForLogQ
#runMake, formatAllOutput, writeFormattedOutput)

# this code snippet goes through every step needed
# to integrate FSPS into Cloudy.
# This example uses stellar pops with a constant SFH
# as the input ionizing source.
# 1. Write an ascii file in Cloudy format with grid
#    of FSPS spectra in all available ages and
#    metallicities
# 2. Compile asii file into binary format required
#    for Cloudy use. Assumes $CLOUDY_EXE is set to
#    your /path/to/cloudy.exe
# 3. Writes Cloudy input files for a subset of grid
#    parameters.
# 4. Runs Cloudy on the *.in files
# 5. Formats the various output files

zsun = 0.019

exec_write_ascii = True
exec_write_input = True
exec_run_cloudy = False
exec_write_output = False
exec_gen_FSPS_grid = False
make_ocelote = True

# Function to write the ascii file.
# This is where you set the properties of the
# ionizing spectrum (SSP/CSFH, IMF, FBHB, etc)

def fbhb_ascii(fileout, **kwargs):
    # change these parameters to modify the ionizing source grid
    # default mode is to produce an ascii grid in age and Z,
    # though different variables and more dimensions are possible.
    fbhb_fracs = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    sp_dict = dict(zcontinuous=1,
                   imf_type=2,
                   sfh=0)
    sp = fsps.StellarPopulation(**sp_dict)
    # all ages, default solar metallicity
    ages = 10.**sp.log_age
    modpars = [(age, fb) for age in ages for fb in fbhb_fracs]
    lam = sp.wavelengths
    all_fluxs = []
    for fb in fbhb_fracs:
        sp.params["fbhb"] = fb
        all_fluxs.append(sp.get_spectrum()[1]) #lsun per hz
    nmod = len(modpars)
    # flatten flux for writing
    flat_flux = np.array([all_fluxs[j][i]
                          for i in range(len(ages))
                          for j in range(len(fbhb_fracs))])
    # this function is flexible, ndim can be 3/4/n.
    writeASCII(fileout, lam, flat_flux, modpars,
               nmod=nmod, ndim=2, npar=2,
               par1='age', par2='fbhb')
    return
#---------------------------------------------------------------------
# ASCII FILE: WRITE AND COMPILE
#---------------------------------------------------------------------
# assumes you have $CLOUDY_EXE and $CLOUDY_DATA_PATH set as sys vars.

# name of ascii file
ascii_file = "FSPS_FBHB.ascii"

# or if there is an already-compiled one you want to use, specify here
compiled_ascii = "{}.mod".format(ascii_file.split(".")[0])

if exec_write_ascii:
    print("Executing write ascii sequence...")
    if not compiledExists(ascii_file):
        print("No compiled model exists...Writing.")
        fbhb_ascii(ascii_file)
        print("Compiling {} with Cloudy".format(ascii_file))
        compileASCII(ascii_file)
        print("Checking to see if compilation was successful...")
        if checkCompiled(ascii_file):
            print("Your model {} is ready to run.".format(compiled_ascii))
        else:
            sys.exit()
    else:
        print("{} already exists.".format(compiled_ascii))

#---------------------------------------------------------------------
# WRITE CLOUDY INPUT
#---------------------------------------------------------------------
# local folder to read and write *.in, *.out files
mod_dir = "/astro/users/ebyler/research/newem/output_fbhb/"
mod_prefix = "FAU"

# GRID PARAMETERS FOR CLOUDY RUN
#--------------
ages = np.array([3.0e9, 5.0e9, 10.0e9])
logUs =  np.array([-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0])
logZs = np.array([0.0])
#
Rinners =  np.array([19.])
nhs = np.array([100.0])
fbhb_fracs = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
#--------------
pars = np.array([(Z, a, U, R, calcForLogQ(logU=U, Rinner=10.0**R, nh=n), n, -1.0, fb)
                 for Z in logZs
                 for a in ages
                 for U in logUs 
                 for R in Rinners
                 for n in nhs
                 for fb in fbhb_fracs])

input_dict = {"r_in_pc":False,
              "set_name":"dopita",
              "dust":True,
              "re_z":False,
              "cloudy_mod":compiled_ascii,
              "extras":'save last physical conditions ".phys"',
              "extra_output":True,
              "to_file":True,
              "verbose":False,
              "par1":"age",
              "par2":"fbhb"
              }

print("{} models".format(len(pars)))
full_model_names = ["{}{}".format(mod_prefix, n+1)
                    for n in range(len(pars))]

printParFile(mod_dir, mod_prefix, pars)

if exec_write_input:
    print("Writing input files...")
    for par, name in zip(pars, full_model_names):
        cloudyInput(mod_dir, name,
                    logZ=par[0], age=par[1], logU=par[2], r_inner=par[3],
                    logQ=par[4], dens=par[5], efrac=par[6], par2val=par[7],
                    **input_dict)
    print("Wrote {} param files".format(len(pars)))
else:
    print("Skipping input writing.")
#---------------------------------------------------------------------
# RUN CLOUDY ON ALL INPUT FILES
#---------------------------------------------------------------------
if exec_run_cloudy:
    print("Running Cloudy....")
    runMake(dir_=mod_dir, n_proc=4, model_name=mod_prefix)
    print("Cloudy finished.")
else:
    print("Not running Cloudy. Skipping to formatting output.")
#---------------------------------------------------------------------
# FORMAT OUTPUT
#---------------------------------------------------------------------
if exec_write_output:
    print("Formatting output files...\n")
    formatAllOutput(mod_dir, mod_prefix)
else:
    print("\n\nNot formatting output. DONE.")
if exec_gen_FSPS_grid:
    print("Creating FSPS input grids...")
    writeFormattedOutput(mod_dir, mod_prefix, "_FBHB")



#-----------------------------------------------------------------------
# ocelote
#-----------------------------------------------------------------------
#set up outfile and essential info
outstr = 'fbhb'
jobfile = '/xdisk/senchp/cloudyfsps/fbhb/cloudy_{0}_jobs.cfg'.format(outstr)
jobfolder = '/xdisk/senchp/cloudyfsps/fbhb/output_{0}/'.format(outstr)

prefix_str = '''#!/bin/bash
#PBS -N cloudyfbhb
#PBS -W group_list=dpstark
#PBS -q oc_windfall
#PBS -l select=1:ncpus=12:mem=30gb
#PBS -M senchp@email.arizona.edu
#PBS -m bea
### #PBS -l walltime=05:00:00
### #PBS -l cput=30:00:00
source activate sci2
cd /xdisk/senchp/cloudyfsps/fbhb/
export CLOUDY_EXE='/home/u7/senchp/bin/cloudyrun'
export CLOUDY_DATA_PATH='/home/u7/senchp/builds/cloudy/c13.04/data/'
'''

#-----------------------------------------------------------------------
f = open(jobfile, 'w')
f.write(prefix_str+'\n')

# array job!
modstr = "python {ex} {dir} {prefix} $PBS_ARRAY_INDEX \n".format(
    ex='/home/u7/senchp/builds/cloudyfsps/scripts/runCloudy.py',
    dir=mod_dir, prefix=mod_prefix)
f.write(modstr+'\n')

f.close()

print('Added {0} jobs to {1}'.format(len(pars), jobfile.split('/')[-1]))
