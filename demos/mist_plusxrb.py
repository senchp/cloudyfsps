#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import, unicode_literals)
import os
import sys
import itertools

import numpy as np
from astropy.analytic_functions import blackbody_nu, blackbody_lambda
import astropy.units as u
import astropy.constants as con

import fsps
from cloudyfsps.ASCIItools import (writeASCII, compileASCII, checkCompiled, compiledExists)
from cloudyfsps.cloudyInputTools import writeParamFiles
from cloudyfsps.generalTools import calcForLogQ

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
# 4. Writes input file for ocelote grid job to run 
#    cloudy and format output.

zsun = 0.0142

exec_write_ascii = True
exec_write_input = True

# Function to write the ascii file.
# This is where you set the properties of the
# ionizing spectrum (SSP/CSFH, IMF, FBHB, etc)

def mist_ascii(fileout, **kwargs):
    # change these parameters to modify the ionizing source grid
    # default mode is to produce an ascii grid in age and Z,
    # though different variables and more dimensions are possible.
    sp_dict = dict(zcontinuous=1,   # ssp interpolated to logzsol val; 
                                    # zmet is ignored
                   imf_type=2,      # kroupa
                   sfh=0,           # 0: ssp
                   const=0.0,
                   sf_start=0.0)
    sp = fsps.StellarPopulation(**sp_dict)
    # all ages and Zs
    ages = 10.**sp.log_age
    logZs = np.log10(sp.zlegend/zsun)
    # xrb_bbnorms = np.arange(0.0, 100.0, 5.0)
    xrb_bbnorms = [0., 1., 5., 10., 20., 30., 40., 50, 75., 100.]
    # fraction of stellar luminosity in hot ~keV blackbody

    modpars = [(age, logZ, xrbn) for age in ages for logZ in logZs
            for xrbn in xrb_bbnorms]

    lam = sp.wavelengths
    all_fluxs = []
    for logZ in logZs:
        sp.params['logzsol'] = logZ
        all_fluxs.append(sp.get_spectrum()[1]) #lsun per hz

    xrb_fluxs = []
    for xrbn in xrb_bbnorms:
        temp = (0.5*u.keV / con.k_B).to(u.K)
        # conversion factor to bump total lum to 1 Lsun
        conv = (1.0*con.L_sun / (con.sigma_sb * temp**4)).to(u.cm**2)

        f = blackbody_nu(lam*u.angstrom, temp)
        f *= np.pi*u.steradian * conv * xrbn

        xrb_fluxs.append((f/con.L_sun).to(1/u.Hz).value)

    # w = np.where( (lam < 228.) )
    # print('fsps fluxes', all_fluxs[0][0][w])
    # print(all_fluxs[0][0][w].shape)
    # print('max', np.max(all_fluxs[0][w]))

    # print('xrb fluxes', xrb_fluxs[0][w])
    # print(xrb_fluxs[0][w].shape)
    # print('max', np.max(xrb_fluxs[0][w]))
    # print('xrb fluxes', xrb_fluxs[4][w])
    # print('max', np.max(xrb_fluxs[4][w]))
    # stop

    nmod = len(modpars)
    # flatten flux for writing
    flat_flux = np.array([all_fluxs[j][i] + xrb_fluxs[k]
                          for i in range(len(ages))
                          for j in range(len(logZs))
                          for k in range(len(xrb_bbnorms))])
    # this function is flexible, ndim can be 3/4/n.
    writeASCII(fileout, lam, flat_flux, modpars,
               nx=len(lam), ndim=3, npar=3, nmod=nmod,
               par1='age', par2='logZ', par3='xrbn')
    return
#---------------------------------------------------------------------
# ASCII FILE: WRITE AND COMPILE
#---------------------------------------------------------------------
# assumes you have $CLOUDY_EXE and $CLOUDY_DATA_PATH set as sys vars.

# name of ascii file
ascii_file = 'FSPS_MIST_SSP_XRB.ascii'

# or if there is an already-compiled one you want to use, specify here
compiled_ascii = '{}.mod'.format(ascii_file.split('.')[0])
if exec_write_ascii:
    print("Executing write ascii sequence...")
    if not compiledExists(ascii_file):
        print("No compiled model exists...Writing.")
        print(ascii_file)
        mist_ascii(ascii_file)
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
mod_dir = '/xdisk/senchp/cloudyfsps/mist_xrb/output_mist_ssp_xrb/'
mod_prefix = 'ZAUX'

# GRID PARAMETERS FOR CLOUDY RUN
#--------------
# ages = np.array([0.5e6, 1.0e6, 2.0e6, 3.0e6, 5.0e6, 7.0e6, 10.0e6])
ages = np.array([0.5e6, 1.0e6, 2.0e6])
# logUs =  np.array([-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0])
logUs =  np.array([-4.0, -3.0, -2.0])
# logZs =  np.array([-2.5, -1.5, -0.75, -0.50, -0.25, 0.0, 0.25, 0.5])
logZs =  np.array([-2.5, -1.5, -0.50])

# xrb_bbnorms = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0])
xrb_bbnorms = np.array([0.0, 1.0, 10.0, 50.0, 100.0])

#
Rinners =  np.array([19.])
nhs = np.array([100.0])
#--------------
pars = np.array([(Z, a, U, R, calcForLogQ(logU=U, Rinner=10.0**R, nh=n), n, -1.0, 0.0, xn)
                 for Z in logZs
                 for a in ages
                 for U in logUs 
                 for R in Rinners
                 for n in nhs
                 for xn in xrb_bbnorms])

if exec_write_input:
    print('Writing input files...')
    writeParamFiles(dir_=mod_dir,
                    model_prefix=mod_prefix,
                    cloudy_mod=compiled_ascii, # file from above
                    run_cloudy=False, #don't run yet
                    ages=ages,
                    logZs=logZs,
                    logUs=logUs,
                    r_inners=Rinners,
                    nhs=nhs,
                    use_Q=True,
                    verbose=False, #don't print output to screen
                    set_name='dopita',
                    dust=False,
                    ntabpars=3,
                    tabpar1='age',
                    tabpar1val=ages,
                    tabpar2='logZ',
                    tabpar2val=logZs,
                    tabpar3='xrbn',
                    tabpar3val=xrb_bbnorms,
                    extra_output=True,
                    extras='set WeakHeatCool 0.001\nsave last cooling each ".cool"')
    print('Wrote {} param files'.format(len(pars)))
else:
    print('Skipping input writing.')


#=======================================================================
#-----------------------------------------------------------------------
#=======================================================================

#-----------------------------------------------------------------------
#print all the jobs you would like to run into myjobs.cfg
#-----------------------------------------------------------------------
#set up outfile and essential info
outstr = 'mist_ssp_xrb'
jobfile = '/xdisk/senchp/cloudyfsps/mist_xrb/cloudy_{0}_jobs.cfg'.format(outstr)
jobfolder = '/xdisk/senchp/cloudyfsps/mist_xrb/output_{0}/'.format(outstr)

prefix_str = '''#!/bin/bash
#PBS -N cloudymist
#PBS -W group_list=dpstark
#PBS -q oc_windfall
#PBS -l select=1:ncpus=12:mem=30gb
#PBS -M senchp@email.arizona.edu
#PBS -m bea
### #PBS -l walltime=05:00:00
### #PBS -l cput=30:00:00
source activate sci2
cd /xdisk/senchp/cloudyfsps/mist_xrb/
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
