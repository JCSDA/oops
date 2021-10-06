#!/usr/bin/env python3
"""
@author: Benjamin Menetrier
@description: print the observation type, location, values and hofx
"""

import os
import argparse
import numpy as np
import netCDF4
import datetime
import pathlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# Parser
parser = argparse.ArgumentParser()
parser.add_argument('filepath', help='File path')
args = parser.parse_args()

# Variables to plot
obs_types = ['Stream', 'Wind', 'WSpeed']

# Plot field
res = netCDF4.Dataset(args.filepath)

for obs_type in obs_types:
    if obs_type in res.groups:
       print('Observation type: ' + obs_type)
       nobs = res.groups[obs_type].dimensions['nobs'].size
       location = res.groups[obs_type].groups['Location'].variables['values'][:,:]
       value = res.groups[obs_type].groups['ObsValue'].variables['values'][:,:]
       hofx = res.groups[obs_type].groups['hofx'].variables['values'][:,:]
       for iobs in range(0, nobs):
           print("   Obs. " + str(iobs+1) + " location / value / hofx: " + str(location[iobs,:]) + " / " + str(value[iobs,:]) + " / " + str(hofx[iobs,:]))
