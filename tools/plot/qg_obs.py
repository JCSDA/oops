#!/usr/bin/env python3
"""
@author: Benjamin Menetrier
@description: print the observation type, location, values and hofx
"""

import os
import sys
import argparse
import numpy as np
import netCDF4
import datetime
import pathlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

def func(args):
    """! Plot obs values"""

    # Variables to plot
    obs_types = ["Stream", "Wind", "WSpeed"]

    # Check file extension
    if not args.filepath.endswith(".nc"):
        print("   Error: filepath extension should be .nc")
        sys.exit(1)

    # Open output file
    if args.output is None:
        txtpath = os.path.splitext(os.path.basename(args.filepath))[0]
    else:
        txtpath = args.output
    txtpath = txtpath + ".txt"
    f = open(txtpath, "w")
    f.write("# location / value / hofx\n")

    # Get data
    res = netCDF4.Dataset(args.filepath)
    for obs_type in obs_types:
        if obs_type in res.groups:
           nobs = res.groups[obs_type].dimensions["nobs"].size
           location = res.groups[obs_type].groups["Location"].variables["values"][:,:]
           value = res.groups[obs_type].groups["ObsValue"].variables["values"][:,:]
           hofx = res.groups[obs_type].groups["hofx"].variables["values"][:,:]
           for iobs in range(0, nobs):
               f.write(str(location[iobs,:]) + " / " + str(value[iobs,:]) + " / " + str(hofx[iobs,:]) + "\n")
    f.close()
    print(" -> Observations values written in " + txtpath)
