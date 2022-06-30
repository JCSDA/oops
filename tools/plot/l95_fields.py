#!/usr/bin/env python3
"""
@author: Mayeul Destouches
@author: modified by Benjamin Menetrier for JEDI
@description: plot fields or increments if two paths are specified
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Read a field from L95 file and return the field and the date
def read_field(filename):
    # Check field file extension
    if not filename.endswith(".l95"):
        print("   Error: filepath extension should be .l95")
        sys.exit(1)
    # Load field
    f = open(filename, "r")
    lines = f.readlines()
    date = lines[1]
    field = []
    for item in lines[2].split():
       field.append(float(item))
    return field, date

# Read observation values from L95 obs files and return obs locations
# and obs values
def read_obs(filename):
    # Check obs file extension
    if not filename.endswith(".obt"):
        print("   Error: obs filepath extension should be .obt")
        sys.exit(1)
    # Load file
    f = open(filename, "r")
    lines = f.read().splitlines()
    # Header: ncolumns; column names; number of obs
    nheader = int(lines[0]) + 1
    # 3d column in the obs file specifies obs location
    obsloc_ind = 2
    # ObsValue's column is 3 + index of "ObsValue" in the header
    obsval_ind = lines.index("ObsValue") + 2
    # read obs locations and values
    obsloc = []
    obsval = []
    for line in lines[nheader+1:len(lines)]:
      items = line.split()
      obsloc.append(float(items[obsloc_ind]))
      obsval.append(float(items[obsval_ind]))
    return obsloc, obsval

def func(args):
    """! Plot fields"""

    # Load analysis field
    field, date = read_field(args.filepath)
    # Variable index
    n = len(field)
    ind = list(range(1,n+1))

    # Plot analysis field
    fig, ax = plt.subplots()
    ax.set_title("Date: " + date)
    ax.set_xlabel("Variable index")
    ax.set_ylabel("Value")
    ax.plot(ind, field, "r.-", label="analysis")
    ax.set_xlim([1,n])
    
    if not args.bgfilepath is None:
        # Load background field
        bgfield, date = read_field(args.bgfilepath)
        ax.plot(ind, bgfield, "b.--", label="background")

    if not args.truthfilepath is None:
        # Load truth field
        truthfield, date = read_field(args.truthfilepath)
        ax.plot(ind, truthfield, "k.-", label="truth")

    if not args.obsfilepath is None:
        # Load observations values
        obsloc, obsval = read_obs(args.obsfilepath)
        # obs locations are relative, multiply by the number of points:
        obsloc = [x*n for x in obsloc]
        ax.plot(obsloc, obsval, "ro", label="observations")

    ax.legend()
    # Save plot
    if args.output is None:
        plotpath = os.path.splitext(os.path.basename(args.filepath))[0]
    else:
        plotpath = args.output
    plotpath = plotpath + ".jpg"
    plt.savefig(plotpath, format="jpg", dpi=300)
    plt.close()
    print(" -> plot produced: " + plotpath)
