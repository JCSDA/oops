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
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

def func(args):
    """! Plot fields"""

    # Check file extension
    if not args.filepath.endswith(".l95"):
        print("   Error: filepath extension should be .l95")
        sys.exit(1)

    # Load field
    f = open(args.filepath, "r")
    lines = f.readlines()
    n = int(lines[0])
    date = lines[1]
    field = []
    for item in lines[2].split():
        field.append(float(item))

    if not args.basefilepath is None:
        # Check base file extension
        if not args.basefilepath.endswith(".l95"):
            print("   Error: filepath extension should be .l95")
            sys.exit(1)

        # Load base field
        fbase = open(args.basefilepath, "r")
        baselines = fbase.readlines()
        basefield = []
        for item in baselines[2].split():
            basefield.append(float(item))

        # Compute increment
        for i in range(len(field)):
            field[i] = field[i]-basefield[i]

    # Variable index
    ind = []
    ind = [i+1 for i in range(n)]

    # Plot cost function
    fig, ax = plt.subplots()
    ax.set_title("Date: " + date)
    ax.set_xlabel("Variable index")
    ax.set_ylabel("Value")
    ax.plot(ind, field, "r.-")
    ax.set_xlim([1,n])

    # Save plot
    if args.output is None:
        plotpath = os.path.splitext(os.path.basename(args.filepath))[0]
    else:
        plotpath = args.output
    if args.basefilepath is None:
        plotpath = plotpath + ".jpg"
    else:
        plotpath = plotpath + "_incr.jpg"
    plt.savefig(plotpath, format="jpg", dpi=300)
    plt.close()
    print(" -> plot produced: " + plotpath)
