#!/usr/bin/env python3
"""
@author: Mayeul Destouches
@author: modified by Benjamin Menetrier and Anna Shlyaeva for JEDI
@description: plot analysis and background (if specified) errors compared to truth
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

def func(args):
    """! Plot errors"""

    # Load analysis field
    field, date = read_field(args.filepath)
    # Load truth field
    truth, date = read_field(args.truthfilepath)

    # Variable index
    n = len(field)
    ind = list(range(1,n+1))

    # Plot analysis error field
    fig, ax = plt.subplots()
    ax.set_title("Date: " + date)
    ax.set_xlabel("Variable index")
    ax.set_ylabel("Value")
    ax.plot(ind, np.subtract(field, truth), "r.-", label="analysis error")
    ax.set_xlim([1,n])

    if not args.bgfilepath is None:
        # Load background field
        bgfield, date = read_field(args.bgfilepath)
        ax.plot(ind, np.subtract(bgfield, truth), "b.--", label="background error")

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
