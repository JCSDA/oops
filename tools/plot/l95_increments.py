#!/usr/bin/env python3
"""
@author: Mayeul Destouches
@author: modified by Benjamin Menetrier and Anna Shlyaeva for JEDI
@description: plot increment (analysis-background) and perfect increment (truth-background)
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
    # Load the background field
    bgfield, date = read_field(args.bgfilepath)

    # Variable index
    n = len(field)
    ind = list(range(1,n+1))
    if not args.recenter:
        indc=list(range(n))
        indlabels=ind
    else:    
        indc=list(range(round(n/2),n))+list(range(0,round(n/2)))
        indlabels=list(range(round(n/2)+1,n+1))+list(range(1,round(n/2)+1))
    # Plot increment (analysis - background) and optinally perfect increment (truth-background)
    fig, ax = plt.subplots()
    if not args.title is None:
        ax.set_title(args.title)
    else:    
        ax.set_title("Date: " + date)
    ax.set_xlabel("Variable index")
    ax.set_ylabel("Value")
    ax.plot(ind, np.subtract(field, bgfield)[indc], "r.-", label="increment")
    ax.set_xlim([1,n])
    plt.xticks(ind,indlabels)
    ax.xaxis.set_major_locator(plt.MaxNLocator(10))
    
    if not args.truthfilepath is None:
        # Load truth field
        truth, date = read_field(args.truthfilepath)
        ax.plot(ind, np.subtract(truth, bgfield)[indc], "b.--", label="perfect increment")

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
