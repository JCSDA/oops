#!/usr/bin/env python3
"""
@author: Benjamin Menetrier
@description: plot the nonlinear cost function at each outer iteration
"""

import os
import sys
import numpy as np
import re
import matplotlib.pyplot as plt

def subfunc(args):
    # Observation types

    # Read log file to find observation types and count outer iterations
    obs_types = []
    niter = 0
    with open(args.filepath, "r") as file:
        for line in file:
            # Get observation types
            pattern = "Nonlinear Jo"
            if re.search(pattern, line):
               obs_type = line.split('(', 1)[1].split(')', 1)[0]
               if not obs_type in obs_types:
                  obs_types.append(obs_type)

            # Get nonlinear J
            pattern = "Nonlinear J ="
            if re.search(pattern, line):
                niter += 1

    # Initialization
    Jb = np.zeros((niter))
    Jo = np.zeros((niter, len(obs_types)))
    J = np.zeros((niter))

    # Read log file to get cost function
    with open(args.filepath, "r") as file:
        iter = 0
        for line in file:
            # Get nonlinear Jb
            pattern = "Nonlinear Jb ="
            if re.search(pattern, line):
                Jb[iter] = float(line.split()[5])

            # Get nonlinear Jo
            for iobs in range(0, len(obs_types)):
               pattern = "Nonlinear Jo\(" + obs_types[iobs] + "\) ="
               if re.search(pattern, line):
                  Jo[iter, iobs] = float(line.split("=", 1)[1].split(",", 1)[0])

            # Get nonlinear J
            pattern = "Nonlinear J ="
            if re.search(pattern, line):
                J[iter] = float(line.split()[4])
                iter += 1

    # Post-processing
    iters = np.linspace(1, niter, niter, dtype=np.int8)
    Jo_tot = J-Jb

    # Plot cost function
    fig, ax = plt.subplots()
    ax.set_title("Nonlinear cost function")
    ax.set_xlabel("Outer iteration")
    ax.set_ylabel("Jb (red), Jo (blue) and J (black)")
    plt.xticks(iters, iters)
    ax.plot(iters, Jb, "r-", iters, Jo_tot, "b-", iters, J, "k-")

    # Save plot
    if args.output is None:
        plotpath = os.path.splitext(os.path.basename(args.filepath))[0]
    else:
        plotpath = args.output
    plotpath = plotpath + ".jpg"
    plt.savefig(plotpath, format="jpg", dpi=300)
    plt.close()
    print(" -> plot produced: " + plotpath)
