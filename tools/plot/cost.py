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

    # Read log file to find observation types and count outer iterations (nouter)
    # and total inner iterations (ninner)
    obs_types = []
    nouter = 0
    ninner = 0
    pattern_test = "Test     :"
    with open(args.filepath, "r") as file:
        for line in file:
            # Get observation types
            pattern = "Nonlinear Jo\("
            if re.search(pattern, line) and not re.search(pattern_test, line):
                obs_type = line.split('(', 1)[1].split(')', 1)[0]
                if not obs_type in obs_types:
                    obs_types.append(obs_type)

            # Get number of outer iterations (+1 with initial value)
            pattern = "Nonlinear J ="
            if re.search(pattern, line) and not re.search(pattern_test, line):
                nouter += 1

            # Get total number of inner iterations
            pattern = "Quadratic cost function: J "
            if re.search(pattern, line) and not re.search(pattern_test, line):
                ninner += 1

    print(nouter-1, "outer iterations with a total of", ninner, "inner iterations")

    # Initialization
    Jb = np.zeros((nouter))
    Jo = np.zeros((nouter, len(obs_types)))
    J = np.zeros((nouter))
    indexJ = []
    quadJ = np.zeros((ninner))
    indexQuadJ = np.zeros((ninner))

    # Read log file to get cost function
    with open(args.filepath, "r") as file:
        iter = 0
        for line in file:
            # Get nonlinear Jb
            pattern = "CostJb   : Nonlinear Jb ="
            if re.search(pattern, line) and not re.search(pattern_test, line):
                Jb[iter] = float(line.split("=")[1])

            # Get nonlinear Jo
            for iobs in range(0, len(obs_types)):
               pattern = "CostJo   : Nonlinear Jo\(" + obs_types[iobs] + "\) ="
               if re.search(pattern, line) and not re.search(pattern_test, line):
                  Jo[iter, iobs] = float(line.split("=", 1)[1].split(",", 1)[0])

            # Get nonlinear J
            pattern = "CostFunction: Nonlinear J ="
            if re.search(pattern, line) and not re.search(pattern_test, line):
                J[iter] = float(line.split("=")[1])
                iter += 1


    with open(args.filepath, "r") as file:
        iter = 0
        for line in file:
            # Get quadratic J
            pattern = "Quadratic cost function: J "
            if re.search(pattern, line) and not re.search(pattern_test, line):
                indexQuadJ[iter] = int(re.split("[()]", line)[1])
                if indexQuadJ[iter] == 1:
                    indexJ.append(iter)
                quadJ[iter] = float(line.split("=")[1])
                iter += 1
        indexJ.append(iter)

    # Post-processing
    iters = np.linspace(1, ninner, ninner, dtype=np.int8)
    Jo_tot = J-Jb

    # Plot cost function
    fig, ax = plt.subplots()
    ax.set_title("Nonlinear and quadratic cost functions")
    ax.set_xlabel("Inner iteration")
    ax.set_ylabel("Cost function")
    plt.xticks(iters, iters)
    ax.plot(indexJ, Jb, "rv", label='Jb')
    ax.plot(indexJ, Jo_tot, "bo", label='Jo', fillstyle="none")
    ax.plot(indexJ, J, "kx", label='J', ms=7)
    ax.plot(iters, quadJ, "k-", label='quad J')
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
