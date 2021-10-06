#!/usr/bin/env python3
"""
@author: Benjamin Menetrier
@description: plot the nonlinear cost function at each outer iteration
"""

import os
import argparse
import numpy as np
import re
import matplotlib.pyplot as plt

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("filepath", help="File path")
args = parser.parse_args()

# Observation types
obs_types = ['Stream', 'Wind' , 'WSpeed']

# Read log file to count outer iterations
niter = 0
with open(args.filepath, "r") as file:
    for line in file:
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
              Jo[iter, iobs] = float(line.split()[5].replace(',', ''))

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
ax.plot(iters, Jb, 'r-', iters, Jo_tot, 'b-', iters, J, 'k-')

# Save plot
plotpath = os.path.splitext(args.filepath)[0] + ".jpg"
plt.savefig(plotpath, format="jpg", dpi=300)
plt.close()
print("Plot produced: " + plotpath)
