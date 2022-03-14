#!/usr/bin/env python3
"""
@author: Mayeul Destouches
@author: modified by Anna Shlyaeva for JEDI
@description: plot time series of RMSE
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
    return np.asarray(field), date


def func(args):
    """! Plot time series of RMSEs"""

    # File path
    filepaths = []
    truthfilepaths = []
    ids = []
    for id in args.times.split(","):
        ids.append(id)
        filepaths.append(args.filepath.replace("%id%", id))
        truthfilepaths.append(args.truthfilepath.replace("%id%", id))

    # RMSE vector
    rmse_fields = []
    dates = []
    # Loop over filepaths
    for anfile, truthfile in zip(filepaths, truthfilepaths):
        # read fields
        field, date      = read_field(anfile)
        # read "truth" fields
        truthfield, date = read_field(truthfile)
        dates.append(date)
        # compute and save RMSE
        rmse_fields.append(np.linalg.norm(field - truthfield) / np.sqrt(len(field)))

    # Plot RMSE
    fig, ax = plt.subplots()
    ax.set_title("RMSE time series")
    ax.set_xlabel("Time")
    ax.set_ylabel("Value")
    ax.plot(dates, rmse_fields, "r.-")
    plt.xticks(rotation=10)

    # Save plot
    if args.output is None:
        plotpath = os.path.splitext(os.path.basename(args.filepath))[0]
    else:
        plotpath = args.output
    plotpath = plotpath.replace("%id%","") + ".jpg"
    plt.savefig(plotpath, format="jpg", dpi=300)
    plt.close()
    print(" -> plot produced: " + plotpath)

