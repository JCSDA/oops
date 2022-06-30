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
    #Read in required first file
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
    print(filepaths)
    for anfile, truthfile in zip(filepaths, truthfilepaths):
        # read fields
        field, date      = read_field(anfile)
        # read "truth" fields
        truthfield, date = read_field(truthfile)
        dates.append(date)
        # compute and save RMSE
        rmse_fields.append(np.linalg.norm(field - truthfield) / np.sqrt(len(field)))
    
    
    #Read in optional second analysis file    
    if not args.filepath2 is None:
        filepaths2 = []
        ids = []
        for id in args.times.split(","):
            ids.append(id)
            filepaths2.append(args.filepath2.replace("%id%", id))
            truthfilepaths.append(args.truthfilepath.replace("%id%", id))

        # RMSE vector
        rmse_fields2 = []
        #dates = []
        # Loop over filepaths
        print(filepaths2)
        for anfile, truthfile in zip(filepaths2, truthfilepaths):
            # read fields
            field2, _      = read_field(anfile)
            #read "truth" fields
            truthfield , _ =read_field(truthfile)
            #dates.append(date)
            # compute and save RMSE
            rmse_fields2.append(np.linalg.norm(field2 - truthfield) / np.sqrt(len(field2)))
            
    #Read in optional third analysis file    
    if not args.filepath3 is None:
        filepaths3 = []
        ids = []
        for id in args.times.split(","):
            ids.append(id)
            filepaths3.append(args.filepath3.replace("%id%", id))
            truthfilepaths.append(args.truthfilepath.replace("%id%", id))

        # RMSE vector
        rmse_fields3 = []
        #dates = []
        # Loop over filepaths
        print(filepaths3)
        for anfile, truthfile in zip(filepaths3, truthfilepaths):
            # read fields
            field3, _      = read_field(anfile)
            #read "truth" fields
            truthfield , _ =read_field(truthfile)
            #dates.append(date)
            # compute and save RMSE
            rmse_fields3.append(np.linalg.norm(field3 - truthfield) / np.sqrt(len(field3)))
    
    #Read in optional fourth analysis file    
    if not args.filepath4 is None:
        filepaths4 = []
        ids = []
        for id in args.times.split(","):
            ids.append(id)
            filepaths4.append(args.filepath4.replace("%id%", id))
            truthfilepaths.append(args.truthfilepath.replace("%id%", id))

        # RMSE vector
        rmse_fields4 = []
        #dates = []
        # Loop over filepaths
        print(filepaths4)
        for anfile, truthfile in zip(filepaths4, truthfilepaths):
            # read fields
            field4, _      = read_field(anfile)
            #read "truth" fields
            truthfield , _ =read_field(truthfile)
            #dates.append(date)
            # compute and save RMSE
            rmse_fields4.append(np.linalg.norm(field4 - truthfield) / np.sqrt(len(field4)))        
                
 

    # Plot RMSE
    fig, ax = plt.subplots()
    ax.set_title("RMSE time series")
    ax.set_xlabel("Time")
    ax.set_ylabel("Value")
    if not args.fileKey is None:
        Key=args.fileKey
    else:
        fieldstr=os.path.basename(args.filepath).split('.')
        Key=fieldstr[0]+'.'+fieldstr[1]+'(f1)'
    ax.plot(dates, rmse_fields, "r.-",label = Key)
    
    if not args.filepath2 is None:
        if not args.file2Key is None:
            Key=args.file2Key
        else:
            field2str=os.path.basename(args.filepath2).split('.')
            Key=field2str[0]+'.'+field2str[1]+'(f2)'
        ax.plot(dates, rmse_fields2, "b.-",label = Key)
   
    if not args.filepath3 is None:
        if not args.file3Key is None:
            Key=args.file3Key
        else:
            field3str=os.path.basename(args.filepath3).split('.')
            Key=field3str[0]+'.'+field3str[1]+'(f3)'
        ax.plot(dates, rmse_fields3, "g.-",label =  Key)
    
    if not args.filepath4 is None:
        if not args.file4Key is None:
            Key=args.file4Key
        else:
            field4str=os.path.basename(args.filepath4).split('.')
            Key=field4str[0]+'.'+field4str[1]+'(f4)'
        ax.plot(dates, rmse_fields4, "m.-",label = Key)        
    plt.xticks(rotation=10)
    ax.legend()

    # Save plot
    if args.output is None:
        plotpath = os.path.splitext(os.path.basename(args.filepath))[0]
    else:
        plotpath = args.output
    plotpath = plotpath.replace("%id%","") + ".jpg"
    plt.savefig(plotpath, format="jpg", dpi=300)
    plt.close()
    print(" -> plot produced: " + plotpath)

