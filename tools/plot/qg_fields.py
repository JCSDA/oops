#!/usr/bin/env python3
"""
@author: Mayeul Destouches
@author: modified by Benjamin Menetrier for JEDI
@description: plot fields or increments if two paths are specified
"""

import os
import sys
import subprocess
import argparse
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

def func(args):
    """! Plot fields"""

    # Geophysical parameters
    req = 6371229.0                      # Earth radius at equator (m)
    domain_zonal = 2 * np.pi * req       # Model domain in zonal direction (m)
    domain_meridional = 0.5 * np.pi *req # Model domain in meridional direction (m)
    plot_width=7 #units in inches
    plot_height=3.5 #units in inches
    vec_len_frac=20 # with --plotwind the longest vector plotted will be (1/vec_len_frac)*plot_width inches in length

    # Variables to plot
    variables = ["x", "q","u","v"]

    # File path
    filepaths = []
    ids = []
    if args.gif is None:
        filepaths.append(args.filepath)
        print(filepaths)
    else:
        check_convert = subprocess.getstatusoutput('convert --help')
        if check_convert[0] != 0:
            print("   Error: convert (imagemagick package) should be available to create an animated gif")
            sys.exit(1)
        if "%id%" in args.filepath:
           for id in args.gif.split(","):
               ids.append(id)
               filepaths.append(args.filepath.replace("%id%", id))
        else:
            print("   Error: filepath should include a %id% pattern for gif generation")
            sys.exit(1)

    # Loop over variables
    for variable in variables:
        # Fields vector
        fields = []
        if args.plotwind:
            fields_u = []
            fields_v = []
        if args.plotObsLocations:
            # Variables to get
            obs_types = ["Stream", "Wind", "WSpeed"]
            # Check file extension
            if not args.plotObsLocations.endswith(".nc"):
                print("   Error: filepath extension should be .nc")
                sys.exit(1)
            # Get data
            res = netCDF4.Dataset(args.plotObsLocations)
            locations = res.groups[list(res.groups.keys())[0]].groups["Location"].variables["values"][:,:]
            obs_lon = locations[:,0]
            obs_lat = locations[:,1]

        # Loop over filepaths
        for filepath in filepaths:
            # Check file extension
            if not filepath.endswith(".nc"):
                print("   Error: filepath extension should be .nc")
                sys.exit(1)

            # Load fields
            fields.append(netCDF4.Dataset(filepath).variables[variable][:])
            if args.plotwind:
                 fields_u.append(netCDF4.Dataset(filepath).variables["u"][:])
                 fields_v.append(netCDF4.Dataset(filepath).variables["v"][:])

        # Plotted fields vector
        fields_plot = []
        if args.plotwind:
            fields_u_plot = []
            fields_v_plot = []

        if args.basefilepath is None:
            for field in fields:
               fields_plot.append(field)
            if args.plotwind:
                for field_u in fields_u:
                    fields_u_plot.append(field_u)
                for field_v in fields_v:
                    fields_v_plot.append(field_v)
        else:
            # Check file extension
            if not args.basefilepath.endswith(".nc"):
                print("   Error: basefilepath extension should be .nc")
                sys.exit(1)

            # Load base fields
            field_base = netCDF4.Dataset(args.basefilepath).variables[variable][:]
            if args.plotwind:
                field_u_base = netCDF4.Dataset(args.basefilepath).variables["u"][:]
                field_v_base = netCDF4.Dataset(args.basefilepath).variables["v"][:]

            # Compute increments
            for field in fields:
                fields_plot.append(field-field_base)
            if args.plotwind:
                for field_u in fields_u:
                    fields_u_plot.append(field_u-field_u_base)
                for field_v in fields_v:
                    fields_v_plot.append(field_v-field_v_base)

        # Get geometry
        nz, ny, nx = fields_plot[0].shape
        levels = list(range(nz))
        z_coord = netCDF4.Dataset(filepaths[0]).variables["z"][:]
        lon_coord = netCDF4.Dataset(filepaths[0]).variables["lon"][:]
        lat_coord = netCDF4.Dataset(filepaths[0]).variables["lat"][:]

        min_lon = np.min(lon_coord)
        max_lon = np.max(lon_coord)
        min_lat = np.min(lat_coord)
        max_lat = np.max(lat_coord)

        # Get obs locations in other unit
        if args.plotObsLocations:
            indexes = []
            for ii in range(len(obs_lon)):
                obs_lat[ii] = 90/(max_lat-min_lat) * (obs_lat[ii] - min_lat)
                if (obs_lon[ii] < min_lon) or (obs_lon[ii] > max_lon) or (obs_lat[ii] < min_lat) or (obs_lat[ii] > max_lat):
                    indexes.append(ii)
            obs_lon = np.delete(obs_lon, indexes)
            obs_lat = np.delete(obs_lat, indexes)

        # Define color levels
        clevels = []
        vmax = 0.0
        for level in levels:
            for field in fields_plot:
                vmax = max(vmax, np.max(np.abs(field[level])))
        for level in levels:
            clevels.append(np.linspace(-vmax, vmax, 30))

        # Define plot
        params = {
            "font.size": 12,
            "text.latex.preamble" : r"\usepackage{amsmath}\usepackage{amsfonts}",
            "ytick.left": False,
            "ytick.labelleft": False,
        }
        plt.rcParams.update(params)
        my_formatter = mticker.FuncFormatter(lambda x, pos:"{:.0f}$\degree$E".format(x).replace("-", "\N{MINUS SIGN}"))
        if args.plotwind:
            # Select scale
            if args.basefilepath is None:
                scale = (vec_len_frac)*np.max(np.sqrt(np.square(fields_u_plot)+np.square(fields_v_plot)),axis=(0,1,2,3))/plot_width
            else:
                scale = (vec_len_frac)*np.max(np.sqrt(np.square(fields_u_plot)+np.square(fields_v_plot)),axis=(0,1,2,3))/plot_width
            dx_quiver = max(nx//20, 1)
            dy_quiver = max(ny//10, 1)

        # Initialize gif command
        if not args.gif is None:
             cmd = "convert -delay 20 -loop 0 "

        for iplot in range(0, len(filepaths)):
            fig, axs = plt.subplots(nrows=2, figsize=(plot_width, plot_height))

            # Loop over levels
            for level, ax in zip(levels, axs[::-1]):
                # Plot variable
                if args.basefilepath is None:
                    im = ax.contourf(lon_coord, lat_coord, fields_plot[iplot][level], cmap="plasma", levels=clevels[level])
                else:
                    im = ax.contourf(lon_coord, lat_coord, fields_plot[iplot][level], cmap="RdYlBu_r", levels=clevels[level])

                if args.plotwind:
                    # Plot wind field
                    ax.quiver(lon_coord[::dy_quiver, ::dx_quiver], lat_coord[::dy_quiver, ::dx_quiver],
                                  fields_u_plot[iplot][level, ::dy_quiver, ::dx_quiver], fields_v_plot[iplot][level, ::dy_quiver, ::dx_quiver],
                                  scale=scale, scale_units="inches")

                if args.plotObsLocations:
                    ax.scatter(obs_lon, obs_lat, marker='x', c='k', s=8, linewidths=0.5)


                # Set plot formatting
                ax.set_aspect("equal")
                cb = fig.colorbar(im, ax=ax, shrink=0.9, format=("%.1e" if variable == "q" else None))
                ax.set_ylabel("Altitude {:.0f}$\,$m".format(z_coord[level]))
                ax.xaxis.set_major_formatter(my_formatter)

                # Set title
                varname = dict(x="Streamfunction", q="Potential vorticity",u="Horizontal Wind", v="Vertical Wind").get(variable)
                unit = dict(x="m$^2$s$^{-1}$", q="s$^{-1}$", u="m/s",v="m/s").get(variable)
                if not args.title is None:
                    fig.suptitle(args.title + " - " + varname + " in " + unit)
                else:
                    fig.suptitle(varname + " in " + unit)
                fig.subplots_adjust(left=0.04, right=0.98, bottom=0.04, top=0.9, hspace=0.01)

            # Save plot
            if args.output is None:
                plotpath = os.path.splitext(os.path.basename(filepaths[iplot]))[0]
            else:
                plotpath = args.output
                if not args.gif is None:
                    plotpath = plotpath.replace("%id%", ids[iplot])
            if args.basefilepath is None:
                plotpath = plotpath + "_" + str(variable) + ".jpg"
            else:
                plotpath = plotpath + "_" + str(variable) + "_diff.jpg"
            if not args.gif is None:
                cmd = cmd + plotpath + " "
                if iplot == 0:
                    gifpath = plotpath.replace(".jpg", ".gif")
            plt.savefig(plotpath, format="jpg", dpi=300)
            plt.close()
            print(" -> plot produced: " + plotpath)

        if not args.gif is None:
            cmd = cmd + gifpath
            os.system(cmd)
            print(" -> gif produced: " + gifpath)
