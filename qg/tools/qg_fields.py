#!/usr/bin/env python3
"""
@author: Mayeul Destouches
@author: modified by Benjamin Menetrier for JEDI
@description: plot fields or increments if two paths are specified
"""

import os
import argparse
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("anafilepath", help="Analysis file path")
parser.add_argument("bkgfilepath", nargs='?', help="Background file path", default=None)
parser.add_argument('--plotwind', dest='plotwind', action='store_true', help="Plot wind")
args = parser.parse_args()

# Geophysical parameters
req = 6371229.0 # Earth radius at equator (m)
domain_zonal = 2 * np.pi * req # Model domain in zonal direction (m)
domain_meridional = 0.5 * np.pi *req # Model domain in meridional direction (m)

# Variables to plot
variables = ['x', 'q']

# Load NetCDF files
anares = netCDF4.Dataset(args.anafilepath)
if args.bkgfilepath != None:
    bkgres = netCDF4.Dataset(args.bkgfilepath)

# Loop over variables in analysis file
for variable in anares.variables:
    if variable in variables:
        if args.bkgfilepath == None:
            # Load field
            field = anares.variables[variable][:]
        else:
            # Load analysis and background
            anafield = anares.variables[variable][:]
            bkgfield = bkgres.variables[variable][:]

            # Compute increment
            field = anafield-bkgfield

        # Get geometry
        nz, ny, nx = field.shape
        levels = list(range(nz))
        z_coord=anares.variables['z'][:]
        dx = domain_zonal / nx # zonal grid cell in km
        dy = domain_meridional / (ny + 1) # meridional cell in km
        x_coord = (np.arange(1, nx+1) - 0.5) * dx
        lon_coord = x_coord / domain_zonal * 360 - 180
        y_coord = np.arange(1, ny+1) * dy
        lat_coord = y_coord / domain_meridional * 90
        xx, yy = np.meshgrid(lon_coord, lat_coord)
   
        # Load wind field
        if args.plotwind:
            dx_quiver = max(nx//20, 1)
            dy_quiver = max(ny//10, 1)
            if args.bkgfilepath == None:
               # Load field
               u = anares.variables['u'][:]
               v = anares.variables['v'][:]
            else:
               # Load analysis and background
               anau = anares.variables['u'][:]
               bkgu = bkgres.variables['u'][:]
               anav = anares.variables['v'][:]
               bkgv = bkgres.variables['v'][:]

               # Compute increment
               u = anau-bkgu
               v = anav-bkgv

        # Define color levels
        mini = np.min(field, axis=(1,2))
        maxi = [np.max(field[level]) for level in levels]
        bounds = [max(np.abs(mn), np.abs(mx)) for mn, mx in zip(mini, maxi)]
        clevels = [np.linspace(-bound, bound, num=30) for bound in bounds]
        
        # Define plot
        params = {
            'font.size': 12,
            'text.latex.preamble' : r'\usepackage{amsmath}\usepackage{amsfonts}',
            'ytick.left': False,
            'ytick.labelleft': False,
        }
        plt.rcParams.update(params)
        my_formatter = mticker.FuncFormatter(lambda x, pos: 
                        '{:.0f}$\degree$E'.format(x).replace('-', '\N{MINUS SIGN}'))

        fig, axs = plt.subplots(nrows=2, figsize=(7,3.5))

        # Loop over levels
        for level, ax in zip(levels, axs[::-1]):
            # Plot variable
            im = ax.contourf(xx, yy, field[level], cmap='plasma', levels=clevels[level])

            if args.plotwind:
                # Select scale
                if args.bkgfilepath == None:
                   scale = 500
                else:
                   scale = 100

                # Plot wind field
                ax.quiver(xx[::dy_quiver, ::dx_quiver], yy[::dy_quiver, ::dx_quiver], 
                          u[level, ::dy_quiver, ::dx_quiver], v[level, ::dy_quiver, ::dx_quiver], 
                          scale=scale, scale_units='inches')

            # Set plot formatting
            ax.set_aspect('equal')
            cb = fig.colorbar(im, ax=ax, shrink=0.9, 
                              format=('%.1e' if variable == 'q' else None))
            ax.set_ylabel('Altitude {:.0f}$\,$m'.format(z_coord[level]))
            ax.xaxis.set_major_formatter(my_formatter)

        # Set title
        varname = dict(x='Streamfunction', q='Potential vorticity').get(variable)
        unit= dict(x='m$^2$s$^{-1}$', q='s$^{-1}$').get(variable)
        fig.suptitle(varname + ' in ' + unit)
        fig.subplots_adjust(left=0.04, right=0.98, bottom=0.01, top=0.9, hspace=0.01)
    
        # Save plot
        if args.bkgfilepath == None:
           plotpath = os.path.splitext(args.anafilepath)[0] + "_" + str(variable) + ".jpg"
        else:
           plotpath = os.path.splitext(args.anafilepath)[0] + "_incr_" + str(variable) + ".jpg"
        plt.savefig(plotpath, format="jpg", dpi=300)
        plt.close()
        print("Plot produced: " + plotpath)
