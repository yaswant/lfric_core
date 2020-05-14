#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2019 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Need to set a non-interactive backend for suites
'''
Plot horizontal and vertical cross sections of desired fields with
fixed contour intervals
'''
import matplotlib
matplotlib.use('Agg')

import sys
import iris

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata
from read_data import read_ugrid_data

# Use magma colormap
from magma import magma_data
from matplotlib.colors import ListedColormap
magma = ListedColormap(magma_data, name='magma')
plt.register_cmap(name='magma', cmap=magma)

iris.FUTURE.netcdf_promote = True

# Size of regular grid
ny, nx = 360, 720


def make_figures(filein, plotpath, fields, vertical_spacing, lid, n_full,
                 figname):

    if vertical_spacing == 'um':
        # um L38 set
        zi_f = np.array([.0, .0005095, .0020380, .0045854, .0081519, .0127373,
                         .0183417, .0249651, .0326074, .0412688, .0509491,
                         .0616485, .0733668, .0861040, .0998603, .1146356,
                         .1304298, .1472430, .1650752, .1839264, .2037966,
                         .2246857, .2465938, .2695209, .2934670, .3184321,
                         .3444162, .3714396, .3998142, .4298913, .4620737,
                         .4968308, .5347160, .5763897, .6230643, .6772068,
                         .7443435, .8383348, 1.000000])*lid

    elif vertical_spacing == 'dcmip':
        # dcmip Stretched grid
        mu = 15.
        zi_f = np.zeros([n_full])
        for k in range(n_full):
            eta = float(k)/float(n_full)
            zi_f[k] = ((np.sqrt(mu*eta**2 + 1.) - 1.)
                      /(np.sqrt(mu + 1.) - 1) * lid)
    else:
        # assume uniform grid
        zmin = 0.0
        zmax = lid
        zi_f = np.linspace(zmin, zmax, n_full)

    zi_h = 0.5*(zi_f[1:] + zi_f[0:n_full-1])

    directions = ['xy', 'yz']

    for t in [-1]:
        interp_fig = plt.figure(figsize=(20, 10))
        if fields is None:
            fields = ['u1', 'u2', 'u3']
        for f, field in enumerate(fields):
            cube = read_ugrid_data(filein, field)
            levels_name = cube.dim_coords[0].name()

            # Set some levels for contours:
            levels = None
            if field == 'theta':
                levels = np.linspace(220, 330, 12)
            if field == 'u1':
                levels = np.linspace(-5, 45, 11)
            if field == 'u2':
                levels = np.linspace(-1.5, 1.5, 11)
            if field == 'u3':
                levels = np.linspace(-0.25, 0.25, 11)
            if field == 'exner':
                # exner will be converted to hPa
                levels = np.linspace(916, 1020, 14)
            if field == 'density':
                levels = np.arange(0, 1.4, 0.05)

            n_levs = len(cube.coord(levels_name).points)

            plot_data = np.zeros((ny, nx, n_levs))

            time = np.around(cube.coord('time').points, decimals=1)

            # Compute the horizontal grid
            x = np.around(cube.coord('longitude').points, decimals=5)
            y = np.around(cube.coord('latitude').points, decimals=5)

            xmin = np.amin(x)
            xmax = np.amax(x)
            ymin = np.amin(y)
            ymax = np.amax(y)

            # Generate a regular grid to interpolate the data.
            xi = np.linspace(xmin, xmax, nx)
            yi = np.linspace(ymin, ymax, ny)

            xf, yf = np.meshgrid(xi, yi)

            # Choose the correct vertical level set
            if n_full == n_levs:
                zi = zi_f
            else:
                zi = zi_h

            # Interpolate using delaunay triangularization
            for p, l in enumerate(range(n_levs)):
                data = cube.data[t, l]
                fi = griddata((x, y), data, (xf, yf), method='linear')
                fi_n = griddata((x, y), data, (xf, yf), method='nearest')
                fi[np.isnan(fi)] = fi_n[np.isnan(fi)]

                plot_data[:, :, l] = fi

                if field == 'exner_pressure':
                    # Convert to hPa
                    rd = 287.05
                    p0 = 100000.0
                    kappa = rd/1005.0
                    plot_data[:, :, l] = 0.01*fi**(1.0/kappa)*p0

            for d, direction in enumerate(directions):

                nxplots = len(fields)
                nyplots = len(directions)
                nplots = nxplots*nyplots

                plot_long = int(nx/2)
                plot_lat = int(ny/2)
                plot_level = n_levs-2

                plotnum = d*nxplots + f + 1
                ax = interp_fig.add_subplot(nyplots, nxplots, plotnum)

                if direction == 'xz':
                    x1, x2 = np.meshgrid(xi, zi)
                    x3 = plot_data[plot_lat, :, :].T
                    plt.title([field, ', lat = ',
                               yi[plot_lat]*360./np.real(nx)])
                if direction == 'yz':
                    x1, x2 = np.meshgrid(yi, zi)
                    x3 = plot_data[:, plot_long, :].T
                    plt.title([field, ', long = ',
                               xi[plot_long]*360./np.real(nx)])
                if direction == 'xy':
                    x2, x1 = np.meshgrid(yi, xi)
                    x3 = plot_data[:, :, plot_level].T
                    plt.title([field, ', Height = ', zi[plot_level]])

                CS = plt.contourf(x1, x2, x3, levels=levels, cmap=magma)
                plt.colorbar(cmap=magma)
                CL = plt.contour(x1, x2, x3,
                                 levels=levels, linewidths=0.5, colors='k')

        plt.tight_layout()

        pngfile = '%s/%s-winds-time%s.png' % (plotpath, figname, time[t])
        plt.savefig(pngfile)
        plt.close()


if __name__ == "__main__":

    try:
        args = sys.argv[:]
        files, plotpath, vertical_grid, lid, n_full, figname = args[1:7]
        field_list = None
        if len(args[:]) > 7:
            field_list = args[7].split(':')
    except ValueError:
        print("Usage: {0} <filein> <plotpath> <vertical_grid> <lid>"
              " <n_full> <figname> [<fields_list>]"
              .format(sys.argv[0]))
        exit(1)

    file_list = files.split(':')
    for filein in file_list:
        make_figures(filein, plotpath, field_list, vertical_grid,
                     int(lid), int(n_full), figname)
