#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
file name: plot
@author: Sacha Gavino
last update: Jan 22
language: PYTHON 3
__________________________________________________________________________________________
short description:  plotting of the disk thermal model
__________________________________________________________________________________________
"""
import glob, sys
from pathlib import Path
import os
import numpy as np
import pandas as pd
from scipy.interpolate import griddata

import re
from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize, LogNorm
from astromugs.constants.constants import autocm


def density2D_grid(path='thermal/', vmin=1e-30, vmax=1e-15, cmap='gnuplot2', dens_type='mass',
                    xlim=None, ylim=None, dust=None, figsize=(10, 14)):
    """Plot all dust species on a single figure with subplots, plus total density."""
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[grid.columns[0]].values

    dens = pd.read_table(path + 'dust_density.inp', engine='python', header=None, skiprows=3)
    dens = dens[0].values
    nspecies = int(len(dens) / (nr * nt))
    dens = np.reshape(dens, (nspecies, nt, nr))

    grid = np.array(grid, copy=True)

    r_edge = grid[:nr+1] / autocm
    theta_edge = grid[nr+1:nr+1+nt+1]
    theta_edge[-1] = np.pi
    rr_edge, tt_edge = np.meshgrid(r_edge, theta_edge)
    R = rr_edge * np.sin(tt_edge)
    Z = rr_edge * np.cos(tt_edge)
    dens = np.array(dens,copy=True)
    dens[dens <= 1e-100] = 1e-100

    # Try to read grain sizes for subplot labels
    import os
    sizes_file = path + 'dust_sizes.inp'
    if os.path.isfile(sizes_file):
        sizes = np.loadtxt(sizes_file)
        sizes = np.atleast_1d(sizes)
    elif dust != None:
        rho_m = dust.rho_m #g.cm3
        sizes = dust.sizes()[0] # microns
        grain_mass = dust.grainmass() # in gram

    else:
        sizes = None

    # Layout: enough panels for all species + 1 total
    npanels = nspecies + 1
    ncols = min(nspecies, 4)
    nrows = int(np.ceil(npanels / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])

    for idx in range(nrows * ncols):
        ax = axes.flat[idx]
        if idx < nspecies:
            if dens_type == 'number':
                im = ax.pcolormesh(R, Z, dens[idx]/grain_mass[idx], cmap=cmap, shading='auto',
                                   norm=LogNorm(vmin=vmin, vmax=vmax))
                fig.colorbar(im, cax=cbar_ax, label=r'n$_\mathrm{d}$ [cm$^{-3}$]')
            if dens_type == 'surface':
                im = ax.pcolormesh(R, Z, 4*np.pi*sizes[idx]*1e-4*dens[idx]/grain_mass[idx], cmap=cmap, shading='auto',
                                   norm=LogNorm(vmin=vmin, vmax=vmax))
                fig.colorbar(im, cax=cbar_ax, label=r'surfaces [cm$^{-1}$]')
            elif dens_type == 'mass':
                im = ax.pcolormesh(R, Z, dens[idx], cmap=cmap, shading='auto',
                                   norm=LogNorm(vmin=vmin, vmax=vmax))
                fig.colorbar(im, cax=cbar_ax, label=r'$\rho_\mathrm{d}$ [g cm$^{-3}$]')
            ax.set_title(f'bin {idx+1}', fontsize=12)
            if sizes is not None and idx < len(sizes):
                s = sizes[idx]
                if s >= 1e3:
                    size_label = f'{s/1e3:.1f} mm'
                else:
                    size_label = f'{s:.2f} ' + r'$\mu$m'
                ax.text(0.05, 0.95, size_label, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top',
                        horizontalalignment='left', bbox=props)
        elif idx == nspecies:
            total = dens.sum(axis=0)
            im = ax.pcolormesh(R, Z, total, cmap=cmap, shading='auto',
                               norm=LogNorm(vmin=vmin, vmax=vmax))
            ax.set_title('total', fontsize=12)
        else:
            ax.set_visible(False)
            continue

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

    # Shared labels
    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel('r [au]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel('z [au]', fontsize=14)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    
    #fig.colorbar(im, cax=cbar_ax, label=r'$\rho_\mathrm{d}$ [g cm$^{-3}$]')

    plt.show()


def density1D_midplane(path='thermal/', vmin=1e-30, vmax=1e-15, dens_type='mass',
                        xlim=None, dust=None, figsize=(12, 8)):
    """
    Plots the 1D dust density profile in the midplane (z=0 / theta=pi/2) 
    as a function of radius for each dust species.

    Parameters:
    -----------
    path : str
        Path to the directory containing RADMC-3D files.
    vmin, vmax : float
        Limits for the Y-axis (density).
    dens_type : str
        Type of density to plot: 'mass' (g/cm^3), 'number' (cm^-3), or 'surface' (cm^-1).
    xlim : tuple/list, optional
        Limits for the X-axis (radius in au).
    dust : object, optional
        An external dust object containing grain sizes and masses if files are missing.
    figsize : tuple
        Size of the output matplotlib figure.
    """
    
    # 1. Read grid structure and dust density data
    # Read the AMR grid file to extract dimensions (nr = radial bins, nt = theta bins)
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = np.array(grid[grid.columns[0]].values, copy=True)
    
    # Read the raw dust density file (flat 1D array of values)
    dens = pd.read_table(path + 'dust_density.inp', engine='python', header=None, skiprows=3)
    dens = dens[0].values
    
    # Deduce the number of dust species and reshape into a 3D array: (species, theta, radius)
    nspecies = int(len(dens) / (nr * nt))
    dens = np.reshape(dens, (nspecies, nt, nr))

    # 2. Extract radial coordinates at cell centers (convert from cm to au)
    # autocm is assumed to be a globally defined constant (1 au = 1.496e13 cm)
    r_edge = grid[:nr+1] / autocm
    r_center = 0.5 * (r_edge[:-1] + r_edge[1:])

    # 3. Identify the midplane index (theta = pi/2)
    # In RADMC-3D spherical coordinates, the equator is exactly at the midpoint of the theta axis
    idx_midplane = nt // 2 

    # 4. Read grain sizes and masses for plotting labels and conversions
    sizes_file = path + 'dust_sizes.inp'
    if os.path.isfile(sizes_file):
        sizes = np.loadtxt(sizes_file)
        sizes = np.atleast_1d(sizes)
    elif dust != None:
        rho_m = dust.rho_m #g.cm3
        sizes = dust.sizes()[0] # microns
        grain_mass = dust.grainmass() # in gram

    else:
        sizes = None

    # 5. Configure the figure layout (Grid of subplots)
    npanels = nspecies + 1  # Number of species + 1 extra panel for the total sum
    ncols = min(npanels, 3) # Maximum of 3 columns
    nrows = int(np.ceil(npanels / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes) # Ensure axes is always a 2D array even for a single row

    # Determine Y-axis label depending on the requested density type
    if dens_type == 'number':
        ylabel = r'$n_\mathrm{d}$ [cm$^{-3}$]'
    elif dens_type == 'surface':
        ylabel = r'Surfaces [cm$^{-1}$]'
    else:
        ylabel = r'$\rho_\mathrm{d}$ [g cm$^{-3}$]'

    # 6. Plotting loop over all available subplot slots
    for idx in range(nrows * ncols):
        ax = axes.flat[idx]
        
        if idx < nspecies:
            # Extract 1D radial profile at the midplane for the current dust species
            profile = dens[idx, idx_midplane, :]
            
            # Apply conversion factors based on the selected density type
            if dens_type == 'number':
                y_data = profile / grain_mass[idx] # Mass density / mass of one grain
            elif dens_type == 'surface':
                # Cross-sectional area calculation (converting size from micron to cm)
                y_data = 4 * np.pi * (sizes[idx] * 1e-4) * profile / grain_mass[idx]
            elif dens_type == 'mass':
                y_data = profile # Default is raw mass density
                
            ax.plot(r_center, y_data, color='darkblue', lw=2)
            ax.set_title(f'Bin {idx+1}', fontsize=12)
            
            # Add text box indicating the grain size for this specific bin
            if sizes is not None and idx < len(sizes):
                s = sizes[idx]
                size_label = f'{s/1e3:.1f} mm' if s >= 1e3 else f'{s:.2f} ' + r'$\mu$m'
                ax.text(0.05, 0.95, size_label, transform=ax.transAxes,
                        fontsize=12, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        elif idx == nspecies:
            # Plot total cumulative density (only relevant/calculated for mass density)
            if dens_type == 'mass':
                total_profile = dens[:, idx_midplane, :].sum(axis=0)
                ax.plot(r_center, total_profile, color='black', lw=2.5, linestyle='--')
                ax.set_title('Total Mass', fontsize=12)
            else:
                ax.axis('off') # Hide total panel if it's not mass density
        else:
            ax.axis('off') # Hide any remaining empty subplots in the grid

        # Configure axes scales and limits
        ax.set_yscale('log')
        ax.set_ylim(vmin, vmax)
        if xlim:
            ax.set_xlim(xlim)
        else:
            ax.set_xscale('log') # Logarithmic scale is standard for protoplanetary disks

    # Add global outer axis labels (only on edge plots thanks to sharex/sharey)
    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel('r [au]', fontsize=12)
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel, fontsize=12)

    fig.tight_layout()
    plt.show()


def density2D_grid_interactive(path='thermal/', vmin=1e-30, vmax=1e-15, cmap='gnuplot2', dens_type='mass',
                                xlim=None, ylim=None, dust=None, figsize=(10, 14)):
    """Interactive version of density2D_grid with sliders for vmin/vmax.
    Requires %matplotlib widget in the notebook."""
    import ipywidgets as widgets
    from IPython.display import display
    import os

    # --- Load data once ---
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[grid.columns[0]].values

    dens = pd.read_table(path + 'dust_density.inp', engine='python', header=None, skiprows=3)
    dens = dens[0].values
    nspecies = int(len(dens) / (nr * nt))
    dens = np.reshape(dens, (nspecies, nt, nr))

    r_edge = grid[:nr+1] / autocm
    theta_edge = grid[nr+1:nr+1+nt+1]
    theta_edge[-1] = np.pi
    rr_edge, tt_edge = np.meshgrid(r_edge, theta_edge)
    R = rr_edge * np.sin(tt_edge)
    Z = rr_edge * np.cos(tt_edge)

    dens[dens <= 1e-100] = 1e-100

    sizes_file = path + 'dust_sizes.inp'
    if os.path.isfile(sizes_file):
        sizes = np.loadtxt(sizes_file)
        sizes = np.atleast_1d(sizes)
    elif dust is not None:
        sizes = dust.sizes()[0]
        grain_mass = dust.grainmass()
    else:
        sizes = None

    # Precompute plot data for each panel
    plot_data = []
    for idx in range(nspecies):
        if dens_type == 'number' and dust is not None:
            plot_data.append(4 * np.pi * sizes[idx] * 1e-4 * dens[idx] / grain_mass[idx])
        else:
            plot_data.append(dens[idx])
    plot_data.append(dens.sum(axis=0))  # total

    npanels = nspecies + 1
    ncols = min(nspecies, 4)
    nrows = int(np.ceil(npanels / ncols))

    # --- Build figure ---
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    meshes = []
    for idx in range(nrows * ncols):
        ax = axes.flat[idx]
        if idx < nspecies:
            im = ax.pcolormesh(R, Z, plot_data[idx], cmap=cmap, shading='auto',
                               norm=LogNorm(vmin=vmin, vmax=vmax))
            meshes.append(im)
            ax.set_title(f'bin {idx+1}', fontsize=12)
            if sizes is not None and idx < len(sizes):
                s = sizes[idx]
                if s >= 1e3:
                    size_label = f'{s/1e3:.1f} mm'
                else:
                    size_label = f'{s:.2f} ' + r'$\mu$m'
                ax.text(0.05, 0.95, size_label, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top',
                        horizontalalignment='left', bbox=props)
        elif idx == nspecies:
            im = ax.pcolormesh(R, Z, plot_data[nspecies], cmap=cmap, shading='auto',
                               norm=LogNorm(vmin=vmin, vmax=vmax))
            meshes.append(im)
            ax.set_title('total', fontsize=12)
        else:
            ax.set_visible(False)
            continue
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel('r [au]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel('z [au]', fontsize=14)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(meshes[-1], cax=cbar_ax, label=r'$\rho_\mathrm{d}$ [g cm$^{-3}$]')

    # --- Sliders ---
    log_vmin = np.log10(vmin)
    log_vmax = np.log10(vmax)

    vmin_slider = widgets.FloatSlider(
        value=log_vmin, min=-50, max=0, step=0.5,
        description='log(vmin)', continuous_update=False,
        style={'description_width': 'initial'}, layout=widgets.Layout(width='400px'))
    vmax_slider = widgets.FloatSlider(
        value=log_vmax, min=-50, max=0, step=0.5,
        description='log(vmax)', continuous_update=False,
        style={'description_width': 'initial'}, layout=widgets.Layout(width='400px'))

    def update_clim(change):
        new_vmin = 10**vmin_slider.value
        new_vmax = 10**vmax_slider.value
        if new_vmin >= new_vmax:
            return
        new_norm = LogNorm(vmin=new_vmin, vmax=new_vmax)
        for m in meshes:
            m.set_norm(new_norm)
        cbar.update_normal(meshes[-1])
        fig.canvas.draw_idle()

    vmin_slider.observe(update_clim, names='value')
    vmax_slider.observe(update_clim, names='value')

    display(widgets.HBox([vmin_slider, vmax_slider]))
    plt.show()




def temperature2D_grid(path='thermal/', vmin=1e0, vmax=1e3, cmap='gnuplot2',
                    xlim=None, ylim=None, figsize=(10, 14)):
    """Plot all dust species on a single figure with subplots, plus total density."""
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[grid.columns[0]].values

    temp = pd.read_table(path + 'dust_temperature.dat', engine='python', header=None, skiprows=3)
    temp = temp[0].values
    nspecies = int(len(temp) / (nr * nt))
    temp = np.reshape(temp, (nspecies, nt, nr))
    grid = np.array(grid,copy=True)
    r_edge = grid[:nr+1] / autocm
    theta_edge = grid[nr+1:nr+1+nt+1]
    theta_edge[-1] = np.pi
    rr_edge, tt_edge = np.meshgrid(r_edge, theta_edge)
    R = rr_edge * np.sin(tt_edge)
    Z = rr_edge * np.cos(tt_edge)

    # Convert edge grid to cell-center grid for contour plotting
    R_center = 0.25 * (R[:-1, :-1] + R[1:, :-1] + R[:-1, 1:] + R[1:, 1:])
    Z_center = 0.25 * (Z[:-1, :-1] + Z[1:, :-1] + Z[:-1, 1:] + Z[1:, 1:])


        # Try to read grain sizes for subplot labels
    import os
    sizes_file = path + 'dust_sizes.inp'
    if os.path.isfile(sizes_file):
        sizes = np.loadtxt(sizes_file)
        sizes = np.atleast_1d(sizes)
    else:
        sizes = None

    # Layout: enough panels for all species + 1 total
    npanels = nspecies + 1
    ncols = min(nspecies, 4)
    nrows = int(np.ceil(npanels / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    levels = np.arange(vmin, vmax + 10, 10)

    for idx in range(nrows * ncols):
        ax = axes.flat[idx]
        if idx < nspecies:
            im = ax.pcolormesh(R, Z, temp[idx], cmap=cmap, shading='auto',
                              norm=LogNorm(vmin=vmin, vmax=vmax))
            #im = ax.contourf(R_center, Z_center, temp[idx], levels=levels, cmap=cmap)
            cs = ax.contour(
                R_center,
                Z_center,
                temp[idx],
                levels=[20],
                colors='black',
                linewidths=2
            )
            #ax.clabel(cs, fmt="20 K", fontsize=10)
            ax.set_title(f'bin {idx+1}', fontsize=12)
            if sizes is not None and idx < len(sizes):
                s = sizes[idx]
                if s >= 1e3:
                    size_label = f'{s/1e3:.1f} mm'
                else:
                    size_label = f'{s:.2f} ' + r'$\mu$m'
                ax.text(0.05, 0.95, size_label, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top',
                        horizontalalignment='left', bbox=props)
        elif idx == nspecies:
            total = temp.sum(axis=0)
            im = ax.pcolormesh(R, Z, total, cmap=cmap, shading='auto',
                               norm=LogNorm(vmin=vmin, vmax=vmax))
            ax.set_title('total', fontsize=12)
        else:
            ax.set_visible(False)
            continue

        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

    # Shared labels
    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel('r [au]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel('z [au]', fontsize=14)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label=r'T [K]')

    plt.show()


def midplane_temp(path='thermal/', xlim=None, ylim=None):
    grid = pd.read_table(path+'amr_grid.inp', engine='python', skiprows=5)
    head = grid.columns
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[head[0]].values
    try:
        temp = pd.read_table(path+'dust_temperature.dat', engine='python', header=None, skiprows=3)
    except IOError:
        print('plot.midplane_temp: the file dust_temperature.dat is not present. Run a dust thermal simulation first.')
        sys.exit(1)
    temp = temp[0].values
    nbspecies = int(len(temp)/(nr*nt))
    temp = np.reshape(temp, (nbspecies, nt, nr))
    grid = np.array(grid,copy=True)
    dist = grid[:nr+1]/autocm
    theta = grid[nr+1:nr+1+nt+1]
    theta[-1] = np.pi
    dist, tt = np.meshgrid(dist, theta)
    rr = dist*np.sin(tt)
    zz = dist*np.cos(tt)
    midtemp = temp[:, 90, :]
    radii = 0.5*(rr[90][0:rr[90].size-1] + rr[90][1:rr[90].size])

    #--PLOT FIGURE--
    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111)
    #-----profiles
    midtemp = pd.DataFrame(data=midtemp.transpose())
    for ispec in range(0, nbspecies):
        ax.plot(radii, midtemp[ispec].rolling(window=6, center=True).mean(), linewidth=2, linestyle='-', label='bin: {}'.format(ispec+1))
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
    ax.set_xlabel(r'r [au]', fontsize = 20)
    ax.set_ylabel(r'T [K]', fontsize = 20)
    ax.legend(fontsize=15)
    ax.tick_params(labelsize=18)
    plt.show()

def vertical_temp(thermpath='thermal/', chempath='chemistry/', r=100):
    grid = pd.read_table(thermpath+'amr_grid.inp', engine='python', skiprows=5)
    head = grid.columns
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[head[0]].values
    temp = pd.read_table(thermpath+'dust_temperature.dat', engine='python', header=None, skiprows=3)
    temp = temp[0].values
    nbspecies = int(len(temp)/(nr*nt))

    if nbspecies == 1:
        try:
            temp = pd.read_table(chempath+str(r)+'AU/1D_static.dat', sep=r"\s+", engine='python', header=None, comment='!')
        except IOError:
            print('plot.vertical_temp: radius {} does not exit in the model or path is not correct.'.format(r))
            sys.exit(1)
        #--PLOT FIGURE--
        fig = plt.figure(figsize=(9.6, 8.2))
        ax = fig.add_subplot(111)
        #-----profiles
        ax.plot(temp[5], temp[0], linewidth=2, linestyle='-', label='{} AU'.format(r))
        # ax.set_ylim(0,60)
        # ax.set_xlim(1,350)
        ax.set_ylabel(r'z [au]', fontsize = 20)
        ax.set_xlabel(r'T$_\mathrm{d}$ [K]', fontsize = 20)
        ax.legend(fontsize=15)
        ax.tick_params(labelsize=18)
        plt.show()
    elif nbspecies > 1:
        try:
            static = pd.read_table(chempath+str(r)+'AU/1D_static.dat', sep=r"\s+", engine='python', header=None, comment='!')
            temp = pd.read_table(chempath+str(r)+'AU/temperatures.dat', sep=r"\s+", engine='python', header=None)
        except IOError:
            print('plot.vertical_temp: radius = {} au does not exit in the model or path is not correct.'.format(r))
            sys.exit(1)
        #--PLOT FIGURE--
        fig = plt.figure(figsize=(9.6, 8.2))
        ax = fig.add_subplot(111)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.91, 0.05, '{} AU'.format(r), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=16, bbox=props)
        #-----profiles
        for ai in range(nbspecies):
            ax.plot(temp[ai], static[0], linewidth=2, linestyle='-', label='bin: {}'.format(ai+1))
        # ax.set_ylim(0,60)
        # ax.set_xlim(1,350)
        ax.set_ylabel(r'z [au]', fontsize = 20)
        ax.set_xlabel(r'T$_\mathrm{d}$ [K]', fontsize = 20)
        ax.legend(fontsize=15)
        ax.tick_params(labelsize=18)
        plt.show()

def avz(chempath='thermal/', r=100):
    static = pd.read_table(chempath+str(r)+'AU/1D_static.dat', sep=r"\s+", engine='python', header=None, comment='!', skiprows=1)
    #--PLOT FIGURE--
    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.91, 0.05, '{} AU'.format(r), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=16, bbox=props)
    #-----profiles
    ax.plot(static[3], static[0], linewidth=2, linestyle='-', label='vertical Av')
    # ax.set_ylim(0,60)
    # ax.set_xlim(1,350)
    ax.set_xlabel(r'z [au]', fontsize = 20)
    ax.set_ylabel(r'A$_\mathrm{\nu}$ [mag]', fontsize = 20)
    ax.legend(fontsize=15)
    ax.tick_params(labelsize=18)
    plt.show()

def opacity(path='thermal/'):
    opaclist = sorted(glob.glob(path+'dustkap*'))

    #---absorption
    fig = plt.figure(figsize=(9.6, 8.2)) #fig = plt.figure(figsize=(9.6, 7.2))
    ax = fig.add_subplot(111) 
    ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize=18)
    ax.set_ylabel(r'$\kappa_\mathrm{abs}$ [cm$^2$/g]', fontsize=18)
    ax.set_xlim(1e-1,1e4)
    ax.set_ylim(1e-2,1e5)
    for opac in opaclist:
        name = opac.split("_")[1].split(".")[0]
        kappa = pd.read_table(opac, sep=r"\s+", comment='#', header=None, skiprows=10)
        ax.loglog(kappa[0], kappa[1], linewidth=2, label=name)
    ax.tick_params(labelsize=22)
    ax.legend(fontsize=15)
    plt.show()

    #---scattering
    fig = plt.figure(figsize=(9.6, 8.2)) #fig = plt.figure(figsize=(9.6, 7.2))
    ax = fig.add_subplot(111) 
    ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize=18)
    ax.set_ylabel(r'$\kappa_\mathrm{scat}$ [cm$^2$/g]', fontsize=18)
    ax.set_xlim(1e-1,1e4)
    ax.set_ylim(1e-2,1e5)
    for opac in opaclist:
        name = opac.split("_")[1].split(".")[0]
        kappa = pd.read_table(opac, sep=r"\s+", comment='#', header=None, skiprows=10)
        ax.loglog(kappa[0], kappa[2], linewidth=2, label=name)
    ax.tick_params(labelsize=22)
    ax.legend(fontsize=15)
    plt.show()

    #---angles
    fig = plt.figure(figsize=(9.6, 8.2)) #fig = plt.figure(figsize=(9.6, 7.2))
    ax = fig.add_subplot(111) 
    ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize=18)
    ax.set_ylabel(r'<cos($\theta$)>', fontsize=18)
    ax.set_xlim(1e-1,1e4)
    #ax.set_ylim(0,1)
    for opac in opaclist:
        name = opac.split("_")[1].split(".")[0]
        kappa = pd.read_table(opac, sep=r"\s+", comment='#', header=None, skiprows=10)
        ax.loglog(kappa[0], kappa[3], linewidth=2, label=name)
    ax.tick_params(labelsize=22)
    ax.legend(fontsize=15)
    plt.show()


def localflux(path='thermal/'):
    #---1/ Get grid shape and reshape the local flux array accordingly
    flux = pd.read_table(path+'mean_intensity.out', sep=r"\s+", comment='#', header=None, skiprows=4)
    grid = pd.read_table(path+'amr_grid.inp', engine='python', skiprows=5)
    lam = pd.read_table(path+'mcmono_wavelength_micron.inp', engine='python', header=None, skiprows=1)
    lam = lam[0].values
    flux = flux[0].values
    
    head = grid.columns
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[head[0]].values
    dist = grid[:nr+1]/autocm
    theta = grid[nr+1:nr+1+nt+1]
    theta[-1] = np.pi
    dist, tt = np.meshgrid(dist, theta)
    rr = dist*np.sin(tt)
    radii = 0.5*(rr[90][0:rr[90].size-1] + rr[90][1:rr[90].size])
    zz = dist*np.cos(tt)
    nlam = int(len(flux)/(nr*nt))
    flux = np.reshape(flux, (nlam, nt, nr))
    midflux = flux[:, 90, :]

    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111)
    #-----profiles
    midflux_df = pd.DataFrame(data=midflux.transpose())
    for ilam in range(0, nlam, 2):
        ax.semilogy(radii, midflux_df[ilam].rolling(window=5, center=True).mean(), linewidth=1, linestyle='-')
    ax.set_xlim(1,200)
    ax.set_ylim(1e-30,1e-10)
    ax.set_xlabel(r'r [au]', fontsize = 20)
    ax.set_ylabel(r'Flux', fontsize = 20)
    ax.grid()
    ax.tick_params(labelsize=22)
    plt.show()


def image(pathfile='thermal/', distance=100, vmin=1e-10, vmax=1e3, cmap='gnuplot2',
          labels=None):
    """Plot a multi-wavelength RADMC3D continuum image (Stokes I).

    Parameters
    ----------
    pathfile : str
        Path to the RADMC3D image.out file.
    distance : float, optional
        Source distance in parsec. Used to convert specific intensity
        [erg/s/cm²/Hz/sr] to flux density [Jy/pixel]. Default is 100 pc.
    vmin, vmax : float, optional
        Color scale limits in Jy/pixel.
    cmap : str, optional
        Colormap name.
    labels : list of str, optional
        Wavelength labels for each panel. If None, labels are generated
        automatically from the wavelengths read in the image header.
    """
    # --- Read RADMC3D image header ---
    with open(pathfile, 'r') as f:
        iformat = int(f.readline())                                 # 1 = I only, 3 = full Stokes
        npix_x, npix_y = [int(x) for x in f.readline().split()]   # image size [pixels]
        nlam = int(f.readline())                                    # number of wavelengths
        pix_cm, _ = [float(x) for x in f.readline().split()]       # pixel size [cm]
        wavelengths = [float(f.readline()) for _ in range(nlam)]   # wavelengths [microns]

    # iformat 3 → full Stokes (I Q U V); anything else → intensity only
    nstokes = 4 if iformat == 3 else 1

    pix_au   = pix_cm / autocm
    box_au   = npix_x * pix_au
    half_box = box_au / 2.0

    # --- Pixel solid angle and Jy/pixel conversion factor ---
    distance_cm = distance * 3.086e18               # pc → cm
    omega_pix   = (pix_cm / distance_cm) ** 2       # sr/pixel
    to_jy       = 1e23 * omega_pix                  # erg/s/cm²/Hz/sr → Jy/pixel

    data = np.loadtxt(pathfile, skiprows=4 + nlam + 1)
    data = np.reshape(data, (nlam, npix_y, npix_x, nstokes))

    extent = [-half_box, half_box, -half_box, half_box]

    # --- Stokes I images in Jy/pixel ---
    imgs = [data[i, :, :, 0] * to_jy for i in range(nlam)]

    # --- Wavelength labels ---
    if labels is None:
        def _wave_label(lam):
            if lam >= 1000:
                return f'Stokes I - {lam/1000:.2g} mm'
            elif lam >= 1:
                return f'Stokes I - {lam:.4g} μm'
            else:
                return f'Stokes I - {lam*1000:.4g} nm'
        labels = [_wave_label(w) for w in wavelengths]

    fig, axes = plt.subplots(1, nlam, figsize=(5*nlam, 5), sharex=True, sharey=True)
    if nlam == 1:
        axes = [axes]

    for i, ax in enumerate(axes):

        im = ax.imshow(
            imgs[i],
            origin='lower',
            extent=extent,
            cmap=cmap,
            norm=LogNorm(vmin=vmin, vmax=vmax),
            interpolation='nearest'
        )

        ax.tick_params(labelsize=17)
        ax.set_xlabel(r'x [au]', fontsize=17)

        if i == 0:
            ax.set_ylabel(r'y [au]', fontsize=17)

        ax.text(
            0.05, 0.95,
            labels[i],
            horizontalalignment='left',
            verticalalignment='top',
            color='red',
            transform=ax.transAxes,
            fontsize=16,
            fontweight='bold'
        )

    # --- Proper colorbar placement ---
    cbar = fig.colorbar(
        im,
        ax=axes,
        location='right',
        fraction=0.025,
        pad=0.02
    )

    cbar.set_label(r'$I_\nu$ [Jy pixel$^{-1}$]', fontsize=17)
    cbar.ax.tick_params(labelsize=14)

    plt.show()


def image_vertical_cut(pathfile='thermal/', distance=100, xlim=None, ylim=None,
                       labels=None, figsize=(9.6, 8.2)):
    """Plot vertical cuts (along y at x=0) of Stokes I for each wavelength.

    Parameters
    ----------
    pathfile : str
        Path to the RADMC3D image.out file.
    distance : float, optional
        Source distance in parsec. Used to convert specific intensity
        [erg/s/cm²/Hz/sr] to flux density [Jy/pixel]. Default is 100 pc.
    xlim, ylim : tuple, optional
        Axis limits.
    labels : list of str, optional
        Wavelength labels. If None, generated automatically from the header.
    figsize : tuple, optional
        Figure size.
    """
    # --- Read RADMC3D image header ---
    with open(pathfile, 'r') as f:
        iformat = int(f.readline())                                 # 1 = I only, 3 = full Stokes
        npix_x, npix_y = [int(x) for x in f.readline().split()]   # image size [pixels]
        nlam = int(f.readline())                                    # number of wavelengths
        pix_cm, _ = [float(x) for x in f.readline().split()]       # pixel size [cm]
        wavelengths = [float(f.readline()) for _ in range(nlam)]   # wavelengths [microns]

    nstokes = 4 if iformat == 3 else 1

    pix_au   = pix_cm / autocm
    box_au   = npix_y * pix_au
    half_box = box_au / 2.0

    # --- Pixel solid angle and Jy/pixel conversion factor ---
    distance_cm = distance * 3.086e18
    omega_pix   = (pix_cm / distance_cm) ** 2
    to_jy       = 1e23 * omega_pix

    data = np.loadtxt(pathfile, skiprows=4 + nlam + 1)
    data = np.reshape(data, (nlam, npix_y, npix_x, nstokes))

    y_au = np.linspace(-half_box, half_box, npix_y)
    ix0  = npix_x // 2  # column at x=0

    if labels is None:
        def _wave_label(lam):
            if lam >= 1000:
                return f'{lam/1000:.2g} mm'
            elif lam >= 1:
                return f'{lam:.4g} μm'
            else:
                return f'{lam*1000:.4g} nm'
        labels = [_wave_label(w) for w in wavelengths]

    fig, ax = plt.subplots(figsize=figsize)

    for i in range(nlam):
        flux_cut = data[i, :, ix0, 0] * to_jy
        ax.semilogy(y_au, flux_cut, linewidth=2, label=labels[i])

    ax.set_xlabel(r'y [au]', fontsize=20)
    ax.set_ylabel(r'$I_\nu$ [Jy pixel$^{-1}$]', fontsize=20)
    ax.legend(fontsize=15)
    ax.tick_params(labelsize=18)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    plt.show()

def numberdens(species='CO', path='thermal/', vmin=1e0, vmax=1e8, cmap='gnuplot2',
               ncols=3, xlim=None, ylim=None, figsize=None,
               save=False, savename='numberdens.pdf'):
    """Plot 2D maps of molecular number densities from ``numberdens_XXX.inp`` files.

    Accepts a single species name or a list of names. Multiple species are
    displayed as a mosaic of subplots sharing the same colour scale.

    Parameters
    ----------
    species : str or list of str, optional
        Species name(s) matching ``numberdens_<species>.inp`` files.
        Default is ``'CO'``.
    path : str, optional
        Path to the directory containing the RADMC-3D files. Default is
        ``'thermal/'``.
    vmin, vmax : float, optional
        Shared colour scale limits [cm :sup:`-3`]. Default is ``1e0``
        and ``1e8``.
    cmap : str, optional
        Colormap name. Default is ``'gnuplot2'``.
    ncols : int, optional
        Maximum number of columns in the mosaic. Default is ``3``.
    xlim, ylim : tuple, optional
        Axis limits (R, Z) in AU applied to every panel.
    figsize : tuple or None, optional
        Figure size. If None, computed automatically from ``ncols`` and
        the number of rows.
    save : bool, optional
        Save the figure to ``savename``. Default is False.
    savename : str, optional
        Output filename when ``save=True``. Default is ``'numberdens.pdf'``.
    """
    species_list = [species] if isinstance(species, str) else list(species)
    nspecies = len(species_list)

    # --- Grid: cell centres (shape nt × nr) — required for shading='gouraud' ---
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = np.array(grid[grid.columns[0]].values, copy=True)

    r_edge     = grid[:nr+1] / autocm
    theta_edge = grid[nr+1:nr+1+nt+1]
    theta_edge[-1] = np.pi
    r_cen     = 0.5 * (r_edge[:-1]     + r_edge[1:])
    theta_cen = 0.5 * (theta_edge[:-1] + theta_edge[1:])
    rr, tt = np.meshgrid(r_cen, theta_cen)          # (nt, nr)
    R = rr * np.sin(tt)
    Z = rr * np.cos(tt)

    # --- Layout ---
    ncols = min(ncols, nspecies)
    nrows = int(np.ceil(nspecies / ncols))
    if figsize is None:
        figsize = (5 * ncols + 1, 4 * nrows)

    norm = LogNorm(vmin=vmin, vmax=vmax)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()

    im = None
    for idx, sp in enumerate(species_list):
        ax = axes[idx]
        nd = pd.read_table(path + f'numberdens_{sp}.inp',
                           engine='python', header=None, skiprows=2)
        nd = nd[0].values.reshape(nt, nr)
        nd = np.where(nd <= 0, 1e-100, nd)

        im = ax.pcolormesh(R, Z, nd, cmap=cmap, shading='gouraud',
                           norm=norm, rasterized=True)
        ax.set_title(sp, fontsize=13)
        ax.tick_params(labelsize=12)
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

    # Hide unused panels
    for idx in range(nspecies, len(axes)):
        axes[idx].set_visible(False)

    # Shared axis labels
    for ax in axes[(nrows - 1) * ncols:]:
        ax.set_xlabel('R [au]', fontsize=13)
    for i in range(nrows):
        axes[i * ncols].set_ylabel('Z [au]', fontsize=13)

    # Single shared colorbar on the right
    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label=r'$n$ [cm$^{-3}$]')

    if save:
        fig.savefig(savename, bbox_inches='tight')

    plt.show()


def plot_velocity_and_temperature(path='thermal/', 
                                  vmin=0.0, vmax=10.0, logscale=False, cmap_v='viridis',
                                  Tmin=10.0, Tmax=1000.0, logscale_T=True, cmap_T='inferno',
                                  xlim=None, ylim=None, figsize=None,
                                  save=False, savename='gas_properties.pdf'):
    """Plot a 2D side-by-side comparison of gas velocity (v_phi) and gas temperature.

    Extracts the azimuthal Keplerian velocity component from a 3D spherical 
    velocity file and the thermal gas structure from a temperature file, both 
    formatted for RADMC-3D line transfer calculations. Displays them as a 
    meridional (R, Z) cross-section slice.

    Parameters
    ----------
    path : str, optional
        Path to the directory containing the RADMC-3D input files. Default is 
        ``'thermal/'``.
    vmin, vmax : float, optional
        Colour scale limits for the azimuthal velocity v_phi [km/s]. Default 
        is ``0.0`` and ``10.0``.
    logscale : bool, optional
        If True, plot the velocity using a logarithmic color scale. Default 
        is False.
    cmap_v : str, optional
        Colormap name for velocity field. Default is ``'viridis'``.
    Tmin, Tmax : float, optional
        Colour scale limits for gas kinetic temperature [K]. Default is 
        ``10.0`` and ``1000.0``.
    logscale_T : bool, optional
        If True, plot the temperature using a logarithmic color scale. Default 
        is True.
    cmap_T : str, optional
        Colormap name for the thermal structure. Default is ``'inferno'``.
    xlim, ylim : tuple, optional
        Spatial axis limits (R, Z) in astronomical units [au], shared by 
        both panels.
    figsize : tuple or None, optional
        Figure dimensions. Default is ``(12, 5)`` for side-by-side layout.
    save : bool, optional
        Save the rendered figure to ``savename``. Default is False.
    savename : str, optional
        Output filename when ``save=True``. Default is ``'gas_properties.pdf'``.
    """
    # --- Grid: cell centres (shape nt × nr) ---
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = np.array(grid[grid.columns[0]].values, copy=True)

    # Conversion of spherical radial grid boundaries from cm to au
    r_edge     = grid[:nr+1] / autocm
    theta_edge = grid[nr+1:nr+1+nt+1]
    theta_edge[-1] = np.pi
    
    # Calculation of cell centres from boundary edges via midpoints
    r_cen     = 0.5 * (r_edge[:-1]     + r_edge[1:])
    theta_cen = 0.5 * (theta_edge[:-1] + theta_edge[1:])
    rr, tt = np.meshgrid(r_cen, theta_cen)          # (nt, nr)
    
    # Transformation from spherical polar (r, theta) to cylindrical (R, Z) coordinates
    R = rr * np.sin(tt)
    Z = rr * np.cos(tt)

    if figsize is None:
        figsize = (12, 5)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, sharex=True, sharey=True)

    if logscale:
        if vmin <= 0: vmin = 1e-2
        norm_v = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm_v = Normalize(vmin=vmin, vmax=vmax)

    # RADMC-3D velocity file contains columns for (v_r, v_theta, v_phi)
    vel_data = pd.read_table(path + 'gas_velocity.inp', engine='python', 
                             header=None, skiprows=2, sep=r'\s+')
    
    # Extracting the 3rd column (index 2) for v_phi and converting from cm/s to km/s
    vphi_all = vel_data[2].values / 1e5
    
    # Determining azimuthal grid size to reconstruct the nested loops (phi, theta, r)
    nphi = len(vphi_all) // (nt * nr)
    vphi_3d = vphi_all.reshape(nphi, nt, nr)
    vphi_2d = vphi_3d[0, :, :]  # Meridional slice at phi = 0

    im1 = ax1.pcolormesh(R, Z, vphi_2d, cmap=cmap_v, shading='gouraud',
                         norm=norm_v, rasterized=True)
    
    title_v = r'Gas Velocity $v_\phi$ (Log)' if logscale else r'Gas Velocity $v_\phi$'
    ax1.set_title(title_v, fontsize=13)
    ax1.set_xlabel('R [au]', fontsize=12)
    ax1.set_ylabel('Z [au]', fontsize=12)
    ax1.tick_params(labelsize=11)

    temp_data = pd.read_table(path + 'gas_temperature.inp', engine='python', 
                              header=None, skiprows=2)
    t_all = temp_data[0].values
    
    # Reconstructing the 3D grid layout matching the nested cell indexing (phi, theta, r)
    t_3d = t_all.reshape(nphi, nt, nr)
    t_2d = t_3d[0, :, :]  # Meridional slice at phi = 0

    if Tmin is None:
        # Exclusion of unphysical values (<= 0) when evaluating bounds under LogNorm
        Tmin = np.min(t_2d[t_2d > 0]) if logscale_T else np.min(t_2d)
    if Tmax is None:
        Tmax = np.max(t_2d)

    if logscale_T:
        if Tmin <= 0: 
            Tmin = 1e-1
        norm_T = LogNorm(vmin=Tmin, vmax=Tmax)
    else:
        norm_T = Normalize(vmin=Tmin, vmax=Tmax)

    im2 = ax2.pcolormesh(R, Z, t_2d, cmap=cmap_T, shading='gouraud',
                         norm=norm_T, rasterized=True)
    
    title_T = 'Gas Temperature'
    ax2.set_title(title_T, fontsize=13)
    ax2.set_xlabel('R [au]', fontsize=12)
    ax2.tick_params(labelsize=11)

    if xlim:
        ax1.set_xlim(xlim)
    if ylim:
        ax1.set_ylim(ylim)

    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.07)
    fig.colorbar(im1, cax=cax1, label=r'$v_\phi$ [km s$^{-1}$]')

    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes("right", size="5%", pad=0.07)
    fig.colorbar(im2, cax=cax2, label=r'$T_{\rm gas}$ [K]')

    fig.tight_layout()

    if save:
        fig.savefig(savename, bbox_inches='tight')

    plt.show()


def static(chempath='chemistry/', column='nH', vmin=1, vmax=50, iso=None, cmap='gnuplot2',
           xlim=None, ylim=None, figsize=(6, 6), save=False, savename='filename.pdf'):
    """Plot a 2D map of a column from the 1D_static.dat files.

    Scans all XXAU/ folders in chempath, reads each 1D_static.dat,
    and builds a 2D (r, z) map of the chosen column.

    Parameters
    ----------
    chempath : str
        Path to the chemistry directory containing the XXAU/ folders.
    column : str
        Column to plot. One of: 'z', 'nH', 'Tg', 'Av', 'diff', 'Td',
        'inv_ab', 'conv_factor', 'a', 'uv'.
    vmin, vmax : float
        Color scale limits.
    iso : float or list, optional
        Draw contour lines of Td at these levels (e.g. 20 or [20, 50]).
    cmap : str
        Colormap name.
    xlim, ylim : tuple, optional
        Axis limits (r, z) in AU.
    figsize : tuple
        Figure size.
    save : bool
        Save the figure to savename.
    savename : str
        Output filename if save is True.
    """
    import os, re

    columns = ['z', 'nH', 'Tg', 'Av', 'diff', 'Td', 'inv_ab', 'conv_factor', 'a', 'uv']

    # Discover XXAU/ folders and extract radii
    folders = [d for d in os.listdir(chempath)
               if os.path.isdir(os.path.join(chempath, d)) and re.match(r'^\d+AU$', d)]
    rchem = sorted([int(d.replace('AU', '')) for d in folders])

    # Read all files at once (nbz may differ per radius after surface truncation)
    all_data = []
    for r in rchem:
        filepath = os.path.join(chempath, f'{r}AU', '1D_static.dat')
        df = pd.read_table(filepath, sep=r"\s+", comment='!', header=None, engine='python')
        df.columns = columns
        all_data.append(df)

    nbz_max = max(len(d) for d in all_data)

    # Build 2D arrays; NaN for cells above the truncation height of each radius
    static_map = np.full((nbz_max, len(rchem)), np.nan)
    temp_map   = np.full((nbz_max, len(rchem)), np.nan)
    zz         = np.zeros((nbz_max, len(rchem)))

    for idx, data in enumerate(all_data):
        nbz_r = len(data)
        start = nbz_max - nbz_r       # top rows belong to the truncated surface
        static_map[start:, idx] = data[column].values
        temp_map[start:, idx]   = data['Td'].values
        z_col = data['z'].values
        zz[start:, idx] = z_col
        # Extrapolate z above the truncation so pcolormesh has valid coordinates
        # (those cells are NaN in data so they will appear transparent)
        if start > 0:
            dz = (z_col[0] - z_col[1]) if nbz_r > 1 else z_col[0] * 0.1
            zz[:start, idx] = z_col[0] + np.arange(start, 0, -1) * dz

    rr, _ = np.meshgrid(rchem, np.arange(nbz_max))

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal', adjustable='box')
    t = ax.pcolormesh(rr, zz, static_map, cmap=cmap, shading='auto',
                      norm=LogNorm(vmin=vmin, vmax=vmax), rasterized=True)
    clr = fig.colorbar(t, pad=0.01)
    clr.set_label(column, fontsize=16)
    clr.ax.tick_params(labelsize=14)

    if iso is not None:
        if not isinstance(iso, (list, tuple)):
            iso = [iso]
        ax.contour(rr, zz, temp_map, iso, colors='black', linewidths=2.5)

    ax.set_xlabel(r'r [au]', fontsize=20)
    ax.set_ylabel(r'z [au]', fontsize=20)
    ax.tick_params(labelsize=15)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    if save:
        fig.savefig(savename, bbox_inches='tight')

    plt.show()


def nmgc_grainsizes(chempath='chemistry/', quantity='Td', vmin=None, vmax=None, cmap='gnuplot2',
                xlim=None, ylim=None, figsize=(14, 10), save=False, savename='grain_sizes.pdf'):
    """Plot 2D (r, z) maps per grain size from the 1D_grain_sizes.in files.

    Each subplot corresponds to one grain size bin.

    Parameters
    ----------
    chempath : str
        Path to the chemistry directory containing the XXAU/ folders.
    quantity : str
        What to plot: 'Td' for dust temperature, 'nd' for dust number density (nH / inv_ab).
    vmin, vmax : float, optional
        Color scale limits. Auto-determined if None.
    cmap : str
        Colormap name.
    xlim, ylim : tuple, optional
        Axis limits (r, z) in AU.
    figsize : tuple
        Figure size.
    save : bool
        Save the figure to savename.
    savename : str
        Output filename if save is True.
    """
    import os, re

    static_columns = ['z', 'nH', 'Tg', 'Av', 'diff', 'Td', 'inv_ab', 'conv_factor', 'a', 'uv']

    # Discover XXAU/ folders and extract radii
    folders = [d for d in os.listdir(chempath)
               if os.path.isdir(os.path.join(chempath, d)) and re.match(r'^\d+AU$', d)]
    rchem = sorted([int(d.replace('AU', '')) for d in folders])

    # Read the first grain_sizes file to determine N (number of grain sizes)
    # and grain radii from the first data line
    first_gs = os.path.join(chempath, f'{rchem[0]}AU', '1D_grain_sizes.in')
    with open(first_gs, 'r') as f:
        for line in f:
            stripped = line.split('!')[0].strip()
            if not stripped:
                continue
            vals = stripped.split()
            ncols = len(vals)
            grain_radii_cm = np.array([float(v) for v in vals[:ncols // 4]])
            break
    # ncols = 4*N (sizes, inv_ab, Td, CR-peak)
    ngrains = ncols // 4
    grain_radii_um = grain_radii_cm * 1e4  # cm to microns

    # Read all static files first (nbz may differ per radius after surface truncation)
    all_static = []
    for r in rchem:
        static_file = os.path.join(chempath, f'{r}AU', '1D_static.dat')
        sd = pd.read_table(static_file, sep=r"\s+", comment='!', header=None, engine='python')
        sd.columns = static_columns
        all_static.append(sd)

    nbz_max = max(len(sd) for sd in all_static)

    # Build arrays; NaN for cells above the truncation height of each radius
    data_map = np.full((ngrains, nbz_max, len(rchem)), np.nan)
    zz = np.zeros((nbz_max, len(rchem)))

    for idx, r in enumerate(rchem):
        static_data = all_static[idx]
        nbz_r = len(static_data)
        start = nbz_max - nbz_r       # top rows belong to the truncated surface
        nH    = static_data['nH'].values
        z_col = static_data['z'].values
        zz[start:, idx] = z_col
        # Extrapolate z above the truncation so pcolormesh has valid coordinates
        if start > 0:
            dz = (z_col[0] - z_col[1]) if nbz_r > 1 else z_col[0] * 0.1
            zz[:start, idx] = z_col[0] + np.arange(start, 0, -1) * dz

        # Read grain_sizes
        gs_file = os.path.join(chempath, f'{r}AU', '1D_grain_sizes.in')
        gs_lines = []
        with open(gs_file, 'r') as f:
            for line in f:
                stripped = line.split('!')[0].strip()
                if not stripped:
                    continue
                gs_lines.append([float(v) for v in stripped.split()])
        gs_array = np.array(gs_lines)  # (nbz_r, 4*ngrains)

        inv_ab = gs_array[:, ngrains:2*ngrains]       # (nbz_r, ngrains)
        Td     = gs_array[:, 2*ngrains:3*ngrains]     # (nbz_r, ngrains)

        if quantity == 'Td':
            for ig in range(ngrains):
                data_map[ig, start:, idx] = Td[:, ig]
        elif quantity == 'nd':
            for ig in range(ngrains):
                data_map[ig, start:, idx] = nH / inv_ab[:, ig]

    rr, _ = np.meshgrid(rchem, np.arange(nbz_max))

    # Layout
    ncols_plot = min(ngrains, 4)
    nrows_plot = int(np.ceil(ngrains / ncols_plot))

    fig, axes = plt.subplots(nrows_plot, ncols_plot, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    for ig in range(nrows_plot * ncols_plot):
        ax = axes.flat[ig]
        if ig < ngrains:
            ax.set_aspect('equal', adjustable='box')
            im = ax.pcolormesh(rr, zz, data_map[ig], cmap=cmap, shading='gouraud',
                               norm=LogNorm(vmin=vmin, vmax=vmax), rasterized=True)
            # Size label
            s = grain_radii_um[ig]
            if s >= 1e3:
                size_label = f'{s/1e3:.1f} mm'
            else:
                size_label = f'{s:.2f} ' + r'$\mu$m'
            ax.text(0.05, 0.95, size_label, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top',
                    horizontalalignment='left', bbox=props)
            if xlim:
                ax.set_xlim(xlim)
            if ylim:
                ax.set_ylim(ylim)
        else:
            ax.set_visible(False)

    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel('r [au]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel('z [au]', fontsize=14)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    if quantity == 'Td':
        label = r'T$_\mathrm{d}$ [K]'
    elif quantity == 'nd':
        label = r'n$_\mathrm{d}$ [cm$^{-3}$]'
    else:
        label = quantity
    fig.colorbar(im, cax=cbar_ax, label=label)

    if save:
        fig.savefig(savename, bbox_inches='tight')

    plt.show()



import os
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

def plot_outputs_nautilus(PIPE,
                          MODEL_NAMES,
                          itime=-1,
                          MODE='chemistry',
                          key_list=['CO'],
                          fracab=True,
                          verbose=True,
                          xlim=None,
                          ylim=None,
                          colormap="plasma",
                          res_colormap="seismic",
                          vmin=None,
                          vmax=None,
                          common_scale=True):
    r"""
    Plots a 2D vertical cross-section (poloidal cut in R, z) of Nautilus simulation outputs.

    This function dynamically handles the comparison and visualization of asymmetric 
    or non-uniform model grids by converting stacked 1D columns into 2D geometric 
    polygons (`PolyCollection`). It automatically builds and scales the subplot 
    canvas layout depending on the number of models provided.

    Subplot Layout Behavior
    -----------------------
    - 1 Model: Displays an optimized standard grid layout (max 3 columns) 
      distributing the requested chemical species or physical keys.
    - 2 Models: Forces a rigid 3-column matrix layout per row (key): 
      [Model 1] [Model 2] [Linear Residual Mapping (Model 1 - Model 2)]
    - >2 Models: Renders a strict multi-column grid where each column represents 
      a single unique model, and each row corresponds to a specific key/species.

    Parameters
    ----------
    PIPE : list
        A collection of model pipe objects containing the structured simulation outputs.
        Each object must feature a `.chemistry` attribute (nested dictionary or 
        xarray Dataset indexed by radius), a `.chempath` attribute (string/Path 
        pointing to raw files), and optionally a `.name` attribute.
    MODEL_NAMES : list of str
        Default display names assigned to the models for plot headers. If the list 
        contains duplicate names, a generic fallback nomenclature ("Model 1", 
        "Model 2", etc.) is automatically generated to guarantee uniqueness.
    itime : int, optional
        Simulation timestep index to extract and visualize. 
        Defaults to -1 (the final recorded timestep).
    MODE : {'chemistry', 'physical'}, optional
        The nature of the variables to extract and map:
        - 'chemistry': Extracts abundances or population densities from the 
          xarray Dataset or dict under the 'abundances' key.
        - 'physical': Extracts direct macroscopic variables (e.g., 'gas_temperature',
          'dust_temperature', 'extinction').
        Defaults to 'chemistry'.
    key_list : str or list of str, optional
        A single identifier string or a list of keys representing the chemical 
        species (e.g., 'CO', 'JCO', 'K2CO') or physical variables to map.
        Defaults to ['CO'].
    fracab : bool, optional
        Used exclusively when MODE='chemistry'. 
        If True, plots fractional abundances relative to total hydrogen (n_X/n_H).
        If False, computes and plots absolute number densities in cm^-3 
        (n_X = abundance * n_H).
        Defaults to True.
    verbose : bool, optional
        If True, prints diagnostic warnings and alerts in the console regarding 
        missing structural configuration files ('1D_static.dat') or shape mismatches 
        between data arrays. Defaults to True.
    xlim : tuple of float, optional
        Custom manual boundaries (min, max) to enforce on the radial axis R [AU]. 
        If None, the axis boundaries dynamically fit up to 102% of the 
        maximum detected radius.
    ylim : tuple of float, optional
        Custom manual boundaries (min, max) to enforce on the altitude axis z [AU]. 
        If None, the axis boundaries dynamically fit up to 107% of the 
        maximum detected altitude.
    colormap : str, optional
        The name of the Matplotlib sequential colormap used to render the main 
        simulation profile panels. Defaults to "plasma".
    res_colormap : str, optional
        The name of the Matplotlib diverging colormap used to render the residual 
        discrepancy panels (only triggered when exactly two models are present). 
        Defaults to "seismic".
    vmin : float, optional
        Forced lower bound for the colorbar sampling scale. If None, the minimum 
        bound is calculated dynamically from the underlying dataset.
    vmax : float, optional
        Forced upper bound for the colorbar sampling scale. If None, the maximum 
        bound is calculated dynamically from the underlying dataset.
    common_scale : bool, optional
        If True, evaluates a unified, global colorbar boundary framework (`vmin`, 
        `vmax`) shared across ALL models AND ALL keys simultaneously. If False, 
        each individual subplot panel normalizes its own independent localized 
        color scale. Defaults to True.

    Returns
    -------
    None
        The function assembles, structures, and renders an interactive Matplotlib 
        figure window (`plt.show()`) without returning any object in memory.

    Notes
    -----
    - **Logarithmic Scaling**: The rendering pipeline automatically switches to a 
      logarithmic normalization (`LogNorm`) if `MODE='chemistry'` or if the keywords 
      "density" or "extinction" are detected anywhere inside the target key string.
    - **2D Mesh Synthesis**: Grid cell edge boundaries are mathematically 
      extrapolated midway between the discrete points defined in the `1D_static.dat` 
      files (altitudes) and the parent radius keys. This guarantees clean, non-biased 
      polygonal boundaries across highly non-uniform scientific grids.
    - **Grain Environment Tracking**: Prefixes 'J' (surface phase) and 'K' (mantle 
      phase) followed by numerical indices are automatically decoded. The routine 
      attempts to open `1D_grain_sizes.in` to dynamically translate the bin index into 
      a physical grain radius ($\mu m$) displayed directly in the plot subtitles.
    - **Residual Engine**: The residual maps generated for the 2-model pipeline compute 
      a raw, linear arithmetic difference ($Model_1 - Model_2$). The resulting 
      diverging color scale is symmetrically locked around zero (`-max_diff` to `+max_diff`).
    """
    # Check shape compatibility
    if len(MODEL_NAMES) != len(PIPE):
        raise ValueError("MODEL_NAMES and PIPE must have the same length")

    # Standardize string inputs to list
    if isinstance(key_list, str):
        key_list = [key_list]

    # Enforce generic fallback nomenclature if duplicate names are provided
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)):
        MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    # --- INTERNAL HELPERS ---
    def title_mol(mol_name, frac, path, verbose):
        """Formats and returns a LaTeX string for chemical species titles and grain environments."""
        m = re.match(r"^([JK])(\d+)", mol_name)
        env = (
            f"{'surface' if m.group(1) == 'J' else 'mantle'} at grain size = {get_grain_size_in_um(Path(path)/'1D_grain_sizes.in', m.group(2), verbose=verbose)} µm"
            if m else "none"
        )
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        
        if env != "none": 
            return f"${f}$ [$n_{{{f}}}/n_H$]\n({env})" if frac else f"${f}$ [$n_{{{f}}}$] [cm$^{{-3}}$]\n({env})"
        else:
            return f"${f}$ [$n_{{{f}}}/n_H$]" if frac else f"${f}$ [$n_{{{f}}}$] [cm$^{{-3}}$]"

    def title_phys(variable):
        """Formats physical variable names and appends units based on string keywords."""
        name = variable.replace("_", " ").title()
        if "temperature" in name.lower(): 
            name += " [K]"
        elif "extinction" in name.lower(): 
            name += " [mag]"
        elif "density" in name.lower(): 
            name += " [$cm^{-3}$]"
        return name

    def get_grain_size_in_um(file_path, bin_index, verbose=False):
        """Parses 1D_grain_sizes.in to retrieve the grain bin radius mapped to micrometers."""
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('!'):
                        continue
                    if '!' in line:
                        line = line.split('!')[0].strip()
                    values = [float(val) for val in line.split()]
                    if not values:
                        continue
                    num_grains = len(values) // 4
                    radii_cm = values[:num_grains]
                    index = int(bin_index) - 1
                    if 0 <= index < num_grains:
                        return radii_cm[index] * 10000.0
            return None
        except FileNotFoundError:
            return None
    
    # --- DATA EXTRACTION & MESH GENERATION ---
    model_data = {}         # Stores structured spatial data per model
    all_global_values = []  # Aggregates values for global color scale calculations
    
    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', f"{MODEL_NAMES[p_idx]}")
        main_output_dict = p.chemistry
        chempath = Path(p.chempath)
        species_data = {k: [] for k in key_list}
        
        # Loop over spatial points (radial grid keys)
        for r_value in main_output_dict.keys():
            folder_name = f"{r_value}AU"
            file_path = os.path.join(chempath, folder_name, "1D_static.dat")
            
            if os.path.exists(file_path):
                try:
                    z_points = np.loadtxt(file_path, comments='!', usecols=0)
                    sub_dict = main_output_dict[r_value]
                    
                    for key in key_list:
                        if MODE == 'physical':
                            full_array = sub_dict[key]
                            v_points = full_array[itime, :].copy()
                        elif MODE == 'chemistry':
                            abundance_array = sub_dict['abundances']
                            v_points = abundance_array.isel(time=itime).sel(species=key).values.copy()
                            if not fracab:
                                nH = sub_dict["H_number_density"][itime, :]
                                v_points *= nH  # Convert fractional to absolute density
                        
                        if len(z_points) == len(v_points):
                            species_data[key].append({
                                'R': float(r_value),  
                                'z': np.array(z_points),
                                'v': np.array(v_points)
                            })
                except Exception as e:
                    if verbose: print(f"Error processing {key} for R={r_value}: {e}")
            else:
                if verbose: print(f"File not found: {file_path}")

        # Reconstruct 2D mesh layout using intermediate edge boundaries
        plot_structures = {}
        for key in key_list:
            columns_data = sorted(species_data[key], key=lambda x: x['R'])
            if not columns_data: continue
                
            polygons = []
            values = []
            radii = [col['R'] for col in columns_data]
            
            # Compute radial boundaries
            if len(radii) > 1:
                r_midshifts = 0.5 * np.diff(radii)
                r_edges = [radii[0] - r_midshifts[0]] + [radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [radii[-1] + r_midshifts[-1]]
            else:
                r_edges = [radii[0] - 0.5, radii[0] + 0.5]
                
            # Compute vertical cell boundaries and construct polygons
            for i, col in enumerate(columns_data):
                r_left, r_right = r_edges[i], r_edges[i+1]
                z_pts, v_pts = col['z'], col['v']
                z_midshifts = 0.5 * np.diff(z_pts)
                z_edges = [z_pts[0] - z_midshifts[0]] + [z_pts[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_pts[-1] + z_midshifts[-1])]
                
                for j in range(len(v_pts)):
                    poly = [(r_left, z_edges[j]), (r_right, z_edges[j]), (r_right, z_edges[j+1]), (r_left, z_edges[j+1])]
                    polygons.append(poly)
                    values.append(v_pts[j])
                    
            vals_array = np.array(values)
            plot_structures[key] = {
                'polygons': polygons, 
                'values': vals_array, 
                'radii': radii, 
                'all_z': np.concatenate([c['z'] for c in columns_data])
            }
            all_global_values.extend(values)
        
        model_data[p_name] = plot_structures
    
    # --- SUBPLOT CANVAS GEOMETRY CONFIGURATION ---
    num_models = len(PIPE)
    model_names = list(model_data.keys())
    
    if num_models == 1:
        cols = min(3, len(key_list))
        rows = (len(key_list) + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4.5 * rows), squeeze=False)
        axes = axes.flatten()
    elif num_models == 2:
        cols = 3  # Columns: [Model 1] [Model 2] [Residuals]
        rows = len(key_list)
        fig, axes = plt.subplots(rows, cols, figsize=(15, 4.5 * rows), squeeze=False)
    else:
        cols = num_models
        rows = len(key_list)
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4.5 * rows), squeeze=False)

    # --- UNIFIED GLOBAL SCALE CALCULATION ---
    has_log_anywhere = MODE == 'chemistry' or any("density" in k.lower() or "extinction" in k.lower() for k in key_list)
    
    if common_scale and all_global_values:
        all_global_values = np.array(all_global_values)
        if has_log_anywhere:
            pos_vals = all_global_values[all_global_values > 0]
            global_vmin = vmin if vmin is not None else (max(1e-15, pos_vals.min()) if len(pos_vals) > 0 else 1e-15)
        else:
            global_vmin = vmin if vmin is not None else all_global_values.min()
        global_vmax = vmax if vmax is not None else all_global_values.max()
        
    # --- RENDERING ENGINE LOOP ---
    for row_idx, key in enumerate(key_list):
        is_log = MODE == 'chemistry' or "density" in key.lower() or "extinction" in key.lower()

        # 1. Render Base Model Panels
        for col_idx, p_name in enumerate(model_names):
            if key not in model_data[p_name]: continue
            
            ax = axes[row_idx, col_idx] if num_models > 1 else axes[row_idx]
            struct = model_data[p_name][key]
            vals = struct['values']
            
            # Determine localized or global color boundaries
            if common_scale:
                actual_vmin, actual_vmax = global_vmin, global_vmax
            else:
                if is_log:
                    pos_vals = vals[vals > 0]
                    actual_vmin = vmin if vmin is not None else (max(1e-15, pos_vals.min()) if len(pos_vals) > 0 else 1e-15)
                else:
                    actual_vmin = vmin if vmin is not None else vals.min()
                actual_vmax = vmax if vmax is not None else vals.max()
            
            color_norm = plt.cm.colors.LogNorm(vmin=actual_vmin, vmax=actual_vmax) if is_log else plt.cm.colors.Normalize(vmin=actual_vmin, vmax=actual_vmax)
            
            # Add polygonal 2D mesh to axes
            coll = PolyCollection(struct['polygons'], array=vals, cmap=colormap, norm=color_norm, edgecolors='none')
            ax.add_collection(coll)
            
            # Format titles and labels
            if MODE == 'physical': lab = title_phys(key)
            elif MODE == 'chemistry': lab = title_mol(key, fracab, Path(PIPE[col_idx].chempath) / f"{int(struct['radii'][0])}AU", verbose=verbose)
            
            sm = plt.cm.ScalarMappable(cmap=colormap, norm=color_norm)
            sm.set_array(vals)
            fig.colorbar(sm, ax=ax, label=lab)
            
            ax.set_title(f"{p_name}\n{lab}")
            ax.set_xlabel('R [AU]')
            ax.set_ylabel('z [AU]')
            ax.set_xlim(xlim if xlim is not None else (0, max(struct['radii']) * 1.02))
            ax.set_ylim(ylim if ylim is not None else (0, max(struct['all_z']) * 1.07))
    
        # 2. Render Residual Panels (Only triggered if num_models == 2)
        if num_models == 2:
            ax_res = axes[row_idx, 2]
            struct1 = model_data[model_names[0]][key]
            struct2 = model_data[model_names[1]][key]
            
            # Arithmetic linear discrepancy (Model 1 - Model 2)
            res_vals = struct1['values'] - struct2['values']
            
            # Enforce symmetric zero-centered bounds
            max_diff = max(abs(res_vals.min()), abs(res_vals.max()))
            if max_diff == 0: max_diff = 1.0
            res_norm = plt.cm.colors.Normalize(vmin=-max_diff, vmax=max_diff)
            
            coll_res = PolyCollection(struct1['polygons'], array=res_vals, cmap=res_colormap, norm=res_norm, edgecolors='none')
            ax_res.add_collection(coll_res)
            
            sm_res = plt.cm.ScalarMappable(cmap=res_colormap, norm=res_norm)
            sm_res.set_array(res_vals)
            fig.colorbar(sm_res, ax=ax_res, label=f"Difference ({model_names[0]} - {model_names[1]})")
            
            ax_res.set_title("Residuals")
            ax_res.set_xlabel('R [AU]')
            ax_res.set_ylabel('z [AU]')
            ax_res.set_xlim(xlim if xlim is not None else (0, max(struct1['radii']) * 1.02))
            ax_res.set_ylim(ylim if ylim is not None else (0, max(struct1['all_z']) * 1.07))
            
    # --- FIG SUPTITLE (EXTRACT METADATA TIME) ---
    try:
        first_p = list(model_data.keys())[0]
        first_key = list(model_data[first_p].keys())[0]
        first_r = model_data[first_p][first_key]['radii'][0]
        time_seconds = PIPE[0].chemistry[first_r]['abundances'].coords['time'].values[itime]
        fig.suptitle(f'Simulation Comparison — $t = {time_seconds/3.156e7:.0f}$ years', fontsize=14, y=0.99)
    except:
        pass
        
    plt.tight_layout()
    plt.show()

import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def plot_midplane_nautilus_multi_comparison(PIPE,
                                           MODEL_NAMES,
                                           itime=-1,
                                           MODE='chemistry',
                                           key_list=['CO'],
                                           fracab=True,
                                           verbose=True,
                                           xlim=None,
                                           ylim=None,
                                           colormap="turbo",
                                           vmin=None,
                                           vmax=None,
                                           figsize=None):
    r"""
    Plots and compares 1D radial profiles of multiple variables or chemical species 
    strictly at the disk midplane (z = 0) across multiple simulation models.

    Colors are assigned to differentiate distinct keys/species, while unique 
    marker styles are mapping individual models to maintain graphical clarity.

    Parameters
    ----------
    PIPE : list
        Collection of model pipe objects to analyze. Each object must feature 
        `.chemistry` and `.chempath` attributes.
    MODEL_NAMES : list of str
        Display names assigned to the models for plot headers and legend mapping. 
        If duplicates are found, generic nomenclature is automatically generated.
    itime : int, optional
        Simulation timestep index to extract and visualize. Defaults to -1.
    MODE : {'chemistry', 'physical'}, optional
        The nature of the variables to extract and map ('chemistry' or 'physical').
    key_list : str or list of str, optional
        A single identifier string or list of target keys/species to plot.
    fracab : bool, optional
        If True, plots fractional abundances ($n_X/n_H$). If False, plots 
        absolute number densities ($cm^{-3}$).
    verbose : bool, optional
        If True, prints structural diagnostics and missing file reports.
    xlim : tuple of float, optional
        Custom manual boundaries (min, max) to enforce on the radial axis R [AU].
    ylim : tuple of float, optional
        Custom manual boundaries (min, max) to enforce on the vertical axis.
    colormap : str, optional
        Matplotlib colormap used to split colors across keys. Defaults to "turbo".
    vmin, vmax : float, optional
        Forced boundary constraints for the vertical profile sampling scale.
    figsize : tuple, optional
        Figure size

    Returns
    -------
    None
        Renders a single unified Matplotlib line canvas.
    """
    # Verify sequence lengths match
    if len(MODEL_NAMES) != len(PIPE):
        raise ValueError("MODEL_NAMES and PIPE must have the same length")
        
    # Standardize string inputs to list
    if isinstance(key_list, str):
        key_list = [key_list]

    # Handle duplicate model names gracefully
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)):
        MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    # --- INTERNAL HELPERS ---
    def title_mol(mol_name, path, verbose):
        """Formats full molecular LaTeX labels including grain environment details."""
        m = re.match(r"^([JK])(\d+)", mol_name)
        env = (
            f"{'surface' if m.group(1) == 'J' else 'mantle'} at {get_grain_size_in_um(Path(path)/'1D_grain_sizes.in', m.group(2), verbose=verbose)} µm"
            if m else "none"
        )
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        return f"${f}$ ({env})" if env != "none" else f"${f}$"

    def clean_molec(mol_name):
        """Strips grain codes and applies standard LaTeX subscript/superscript formatting."""
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        return f"${f}$"

    def title_phys(variable):
        """Converts snake_case physical keys to titlecase and attaches specific scientific units."""
        name = variable.replace("_", " ").title()
        if "temperature" in name.lower(): name += " [K]"
        elif "extinction" in name.lower(): name += " [mag]"
        elif "density" in name.lower(): name += " [$cm^{-3}$]"
        return name

    def get_grain_size_in_um(file_path, bin_index, verbose=False):
        """Retrieves and translates micro-meter size thresholds from grain profile config lines."""
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('!'): continue
                    if '!' in line: line = line.split('!')[0].strip()
                    values = [float(val) for val in line.split()]
                    if not values: continue
                    num_grains = len(values) // 4
                    radii_cm = values[:num_grains]
                    index = int(bin_index) - 1
                    if 0 <= index < num_grains: return radii_cm[index] * 10000.0
            return None
        except FileNotFoundError:
            return None

    # --- PLOT INITIALIZATION ---
    fig, ax = plt.subplots(figsize=(10, 6)) if figsize == None else plt.subplots(figsize=figsize)
    
    # Configure the categorical color selection profile
    if len(key_list) == 1:
        colors = [plt.colormaps[colormap](0.5)]
    else:
        colors = plt.colormaps[colormap](np.linspace(0, 0.9, len(key_list)))

    # Cycle configurations for tracing distinct models visually
    marker_pool = ['o', 's', '^', 'D', 'v', 'p', '*', 'X', 'h']
    linestyle_pool = ['-', '--', ':', '-.']

    all_valid_mins = []
    all_valid_maxs = []
    legend_labels = {}
    time_years_string = None

    # --- CORE EXTRACTION ENGINE ---
    # Outer Loop: unique color assigned per mapped variable/species
    for k_idx, key in enumerate(key_list):
        
        # Inner Loop: unique lines/markers assigned per simulation model
        for p_idx, p in enumerate(PIPE):
            p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
            main_output_dict = p.chemistry
            
            radii_list = []
            values_list = []

            # Extract data from the midplane index (last grid cell of each column)
            for r_value in main_output_dict.keys():
                try:
                    sub_dict = main_output_dict[r_value]
                    MIDPLANE_INDEX = -1 
                    
                    if MODE == 'physical':
                        v_midplane = sub_dict[key][itime, MIDPLANE_INDEX]
                    elif MODE == 'chemistry':
                        abundance_array = sub_dict['abundances']
                        v_midplane = float(abundance_array.isel(time=itime).sel(species=key).values[MIDPLANE_INDEX])
                        if not fracab:
                            v_midplane *= sub_dict["H_number_density"][itime, MIDPLANE_INDEX]
                    
                    radii_list.append(float(r_value))
                    values_list.append(v_midplane)
                    
                    # Capture time tracking string metadata safely once
                    if time_years_string is None:
                        try:
                            t_sec = sub_dict['abundances'].coords['time'].values[itime]
                            time_years_string = f"{t_sec / 3.156e7:.0f}"
                        except: pass

                except Exception as e:
                    if verbose:
                        print(f"[{p_name}] Error processing R={r_value}, Key={key}: {e}")

            if not radii_list:
                continue

            # Ensure data arrays are strictly ordered by increasing radius
            sort_indices = np.argsort(radii_list)
            radii_arr = np.array(radii_list)[sort_indices]
            values_arr = np.array(values_list)[sort_indices]

            # Filter out non-positive data boundaries for robust logarithmic evaluation
            pos_values = values_arr[values_arr > 0]
            if len(pos_values) > 0:
                all_valid_mins.append(pos_values.min())
            all_valid_maxs.append(values_arr.max())

            # Synthesize human-readable string labels
            if MODE == 'physical':
                clean_label = title_phys(key)
            else:
                clean_label = clean_molec(key)
                
            if len(key_list) == 1:
                full_label = f"{p_name}"
            else:
                full_label = f"{clean_label} ({p_name})"

            if k_idx == 0:
                legend_labels[p_idx] = clean_label

            marker_style = marker_pool[p_idx % len(marker_pool)]
            line_style = linestyle_pool[p_idx % len(linestyle_pool)]

            # Generate the specific dataset line trace on the shared canvas
            ax.plot(radii_arr, values_arr, 
                    color=colors[k_idx], 
                    linestyle=line_style, 
                    marker=marker_style, 
                    markersize=5, 
                    linewidth=1.5, 
                    label=full_label)

    if not all_valid_maxs:
        if verbose: print("No plottable datasets parsed across pipelines.")
        return

    # --- CANVAS POST-PROCESSING & NORMS ---
    is_log = MODE == 'chemistry' or any("density" in k.lower() or "extinction" in k.lower() for k in key_list)
    
    if is_log:
        ax.set_yscale('log')
        global_min = min(all_valid_mins) if all_valid_mins else 1e-15
        global_max = max(all_valid_maxs) if all_valid_maxs else 1.0
        actual_vmin = vmin if vmin is not None else max(1e-15, global_min)
        actual_vmax = vmax if vmax is not None else global_max * 2
        ax.set_ylim(actual_vmin, actual_vmax)
    else:
        if vmin is not None or vmax is not None:
            ax.set_ylim(vmin, vmax)

    ax.set_xlabel('Radius R [AU]')
    
    # Context-dependent vertical axis label population
    if len(key_list) == 1:
        if MODE == 'chemistry':
            ax.set_ylabel(f"{clean_label} — " + ("Fractional Abundance [$n_X/n_H$]" if fracab else "Number Density [$cm^{-3}$]"))
        else:
            ax.set_ylabel(clean_label)
    else:
        if MODE == 'chemistry':
            ax.set_ylabel("Fractional Abundance [$n_X/n_H$]" if fracab else "Number Density [$cm^{-3}$]")
        else:
            ax.set_ylabel("Physical Units (See Legend)")

    # Set canvas descriptive header
    if time_years_string:
        ax.set_title(f'Midplane ($z = 0$) Radial Profile Comparison — $t = {time_years_string}$ years')
    else:
        ax.set_title('Midplane ($z = 0$) Radial Profile Comparison')

    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)

    ax.legend(loc='best', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.5)
    
    plt.tight_layout()
    plt.show()


import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def plot_grain_surface_midplane_comparison(PIPE,
                                          MODEL_NAMES,
                                          itime=-1,
                                          verbose=True,
                                          xlim=None,
                                          ylim=None,
                                          colormap="viridis"):
    r"""
    Plots and compares the total available grain surface area strictly at the 
    disk midplane (z = 0) as a function of Radius across multiple models.

    The total physical grain surface area per unit volume ($cm^2/cm^3$) is calculated 
    by performing a integrated vectorized sum across all available size bins:
    $$ \Sigma_{\text{grain}} = 4\pi \cdot n_H \cdot \sum_{i=1}^{N} \frac{a_i^2}{\text{GTODN}_i} $$

    Parameters
    ----------
    PIPE : list
        Collection of model pipe objects to analyze. Each object must contain 
        `.chemistry` (main output dictionary) and `.chempath` attributes.
    MODEL_NAMES : list of str
        Display names assigned to the models for the plot legend. If the list 
        contains duplicate names, generic fallbacks are automatically generated.
    itime : int, optional
        Simulation timestep index to extract and visualize.
        Defaults to -1 (the final recorded timestep).
    verbose : bool, optional
        If True, prints diagnostic warnings regarding missing configuration files 
        ('1D_grain_sizes.in') or calculation exceptions. Defaults to True.
    xlim : tuple of float, optional
        Custom manual boundaries (min, max) to enforce on the radial axis R [AU].
    ylim : tuple of float, optional
        Custom manual boundaries (min, max) to enforce on the vertical axis [cm^2/cm^3].
    colormap : str, optional
        The name of the Matplotlib colormap used to automatically differentiate 
        the lines if multiple models are passed. Defaults to "viridis".

    Returns
    -------
    None
        Renders a single combined Matplotlib line plot canvas.
    """
    # Verify input list matching sizes
    if len(MODEL_NAMES) != len(PIPE):
        raise ValueError("MODEL_NAMES and PIPE must have the same length")
        
    # Generate distinct generic model names if duplicates exist
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)):
        MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    # Setup the single figure canvas with log-scaled vertical axis
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.set_yscale('log')
    
    # Generate a discrete sampling sequence from the colormap
    colors = plt.colormaps[colormap].resampled(len(PIPE))(np.linspace(0, 0.95, len(PIPE)))
    
    global_max_val = -np.inf
    global_min_pos_val = np.inf
    time_years_string = None

    # --- LOOP OVER ALL MODELS IN PIPE ---
    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        main_output_dict = p.chemistry
        chempath = Path(p.chempath)
        
        radii_list = []
        surface_list = []

        # Loop over spatial radius keys for the current model
        for r_value in main_output_dict.keys():
            folder_name = f"{r_value}AU"
            file_path = os.path.join(chempath, folder_name, "1D_grain_sizes.in")
            
            if os.path.exists(file_path):
                try:
                    # 1. Load data from 1D_grain_sizes.in (ignoring comments)
                    grain_data = np.loadtxt(file_path, comments='!')
                    
                    # Extract only the midplane layer (the last row, index -1)
                    midplane_row = grain_data[-1, :]
                    
                    # 2. Determine the number of grain bins (N)
                    total_columns = len(midplane_row)
                    N = int(total_columns / 4)
                    
                    # 3. Extract group 1 (grain radii 'a') and group 2 ('GTODN')
                    a_array = midplane_row[0:N]          # First N columns [cm]
                    gtodn_array = midplane_row[N:2*N]    # Next N columns
                    
                    # 4. Extract Hydrogen number density (nH) at the midplane
                    sub_dict = main_output_dict[r_value]
                    nH_midplane = sub_dict["H_number_density"][itime, -1]
                    
                    # 5. Compute the physical formula (Vectorized sum over all bins)
                    total_surface = 4 * np.pi * nH_midplane * np.sum((a_array**2) / gtodn_array)
                    
                    radii_list.append(float(r_value))
                    surface_list.append(total_surface)
                    
                    # Safely grab simulation time for subtitle from the first valid entry
                    if time_years_string is None:
                        try:
                            t_sec = sub_dict['abundances'].coords['time'].values[itime]
                            time_years_string = f"{t_sec / 3.156e7:.0f}"
                        except:
                            pass
                            
                except Exception as e:
                    if verbose: 
                        print(f"[{p_name}] Error processing grain data for R={r_value}: {e}")
            else:
                if verbose: 
                    print(f"[{p_name}] File not found: {file_path}")

        # If data was successfully collected for this model, sort and plot it
        if radii_list:
            sort_indices = np.argsort(radii_list)
            radii_arr = np.array(radii_list)[sort_indices]
            surface_arr = np.array(surface_list)[sort_indices]
            
            # Update global scale trackers for safe auto-scaling later
            global_max_val = max(global_max_val, surface_arr.max())
            pos_vals = surface_arr[surface_arr > 0]
            if len(pos_vals) > 0:
                global_min_pos_val = min(global_min_pos_val, pos_vals.min())

            # Plot model data on the shared axis
            ax.plot(radii_arr, surface_arr, 
                    color=colors[p_idx], 
                    linestyle='-', 
                    marker='s', 
                    markersize=4, 
                    linewidth=1.5, 
                    label=p_name)
        else:
            if verbose:
                print(f"[{p_name}] No data collected. Skipped from plotting canvas.")

    # --- CANVAS POST-PROCESSING & FORMATTING ---
    ax.set_xlabel('Radius R [AU]')
    ax.set_ylabel(r'Total Grain Surface Area [$\text{cm}^{2}/\text{cm}^{3}$]')
    
    # Title formatting based on retrieved time metadata
    if time_years_string:
        ax.set_title(r'Total Grain Surface Area at Midplane ($z=0$) - $t = ' + time_years_string + r'$ yr')
    else:
        ax.set_title(r'Total Grain Surface Area at Midplane ($z=0$)')

    # Apply boundary limits (User forced vs Dynamic safe bounds)
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: 
        ax.set_ylim(ylim)
    else:
        # Dynamic log safe fallback bounding
        actual_vmin = max(1e-25, global_min_pos_val if global_min_pos_val != np.inf else 1e-25)
        actual_vmax = global_max_val * 2 if global_max_val != -np.inf else 1e-10
        ax.set_ylim(actual_vmin, actual_vmax)

    ax.grid(True, linestyle=':', alpha=0.5)
    ax.legend(loc='best', frameon=True)
    
    plt.tight_layout()
    plt.show()


import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_vertical_cut_nautilus_comparison(PIPE,
                                         MODEL_NAMES,
                                         R,
                                         species='CO',
                                         itime=-1,
                                         fracab=True,
                                         colormap="turbo",
                                         xlim=None,
                                         ylim=None,
                                         xscale="linear",
                                         yscale="linear",
                                         verbose=True):
    r"""
    Plots and compares 1D vertical cut profiles (abundance or density vs. height z) 
    across multiple simulation models using NAUTILUS data outputs.

    Supports tracing either multiple radii for a single species, or multiple species 
    at a single radius. Distinct targets are mapped to colors, whereas individual 
    simulation models are distinguished using unique markers and linestyles.

    Parameters
    ----------
    PIPE : list
        Collection of model pipe objects to analyze. Each object must feature 
        `.chemistry` and `.chempath` attributes.
    MODEL_NAMES : list of str
        Display names assigned to the models for legends. If duplicate names 
        are found, a generic fallback nomenclature is generated.
    R : int, float or list
        Target radius or list of radii (in AU) to extract data for.
    species : str or list of str, optional
        Target species formula string or list of species to plot. Defaults to 'CO'.
    itime : int, optional
        Simulation timestep index to visualize. Defaults to -1 (final timestep).
    fracab : bool, optional
        If True, plots fractional abundances ($n_X/n_H$). If False, plots 
        absolute number densities ($cm^{-3}$). Defaults to True.
    colormap : str, optional
        Matplotlib colormap string used to split colors across profile keys. Defaults to "turbo".
    xlim, ylim : tuple of float, optional
        Custom manual boundaries to enforce on the horizontal and vertical axes.
    xscale, yscale : {'linear', 'log'}, optional
        Matplotlib axis scale configurations. Defaults to 'linear'.
    verbose : bool, optional
        If True, prints missing file or parsing diagnostic messages.

    Returns
    -------
    None
        Renders a single unified Matplotlib line canvas.
    """
    chempath_ref = Path(PIPE[0].chempath)

    # Verify input sequence lengths match
    if len(MODEL_NAMES) != len(PIPE):
        raise ValueError("MODEL_NAMES and PIPE must have the same length")

    # Standardize inputs to lists for uniform loop processing
    r_list = [R] if not isinstance(R, list) else R
    species_list = [species] if not isinstance(species, list) else species
    r_list = [int(r) for r in r_list]

    # Enforce mutual exclusivity restriction between radii and species lists
    if len(r_list) > 1 and len(species_list) > 1:
        raise ValueError("Cannot supply multiple values for both 'R' and 'species' simultaneously. One parameter must be a list of length 1.")

    # Prevent layout legend collisions by forcing unique fallback nomenclature
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)):
        MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    # --- INTERNAL HELPERS ---
    def title_mol(mol_name, path, verbose):
        """Formats chemical labels with embedded physical grain size mapping metadata."""
        m = re.match(r"^([JK])(\d+)", mol_name)
        env = (
            f"{'surface' if m.group(1) == 'J' else 'mantle'} at {get_grain_size_in_um(Path(path)/'1D_grain_sizes.in', m.group(2), verbose=verbose)} µm"
            if m else "none"
        )
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        return f"${f}$ ({env})" if env != "none" else f"${f}$"

    def clean_molec(mol_name):
        """Strips numerical grain phases to generate standardized chemical LaTeX subtitles."""
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        return f"${f}$"

    def get_grain_size_in_um(file_path, bin_index, verbose=False):
        """Parses configuration file lines to determine mapped grain radius configurations."""
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('!'): continue
                    if '!' in line: line = line.split('!')[0].strip()
                    values = [float(val) for val in line.split()]
                    if not values: continue
                    num_grains = len(values) // 4
                    radii_cm = values[:num_grains]
                    index = int(bin_index) - 1
                    if 0 <= index < num_grains: return radii_cm[index] * 10000.0
            return None
        except FileNotFoundError:
            return None

    # --- PLOT INITIALIZATION ---
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Calculate colors based on total properties traced per model
    total_distinct_targets = max(len(r_list), len(species_list))
    if total_distinct_targets == 1:
        colors = [plt.colormaps[colormap](0.5)]
    else:
        colors = plt.colormaps[colormap](np.linspace(0, 0.9, total_distinct_targets))

    # Cyclic layout pools to differentiate distinct simulation data pipelines
    marker_pool = ['o', 's', '^', 'D', 'v', 'p', '*', 'X']
    linestyle_pool = ['-', '--', ':', '-.']

    clean_label_ref = ""
    total_successful_plots = 0
    time_years_string = None

    # --- CORE EXTRACTION ENGINE ---
    # Outer Loop: Targets (Radii or Species) -> Mapped to distinct Colors
    target_idx = 0
    for r_val in r_list:
        for spec_val in species_list:
            
            # Inner Loop: Models -> Mapped to distinct Styles (Linestyles / Markers)
            for p_idx, p in enumerate(PIPE):
                p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
                main_output_dict = p.chemistry
                chempath = Path(p.chempath)
                
                if r_val not in main_output_dict:
                    if verbose: print(f"[{p_name}] Radius {r_val} AU missing inside data structure. Skipping.")
                    continue

                try:
                    # 1. Load the localized vertical grid coordinates (z)
                    static_file = chempath / f"{r_val}AU" / "1D_static.dat"
                    static = pd.read_table(static_file, sep=r'\s+', comment='!', header=None, engine='python')
                    z = static[0].values  
                    
                    # 2. Extract abundance values array
                    ab = main_output_dict[r_val]['abundances']
                    sp_arr = ab.isel(time=itime).sel(species=spec_val).values
                    
                    # 3. Apply profile conversion factor if looking for absolute density
                    if not fracab:
                        nH = main_output_dict[r_val]['H_number_density'][itime, :]
                        n_plot = nH * sp_arr
                    else:
                        n_plot = sp_arr

                    # 4. Label generation strings
                    label_mol_text = title_mol(spec_val, chempath / f"{int(r_val)}AU", verbose=verbose)
                    clean_label_ref = clean_molec(spec_val)
                    
                    if len(PIPE) > 1:
                        if len(r_list) > 1:
                            legend_string = f"{label_mol_text} @ {r_val} AU ({p_name})"
                        elif len(species_list) > 1:
                            legend_string = f"{label_mol_text} ({p_name})"
                        else:
                            legend_string = f"{p_name}"
                    else:
                        if len(r_list) > 1:
                            legend_string = f"{label_mol_text} @ {r_val} AU"
                        elif len(species_list) > 1:
                            legend_string = f"{label_mol_text}"
                        else:
                            legend_string = f"{p_name}"

                    # Retrieve simulation time metadata safely once
                    if time_years_string is None:
                        try:
                            t_sec = ab.coords['time'].values[itime]
                            time_years_string = f"{t_sec / 3.156e7:.0f}"
                        except: pass

                    # 5. Extract unique styles for the model row
                    marker_style = marker_pool[p_idx % len(marker_pool)]
                    line_style = linestyle_pool[p_idx % len(linestyle_pool)]

                    # 6. Render step plot trace on shared canvas axes
                    ax.scatter(n_plot, z, color=colors[target_idx], s=15, marker=marker_style)
                    ax.plot(n_plot, z, 
                            color=colors[target_idx], 
                            linestyle=line_style, 
                            linewidth=1.5, 
                            label=legend_string)
                    
                    total_successful_plots += 1

                except Exception as e:
                    if verbose: 
                        print(f"[{p_name}] Error parsing vertical cut for R={r_val} AU, species={spec_val}: {e}")
            
            # Increment color index per unique species/radius step
            total_distinct_plots_per_model = max(len(r_list), len(species_list))
            if total_distinct_plots_per_model > 1:
                target_idx += 1

    if total_successful_plots == 0:
        if verbose: print("No profile data successfully parsed across paths.")
        return

    # --- CANVAS POST-PROCESSING & LABELS ---
    ax.set_ylabel("z [AU]")
    
    if len(species_list) == 1 and len(r_list) == 1:
        ax.set_xlabel(f"{clean_label_ref} — " + ("Fractional Abundance [$n_X/n_H$]" if fracab else "Number Density [$cm^{-3}$]"))
    else:
        ax.set_xlabel("Fractional Abundance [$n_X/n_H$]" if fracab else "Number Density [$cm^{-3}$]")

    # Context-aware titles structure formatting
    if time_years_string:
        if len(PIPE) == 1:
            if len(r_list) == 1:
                ax.set_title(f"Vertical Cut Profile Comparison @ R = {r_list[0]} AU — $t = {time_years_string}$ years \n{MODEL_NAMES[0]}")
            else:
                ax.set_title(f"Multi-Radius Vertical Cut Profile Comparison — $t = {time_years_string}$ years\n{MODEL_NAMES[0]}")
        else:
            if len(r_list) == 1:
                ax.set_title(f"Vertical Cut Profile Comparison @ R = {r_list[0]} AU — $t = {time_years_string}$ years")
            else:
                ax.set_title(f"Multi-Radius Vertical Cut Profile Comparison — $t = {time_years_string}$ years")
    else:
        ax.set_title("Vertical Cut Profile Comparison")

    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    
    ax.legend(loc='best', frameon=True, ncol=len(species) if len(PIPE) > 2 else 2)
    ax.grid(True, linestyle=':', alpha=0.5)
    
    plt.tight_layout()
    plt.show()


import os
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

def plot_atom_ratio_nautilus_comparison(PIPE,
                                       MODEL_NAMES,
                                       ratio_list=['C/O'],
                                       itime=-1,
                                       verbose=True,
                                       xlim=None,
                                       ylim=None,
                                       colormap="gnuplot",
                                       res_colormap="seismic",
                                       common_scale=True):
    r"""
    Plots 2D vertical cross-sections (poloidal cuts) of elemental abundance ratios 
    for specified atoms, dynamically managing multi-model structural comparisons.

    Dynamically manages structural subplots depending on the number of models in PIPE:
    - 1 Model: Displays standard grid layout for selected ratios (max 3 columns).
    - 2 Models: Displays a 3-column layout (Model 1, Model 2, and Model 1 - Model 2 Residuals).
    - >2 Models: Displays a multi-column grid where each column represents a single model.

    Parameters
    ----------
    PIPE : list
        Collection of model pipe objects to analyze. Each object must feature 
        `.chemistry` and `.chempath` attributes.
    MODEL_NAMES : list of str
        Display names assigned to the models for plot headers. If duplicate names 
        are found, a generic fallback nomenclature is generated.
    ratio_list : str or list of str, optional
        Standalone string or list of targeted element division pairs (e.g., ['C/O', 'Mg/Si']).
        Defaults to ['C/O'].
    itime : int, optional
        Simulation timestep index to visualize. Defaults to -1 (final timestep).
    verbose : bool, optional
        If True, prints missing file, axis size mismatch, or null-denominator cell warnings.
    xlim, ylim : tuple of float, optional
        Custom (min, max) boundaries for the horizontal Radius and vertical Altitude axes.
    colormap : str, optional
        Matplotlib colormap string used for main models. Defaults to "gnuplot".
    res_colormap : str, optional
        Matplotlib colormap string used for the residual map panels. Defaults to "seismic".
    common_scale : bool, optional
        If True, shares an identical unscaled global colorbar framework across ALL 
        models for a given ratio target. Defaults to True.

    Returns
    -------
    None
        Renders a structured Matplotlib figure window.
    """
    # Check length compatibility between names and pipes
    if len(MODEL_NAMES) != len(PIPE):
        raise ValueError("MODEL_NAMES and PIPE must have the same length")
        
    # Standardize string inputs to a list
    if isinstance(ratio_list, str):
        ratio_list = [ratio_list]

    # Authorized chemical elements list
    elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    
    # Parse and validate elemental ratio pairs
    parsed_ratios = []
    for item in ratio_list:
        if '/' not in item:
            raise ValueError(f"Invalid ratio token entry: '{item}'. It must contain a '/' divider symbol.")
        s1, s2 = item.split('/')
        if s1 not in elements or s2 not in elements:
            raise ValueError(f"Invalid element pair requested in '{item}'. Must belong to: {elements}")
        parsed_ratios.append((s1, s2, item))

    # Fallback to unique generic model names if duplicates are found
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)):
        MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    # --- INTERNAL HELPERS ---
    def count_species_elements(species_name, element1, element2):
        """Uses regex to count the occurrences of element1 and element2 in a chemical formula."""
        if species_name == 'e-': return {element1: 0, element2: 0}
        formula = species_name.replace('c-', '').replace('l-', '')
        if formula.endswith('+') or formula.endswith('-'):
            formula = formula[:-1]
        pattern = re.compile(r'([A-Z][a-z]?)(-?\d*)')
        composition = {}
        for atom, n in pattern.findall(formula):
            count = int(n) if n else 1
            composition[atom] = composition.get(atom, 0) + count
        return {
            element1: composition.get(element1, 0),
            element2: composition.get(element2, 0)
        }

    def keep_gas_species_only(species):
        """Filters out grain surface and mantle species prefixes to isolate gas-phase items."""
        motif = re.compile(r'^[JK]\d{2}|^GRAIN')
        return [e for e in species if not motif.match(e)]

    # --- MULTI-MODEL DATA ACQUISITION & MESH GENERATION ---
    model_data = {}  
    atom_cache = {}  # Cache structure to optimize element counting calculations
    time_years_string = None

    # Loop over pipelines
    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        main_output_dict = p.chemistry
        chempath = Path(p.chempath)
        
        ratio_database = {token: [] for _, _, token in parsed_ratios}

        # Scan through 1D spatial points directories
        for r_value in main_output_dict.keys():
            folder_name = f"{r_value}AU"
            file_path = os.path.join(chempath, folder_name, "1D_static.dat")
            
            if os.path.exists(file_path):
                try:
                    z_points = np.loadtxt(file_path, comments='!', usecols=0)
                    sub_dict = main_output_dict[r_value]
                    abundance_array = sub_dict['abundances']
                    
                    # Isolate pure gas phase components
                    local_species_list = keep_gas_species_only(list(abundance_array.coords['species'].values))
                    sliced_abundances = abundance_array.isel(time=itime).sel(species=local_species_list).values
                    
                    # Extract timestamp tracking metadata once
                    if time_years_string is None:
                        try:
                            t_sec = abundance_array.coords['time'].values[itime]
                            time_years_string = f"{t_sec / 3.156e7:.0f}"
                        except: pass

                    # Process each targeted atomic ratio
                    for s1, s2, token in parsed_ratios:
                        s1_coeffs = []
                        s2_coeffs = []
                        
                        # Populate stoichiometric coefficient weights
                        for species in local_species_list:
                            cache_key = f"{species}_{s1}_{s2}"
                            if cache_key not in atom_cache:
                                counts = count_species_elements(species, s1, s2)
                                atom_cache[cache_key] = (counts[s1], counts[s2])
                            
                            c1, c2 = atom_cache[cache_key]
                            s1_coeffs.append(c1)
                            s2_coeffs.append(c2)
                        
                        s1_coeffs = np.array(s1_coeffs)[:, np.newaxis]
                        s2_coeffs = np.array(s2_coeffs)[:, np.newaxis]
                        
                        # Calculate total abundance per element
                        total_s1 = np.sum(sliced_abundances * s1_coeffs, axis=0)
                        total_s2 = np.sum(sliced_abundances * s2_coeffs, axis=0)

                        if verbose and total_s2.min() == 0:
                            zero_cells = np.where(total_s2 == 0)[0]
                            print(f"[{p_name} - {token} @ R={r_value} AU] Denominator {s2} null inside {len(zero_cells)} cells.")
                        
                        # Safe division handling to prevent NaN runtime crashes
                        with np.errstate(divide='ignore', invalid='ignore'):
                            v_points = np.where(total_s2 > 0, total_s1 / total_s2, 0.0)
                        
                        if len(z_points) == len(v_points):
                            ratio_database[token].append({
                                'R': float(r_value),
                                'z': np.array(z_points),
                                'v': np.array(v_points)
                            })
                except Exception as e:
                    if verbose: print(f"[{p_name}] Error compiling ratios for R={r_value}: {e}")
            else:
                if verbose: print(f"[{p_name}] File not found: {file_path}")

        # Construct spatial 2D grid polygons for the current model
        plot_structures = {}
        for s1, s2, token in parsed_ratios:
            columns_data = sorted(ratio_database[token], key=lambda x: x['R'])
            if not columns_data: continue
                
            polygons = []
            values = []
            radii = [col['R'] for col in columns_data]
            
            # Form radial boundaries
            if len(radii) > 1:
                r_midshifts = 0.5 * np.diff(radii)
                r_edges = [radii[0] - r_midshifts[0]] + [radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [radii[-1] + r_midshifts[-1]]
            else:
                r_edges = [radii[0] - 0.5, radii[0] + 0.5]
                
            # Form altitude boundaries and assemble discrete cell loops
            for i, col in enumerate(columns_data):
                r_left, r_right = r_edges[i], r_edges[i+1]
                z_pts, v_pts = col['z'], col['v']
                z_midshifts = 0.5 * np.diff(z_pts)
                z_edges = [z_pts[0] - z_midshifts[0]] + [z_pts[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_pts[-1] + z_midshifts[-1])]
                
                for j in range(len(v_pts)):
                    poly = [(r_left, z_edges[j]), (r_right, z_edges[j]), (r_right, z_edges[j+1]), (r_left, z_edges[j+1])]
                    polygons.append(poly)
                    values.append(v_pts[j])
                    
            plot_structures[token] = {
                'polygons': polygons,
                'values': np.array(values),
                'radii': radii,
                'all_z': np.concatenate([c['z'] for c in columns_data])
            }
        model_data[p_name] = plot_structures

    # --- CANVAS GEOMETRY LAYOUT ENGINE ---
    num_models = len(PIPE)
    model_names = list(model_data.keys())
    num_ratios = len(parsed_ratios)

    if num_models == 1:
        cols = min(3, num_ratios)
        rows = (num_ratios + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4.5 * rows), squeeze=False)
        axes = axes.flatten()
    elif num_models == 2:
        cols = 3  # Columns: [Model 1] [Model 2] [Residuals]
        rows = num_ratios
        fig, axes = plt.subplots(rows, cols, figsize=(15, 4.5 * rows), squeeze=False)
    else:
        cols = num_models
        rows = num_ratios
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4.5 * rows), squeeze=False)

    # --- RENDERING GENERATION ENGINE ---
    for row_idx, (_, _, token) in enumerate(parsed_ratios):
        
        # Pre-calculate uniform normalization boundaries if common scaling is enabled
        if common_scale:
            all_row_vals = []
            for p_name in model_names:
                if token in model_data[p_name]:
                    all_row_vals.extend(model_data[p_name][token]['values'])
            all_row_vals = np.array(all_row_vals) if all_row_vals else np.array([0.0, 1.0])
            row_pos_vals = all_row_vals[all_row_vals > 0]
            global_row_vmin = row_pos_vals.min() if len(row_pos_vals) > 0 else 1e-15
            global_row_vmax = all_row_vals.max()
            global_row_is_log = len(row_pos_vals) > 0 and (global_row_vmax / global_row_vmin) > 10.0

        # 1. Render Main Simulation Panels
        for col_idx, p_name in enumerate(model_names):
            if token not in model_data[p_name]: continue
            
            ax = axes[row_idx, col_idx] if num_models > 1 else axes[row_idx]
            struct = model_data[p_name][token]
            vals = struct['values']
            
            if common_scale:
                actual_vmin, actual_vmax = global_row_vmin, global_row_vmax
                is_log = global_row_is_log
            else:
                pos_vals = vals[vals > 0]
                actual_vmin = pos_vals.min() if len(pos_vals) > 0 else 1e-15
                actual_vmax = vals.max()
                is_log = len(pos_vals) > 0 and (actual_vmax / actual_vmin) > 10.0
                
            color_norm = plt.cm.colors.LogNorm(vmin=actual_vmin, vmax=actual_vmax) if is_log else plt.cm.colors.Normalize(vmin=actual_vmin, vmax=actual_vmax)
            
            # Map values to PolyCollection
            coll = PolyCollection(struct['polygons'], array=vals, cmap=colormap, norm=color_norm, edgecolors='none')
            ax.add_collection(coll)
            
            sm = plt.cm.ScalarMappable(cmap=colormap, norm=color_norm)
            sm.set_array(vals)
            fig.colorbar(sm, ax=ax, label=f"Ratio [{token}]")
            
            ax.set_title(f"{p_name}\nAtomic Ratio: {token}")
            ax.set_xlabel('R [AU]')
            ax.set_ylabel('z [AU]')
            ax.set_xlim(xlim if xlim is not None else (0, max(struct['radii']) * 1.02))
            ax.set_ylim(ylim if ylim is not None else (0, max(struct['all_z']) * 1.07))

        # 2. Render Residual Panels (Only triggered if num_models == 2)
        if num_models == 2:
            ax_res = axes[row_idx, 2]
            struct1 = model_data[model_names[0]][token]
            struct2 = model_data[model_names[1]][token]
            
            # Linear discrepancy layout calculations (Model 1 - Model 2)
            res_vals = struct1['values'] - struct2['values']
            max_diff = max(abs(res_vals.min()), abs(res_vals.max()))
            if max_diff == 0: max_diff = 1.0
            res_norm = plt.cm.colors.Normalize(vmin=-max_diff, vmax=max_diff)
            
            coll_res = PolyCollection(struct1['polygons'], array=res_vals, cmap=res_colormap, norm=res_norm, edgecolors='none')
            ax_res.add_collection(coll_res)
            
            sm_res = plt.cm.ScalarMappable(cmap=res_colormap, norm=res_norm)
            sm_res.set_array(res_vals)
            fig.colorbar(sm_res, ax=ax_res, label=f"Difference ({model_names[0]} - {model_names[1]})")
            
            ax_res.set_title(f"Residuals: {token}")
            ax_res.set_xlabel('R [AU]')
            ax_res.set_ylabel('z [AU]')
            ax_res.set_xlim(xlim if xlim is not None else (0, max(struct1['radii']) * 1.02))
            ax_res.set_ylim(ylim if ylim is not None else (0, max(struct1['all_z']) * 1.07))

    # Clean up unoccupied layout grids in single model plots
    if num_models == 1 and num_ratios > 1:
        for idx in range(num_ratios, len(axes)):
            fig.delaxes(axes[idx])

    # Global window suptitle management 
    if time_years_string:
        if num_models == 1 and num_ratios == 1:
            axes[0].set_title(f"{axes[0].get_title()} \n $t = {time_years_string}$ years")
        else:
            fig.suptitle(f'Gas Phase Elemental Ratios — $t = {time_years_string}$ years', fontsize=14, y=0.99)

    plt.tight_layout()
    plt.show()



def plot_top_contributing_species(chempath,
                                  main_output_dict,
                                  target_atom="C",
                                  itime=-1,
                                  verbose=True,
                                  spnumber=5,
                                  color="darkred",
                                  phase="gas",
                                  grain_bin=None):
    """
    Computes and plots the top chemical species contributing to the global, 
    volume-integrated total budget of a target chemical element within a protoplanetary disk.

    Parameters:
    -----------
    chempath : str
        Path to the directory containing spatial grid subfolders (e.g., '10AU/', '100AU/').
    main_output_dict : dict
        Nested dictionary mapping radial keys to data sub-structures containing 'abundances' 
        (xarray.DataArray) and 'H_number_density' spatial profiles.
    target_atom : str, default "C"
        The chemical element symbol whose structural reservoir distribution is evaluated.
    itime : int, default -1
        Time index to slice from the multi-epoch data arrays.
    verbose : bool, default True
        If True, outputs warning and missing file notifications to the console.
    spnumber : int, default 5
        Number of top contributing chemical species to display on the bar chart.
    color : str, default "darkred"
        Visual fill color of the generated bar elements.
    phase : str, default "gas"
        The chemical phase environment to isolate ("gas", "surface", "mantle", "grain", or "all").
    grain_bin : int or str, optional
        Specific grain size category to filter when analyzing surface or mantle ice matrices.
    """
                                    
    allowed_elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    if target_atom not in allowed_elements:
        raise ValueError(f"Target atom '{target_atom}' is not recognized in the chemical network.")
    
    valid_phases = ["gas", "surface", "mantle", "grain", "all"]
    if phase not in valid_phases:
        raise ValueError(f"Phase '{phase}' unrecognized. Choose among {valid_phases}")
    if grain_bin is not None and phase in ["gas","all"] : raise ValueError("grain_bin and gas phase can not be defined simultaneously")

    chempath = Path(chempath)

    # --- INTERNAL HELPERS ---
    def get_grain_size_in_um(file_path, bin_index):
        """Parses 1D_grain_sizes.in to retrieve the grain bin radius mapped to micrometers."""
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('!'):
                        continue
                    if '!' in line:
                        line = line.split('!')[0].strip()
                    values = [float(val) for val in line.split()]
                    if not values:
                        continue
                    num_grains = len(values) // 4
                    radii_cm = values[:num_grains]
                    index = int(bin_index) - 1
                    if 0 <= index < num_grains:
                        return radii_cm[index] * 10000.0
            return None
        except FileNotFoundError:
            return None

    def parse_species(species_name):
        if species_name == 'e-':
            return "gas", None, "e-"
            
        grain_match = re.match(r'^([JK])(\d+)(.+)', species_name)
        
        if grain_match:
            p_code, g_bin, raw_formula = grain_match.groups()
            sp_phase = "surface" if p_code == 'J' else "mantle"
            sp_bin = g_bin
        else:
            sp_phase = "gas"
            sp_bin = None
            raw_formula = species_name
        
        clean_formula = raw_formula.replace('c-', '').replace('l-', '')
        return sp_phase, sp_bin, clean_formula

    def count_target_atom(clean_formula, target):
        if clean_formula == 'e-': return 0
        
        calc_formula = clean_formula
        if calc_formula.endswith('+') or calc_formula.endswith('-'):
            calc_formula = calc_formula[:-1]
            
        calc_formula = calc_formula.replace('-', '')
            
        pattern = re.compile(r'([A-Z][a-z]?)(\d*)') 
        composition = {}
        for atom, n in pattern.findall(calc_formula):
            try:
                count = int(n) if n else 1
            except ValueError:
                count = 1
            composition[atom] = composition.get(atom, 0) + count
        return composition.get(target, 0)

    def filter_and_cache_species(species_list):
        valid_list = []
        coeffs = []
        
        for sp in species_list:
            if "GRAIN" in sp:
                continue
                
            sp_phase, sp_bin, clean_formula = parse_species(sp)
            
            if phase == "all":
                pass
            elif phase == "grain":
                if sp_phase not in ["surface", "mantle"]:
                    continue
            elif sp_phase != phase:
                continue
                
            if grain_bin is not None and sp_phase in ["surface", "mantle"]:
                if sp_bin != str(grain_bin):
                    continue
            
            coef = count_target_atom(clean_formula, target_atom)
            if coef > 0:
                valid_list.append(sp)
                coeffs.append(coef)
                
        return valid_list, np.array(coeffs)

    def clean_molec(mol_name):
        """Cleans and isolates LaTeX subscripts/superscripts for the chemical formula without grain environments."""
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        return f"${f}$"

    global_species_contributions = {}
    AU_to_cm = 1.496e13  

    # --- 1. RECONSTRUCT RADIAL EDGES FOR VOLUME ---
    radii_map = {}
    for original_key in main_output_dict.keys():
        digits = re.findall(r'\d+', str(original_key))
        if digits:
            try:
                rad_int = int(digits[0])
                radii_map[rad_int] = original_key
            except ValueError:
                continue
                
    radii = sorted(list(radii_map.keys()))
    
    if len(radii) == 0:
        if verbose: print("Main output dictionary is empty or contains no valid numerical keys.")
        return
        
    if len(radii) > 1:
        r_midshifts = 0.5 * np.diff(radii)
        r_edges = [radii[0] - r_midshifts[0]] + [radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [radii[-1] + r_midshifts[-1]]
    else:
        r_edges = [radii[0] * 0.9, radii[0] * 1.1]

    # --- 2. DATA ACCUMULATION WITH PHYSICAL VOLUMES ---
    for i, r_value in enumerate(radii):
        folder_name = f"{r_value}AU"
        file_path = os.path.join(chempath, folder_name, "1D_static.dat")
        
        if os.path.exists(file_path):
            try:
                z_points = np.loadtxt(file_path, comments='!', usecols=0)
                orig_key = radii_map[r_value]
                sub_dict = main_output_dict[orig_key]
                
                abundance_array = sub_dict['abundances']
                nH_profile = sub_dict["H_number_density"][itime,:]
                
                if len(z_points) > 1:
                    z_midshifts = 0.5 * np.diff(z_points)
                    z_edges = [z_points[0] - z_midshifts[0]] + [z_points[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_points[-1] + z_midshifts[-1])]
                    dz = np.abs(np.diff(z_edges))
                else:
                    dz = np.array([z_points[0] if z_points[0] > 0 else 1.0])
                
                r_left = r_edges[i] * AU_to_cm
                r_right = r_edges[i+1] * AU_to_cm
                dR = r_right - r_left
                R_center = float(r_value) * AU_to_cm
                
                cell_volumes = 2 * np.pi * R_center * dR * (dz * AU_to_cm) 
                
                raw_species = list(abundance_array.coords['species'].values)
                local_species_list, target_coeffs = filter_and_cache_species(raw_species)
                
                if not local_species_list:
                    continue 
                
                y_abundances = abundance_array.isel(time=itime).sel(species=local_species_list).values
                
                nH_2d = nH_profile[np.newaxis, :]     
                volumes_2d = cell_volumes[np.newaxis, :] 
                
                physical_atoms = y_abundances * nH_2d * volumes_2d
                total_column_atoms = np.sum(physical_atoms, axis=1)
                contributions_per_species = total_column_atoms * target_coeffs
                
                for idx, species in enumerate(local_species_list):
                    if contributions_per_species[idx] > 0:
                        _, _, clean_formula = parse_species(species)
                        global_species_contributions[clean_formula] = global_species_contributions.get(clean_formula, 0.0) + contributions_per_species[idx]
                        
            except Exception as e:
                if verbose: print(f"Error processing data for R={r_value}: {e}")
        elif verbose: 
            print(f"File not found: {file_path}")

    # --- 3. STATISTICAL PROCESSING & PLOTTING ---
    total_atoms_sum = sum(global_species_contributions.values())
    if total_atoms_sum == 0:
        bin_str = f" (Bin: {grain_bin})" if grain_bin else ""
        print(f"No {target_atom} atoms detected within the chosen criteria [Phase: {phase}{bin_str}].")
        return

    species_percentages = {sp: (val / total_atoms_sum) * 100 for sp, val in global_species_contributions.items()}
    sorted_species = sorted(species_percentages.items(), key=lambda x: x[1], reverse=True)
    
    top = sorted_species[:spnumber]
    others_sum = 100.0 - sum(val for sp, val in top)

    labels = [clean_molec(item[0]) for item in top]
    percentages = [item[1] for item in top]

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(labels, percentages, color=color, edgecolor='grey', alpha=0.85)
    
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.4f}%',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 4),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Reconstruct localized grain properties if an ice-phase bin value exists
    is_ice_phase = phase in ["surface", "mantle", "grain"]
    if grain_bin and is_ice_phase:
        first_r = radii[0]
        size_um = get_grain_size_in_um(chempath / f"{int(first_r)}AU" / "1D_grain_sizes.in", grain_bin)
        bin_title = f" (Grain Size = {size_um:.1f} µm)" if size_um is not None else f" (Bin {grain_bin})"
    else:
        bin_title = " (All Grains)" if is_ice_phase else ""
        
    phase_title = f"Phase: {phase.upper()}{bin_title}"

    ax.set_ylabel(f'Global {phase_title} Budget Contribution (%)', fontsize=11)
    ax.set_xlabel(f'Chemical Species (Top {spnumber})', fontsize=12)
    
    try:
        first_r = radii[0]
        time_seconds = main_output_dict[radii_map[first_r]]['abundances'].coords['time'].values[itime]
        ax.set_title(f'Top {spnumber} Reservoirs for Element [{target_atom}] — {phase_title}\n$t = {time_seconds/3.156e7:.0f}$ years', fontsize=13, pad=15)
    except:
        ax.set_title(f'Top {spnumber} Reservoirs for Element [{target_atom}] — {phase_title}', fontsize=13, pad=15)
        
    ax.set_ylim(0, max(percentages) * 1.15 if percentages else 100)
    ax.grid(axis='y', linestyle='--', alpha=0.4)

    if others_sum > 0.0001:
        plt.figtext(0.15, -0.01, f"* Omitted minor species account for {others_sum:.4f}% of the filtered physical {target_atom} budget.", 
                    fontsize=10, style='italic')

    plt.tight_layout()
    plt.show()

def plot_top_species_per_radius(chempath,
                                main_output_dict,
                                target_atom="C",
                                itime=-1,
                                verbose=True,
                                spnumber=5,
                                phase="gas",
                                grain_bin=None,
                                cmap_name="tab10",
                                rmin=None,
                                rmax=None):
    """
    Computes and plots the top N contributing chemical species for a target element
    as a function of disk radius within an optional [rmin, rmax] radial range. 

    Displays a horizontal bar chart grouped by radius on the Y-axis, with local budget 
    percentages on the X-axis. Each unique chemical species is mapped to a distinct 
    color from the chosen colormap, formatting molecule labels inside LaTeX blocks.

    Parameters:
    -----------
    chempath : str
        Path to the directory containing spatial grid subfolders (e.g., '10AU/', '100AU/').
    main_output_dict : dict
        Nested dictionary mapping radial keys to data sub-structures containing 'abundances' 
        and 'H_number_density' spatial profiles.
    target_atom : str, default "C"
        The chemical element symbol whose structural reservoir distribution is evaluated.
    itime : int, default -1
        Time index to slice from the multi-epoch data arrays.
    verbose : bool, default True
        If True, outputs warning and missing file notifications to the console.
    spnumber : int, default 5
        Number of top contributing chemical species to display for each radius.
    phase : str, default "gas"
        The chemical phase environment to isolate ("gas", "surface", "mantle", "grain", or "all").
    grain_bin : int or str, optional
        Specific grain size category to filter when analyzing surface or mantle ice matrices.
    cmap_name : str, default "tab10"
        Name of the Matplotlib colormap used to dynamically assign unique colors to species.
    rmin : float or int, optional
        Minimum inclusive radius boundary (in AU) to process. If None, defaults to the disk innermost grid point.
    rmax : float or int, optional
        Maximum inclusive radius boundary (in AU) to process. If None, defaults to the disk outermost grid point.
    """
    
    # Allowed elements check
    allowed_elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    if target_atom not in allowed_elements:
        raise ValueError(f"Target atom '{target_atom}' is not recognized in the chemical network.")
    
    # Phase validation
    valid_phases = ["gas", "surface", "mantle", "grain", "all"]
    if phase not in valid_phases:
        raise ValueError(f"Phase '{phase}' unrecognized. Choose among {valid_phases}")
    if grain_bin is not None and phase in ["gas","all"] : raise ValueError("grain_bin and gas phase can not be defined simultaneously")

    chempath = Path(chempath)

    # --- INTERNAL HELPERS ---
    def get_grain_size_in_um(file_path, bin_index):
        """Parses 1D_grain_sizes.in to retrieve the grain bin radius mapped to micrometers."""
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('!'):
                        continue
                    if '!' in line:
                        line = line.split('!')[0].strip()
                    values = [float(val) for val in line.split()]
                    if not values:
                        continue
                    num_grains = len(values) // 4
                    radii_cm = values[:num_grains]
                    index = int(bin_index) - 1
                    if 0 <= index < num_grains:
                        return radii_cm[index] * 10000.0
            return None
        except FileNotFoundError:
            return None

    def parse_species(species_name):
        if species_name == 'e-':
            return "gas", None, "e-"
            
        grain_match = re.match(r'^([JK])(\d+)(.+)', species_name)
        
        if grain_match:
            p_code, g_bin, raw_formula = grain_match.groups()
            sp_phase = "surface" if p_code == 'J' else "mantle"
            sp_bin = g_bin
        else:
            sp_phase = "gas"
            sp_bin = None
            raw_formula = species_name
        
        clean_formula = raw_formula.replace('c-', '').replace('l-', '')
        return sp_phase, sp_bin, clean_formula

    def count_target_atom(clean_formula, target):
        if clean_formula == 'e-': return 0
        
        calc_formula = clean_formula
        if calc_formula.endswith('+') or calc_formula.endswith('-'):
            calc_formula = calc_formula[:-1]
            
        calc_formula = calc_formula.replace('-', '')
            
        pattern = re.compile(r'([A-Z][a-z]?)(\d*)') 
        composition = {}
        for atom, n in pattern.findall(calc_formula):
            try:
                count = int(n) if n else 1
            except ValueError:
                count = 1
            composition[atom] = composition.get(atom, 0) + count
        return composition.get(target, 0)

    def filter_and_cache_species(species_list):
        valid_list = []
        coeffs = []
        for sp in species_list:
            if "GRAIN" in sp:
                continue
            
            sp_phase, sp_bin, clean_formula = parse_species(sp)
            
            # Phase filtering
            if phase == "all":
                pass
            elif phase == "grain":
                if sp_phase not in ["surface", "mantle"]:
                    continue
            elif sp_phase != phase:
                continue
                
            # Size bin filtering
            if grain_bin is not None and sp_phase in ["surface", "mantle"]:
                if sp_bin != str(grain_bin):
                    continue
                    
            coef = count_target_atom(clean_formula, target_atom)
            if coef > 0:
                valid_list.append(sp)
                coeffs.append(coef)
        return valid_list, np.array(coeffs)

    def clean_molec(mol_name):
        """Cleans and isolates LaTeX subscripts/superscripts for the chemical formula without grain environments."""
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        return f"${f}$"

    AU_to_cm = 1.496e13  

    # --- 1. RECONSTRUCT RADIAL EDGES FOR VOLUME ---
    radii_map = {}
    for original_key in main_output_dict.keys():
        digits = re.findall(r'\d+', str(original_key))
        if digits:
            try:
                rad_int = int(digits[0])
                radii_map[rad_int] = original_key
            except ValueError:
                continue
                
    extracted_radii = sorted(list(radii_map.keys()))
    
    # --- RADIAL LIMITS FILTERING ---
    radii = []
    for r in extracted_radii:
        if rmin is not None and r < rmin:
            continue
        if rmax is not None and r > rmax:
            continue
        radii.append(r)
    
    if len(radii) == 0:
        if verbose: print(f"No grid radii found within the specified boundaries [rmin={rmin}, rmax={rmax}].")
        return
        
    if len(radii) > 1:
        r_midshifts = 0.5 * np.diff(radii)
        r_edges = [radii[0] - r_midshifts[0]] + [radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [radii[-1] + r_midshifts[-1]]
    else:
        r_edges = [radii[0] * 0.9, radii[0] * 1.1]

    # Structure to hold plot data per radius
    radial_plot_data = {}
    all_encountered_species = set()

    # --- 2. DATA ACCUMULATION PER CELL SHELL ---
    for i, r_value in enumerate(radii):
        folder_name = f"{r_value}AU"
        file_path = os.path.join(chempath, folder_name, "1D_static.dat")
        
        if os.path.exists(file_path):
            try:
                z_points = np.loadtxt(file_path, comments='!', usecols=0)
                orig_key = radii_map[r_value]
                sub_dict = main_output_dict[orig_key]
                
                abundance_array = sub_dict['abundances']
                nH_profile = sub_dict["H_number_density"][itime,:]
                
                if len(z_points) > 1:
                    z_midshifts = 0.5 * np.diff(z_points)
                    z_edges = [z_points[0] - z_midshifts[0]] + [z_points[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_points[-1] + z_midshifts[-1])]
                    dz = np.abs(np.diff(z_edges))
                else:
                    dz = np.array([z_points[0] if z_points[0] > 0 else 1.0])
                
                r_left = r_edges[i] * AU_to_cm
                r_right = r_edges[i+1] * AU_to_cm
                dR = r_right - r_left
                R_center = float(r_value) * AU_to_cm
                
                cell_volumes = 2 * np.pi * R_center * dR * (dz * AU_to_cm) 
                
                raw_species = list(abundance_array.coords['species'].values)
                local_species_list, target_coeffs = filter_and_cache_species(raw_species)
                
                if not local_species_list:
                    continue 
                
                y_abundances = abundance_array.isel(time=itime).sel(species=local_species_list).values
                
                nH_2d = nH_profile[np.newaxis, :]     
                volumes_2d = cell_volumes[np.newaxis, :] 
                
                physical_atoms = y_abundances * nH_2d * volumes_2d
                total_column_atoms = np.sum(physical_atoms, axis=1)
                contributions_per_species = total_column_atoms * target_coeffs
                
                local_radius_contributions = {}
                for idx, species in enumerate(local_species_list):
                    if contributions_per_species[idx] > 0:
                        _, _, clean_formula = parse_species(species)
                        local_radius_contributions[clean_formula] = local_radius_contributions.get(clean_formula, 0.0) + contributions_per_species[idx]
                
                radius_total_sum = sum(local_radius_contributions.values())
                if radius_total_sum > 0:
                    species_percentages = {sp: (val / radius_total_sum) * 100 for sp, val in local_radius_contributions.items()}
                    sorted_species = sorted(species_percentages.items(), key=lambda x: x[1], reverse=True)
                    top_entries = sorted_species[:spnumber]
                    radial_plot_data[r_value] = top_entries
                    
                    for species_name, _ in top_entries:
                        all_encountered_species.add(species_name)
                        
            except Exception as e:
                if verbose: print(f"Error processing data for R={r_value}: {e}")
        elif verbose: 
            print(f"File not found: {file_path}")

    if not radial_plot_data:
        print(f"No data available to plot for target atom {target_atom} within the specified radius range.")
        return

    # --- 3. DYNAMIC COLOR ASSIGNMENT VIA COLORMAP ---
    unique_species_list = sorted(list(all_encountered_species))
    num_unique_species = len(unique_species_list)
    
    try:
        cmap = plt.get_cmap(cmap_name)
    except ValueError:
        if verbose: print(f"Colormap '{cmap_name}' not found. Falling back to default 'tab10'.")
        cmap = plt.get_cmap("tab10")
        
    if hasattr(cmap, 'colors') and len(cmap.colors) >= num_unique_species:
        species_color_mapping = {sp: cmap(idx) for idx, sp in enumerate(unique_species_list)}
    else:
        color_indices = np.linspace(0, 1, num_unique_species) if num_unique_species > 1 else [0.0]
        species_color_mapping = {sp: cmap(color_indices[idx]) for idx, sp in enumerate(unique_species_list)}

    # --- 4. DYNAMIC GRAPHICAL CONSTRUCTION ---
    y_labels = []
    x_percentages = []
    bar_colors = []
    
    reversed_radii = sorted(list(radial_plot_data.keys()), reverse=True)
    group_edges = []
    current_index = 0

    for r_val in reversed_radii:
        top_entries = radial_plot_data[r_val]
        for species_name, pct_val in reversed(top_entries):
            # Parse species keys inside math/latex formatting text blocks on the horizontal axis labels
            latex_formula = clean_molec(species_name)
            y_labels.append(f"{r_val} AU — {latex_formula}")
            x_percentages.append(pct_val)
            bar_colors.append(species_color_mapping[species_name])
            current_index += 1
        group_edges.append(current_index)

    fig_height = max(4, len(y_labels) * 0.45)
    fig, ax = plt.subplots(figsize=(11, fig_height))
    
    bars = ax.barh(y_labels, x_percentages, color=bar_colors, edgecolor='grey', alpha=0.85, height=0.7)
    
    for bar in bars:
        width = bar.get_width()
        ax.annotate(f' {width:.2f}%',
                    xy=(width, bar.get_y() + bar.get_height() / 2),
                    xytext=(4, 0),
                    textcoords="offset points",
                    ha='left', va='center', fontsize=9, fontweight='bold')

    for edge in group_edges[:-1]:
        ax.axhline(y=edge - 0.5, color='black', linestyle='-', alpha=0.3, linewidth=1.2)

    # Reconstruct localized grain parameters safely if an ice-phase bin is supplied
    is_ice_phase = phase in ["surface", "mantle", "grain"]
    if grain_bin and is_ice_phase:
        first_r = radii[0]
        size_um = get_grain_size_in_um(chempath / f"{int(first_r)}AU" / "1D_grain_sizes.in", grain_bin)
        bin_title = f" (Grain Size = {size_um:.1f} µm)" if size_um is not None else f" (Bin {grain_bin})"
    else:
        bin_title = " (All Grains)" if is_ice_phase else ""
        
    phase_title = f"Phase: {phase.upper()}{bin_title}"

    ax.set_xlabel('Local Radial Budget Contribution (%)', fontsize=12)
    ax.set_ylabel('Disk Radius & Associated Chemical Carrier', fontsize=12)
    
    try:
        first_r = radii[0]
        time_seconds = main_output_dict[radii_map[first_r]]['abundances'].coords['time'].values[itime]
        ax.set_title(f'Top {spnumber} Radial Reservoirs for Element [{target_atom}] — {phase_title}\n$t = {time_seconds/3.156e7:.0f}$ years', fontsize=13, pad=15)
    except:
        ax.set_title(f'Top {spnumber} Radial Reservoirs for Element [{target_atom}] — {phase_title}', fontsize=13, pad=15)
        
    # Bounds configuration
    max_pct = max(x_percentages) if x_percentages else 100.0
    ax.set_xlim(0, min(max_pct * 1.05, 100) + 8.0)  
    ax.set_ylim(-0.6, len(y_labels) - 0.4)
    
    ax.grid(axis='x', linestyle='--', alpha=0.4)

    plt.tight_layout()
    plt.show()

# def plot_disk_atomic_composition(chempath,
#                                  main_output_dict,
#                                  itime=-1,
#                                  verbose=True,
#                                  grain_bin=None,
#                                  cmap_name="Set1"):
#     """
#     Computes and plots the total volume-integrated atomic composition of the protoplanetary 
#     disk relative to Hydrogen within each specific phase (X_phase / H_phase). 
#     Elements are sorted by atomic mass on the X-axis (excluding 'X').
#     The Y-axis displays the local abundance ratio in log scale across 5 chemical phases.

#     Parameters:
#     -----------
#     chempath : str
#         Path to the directory containing spatial grid subfolders (e.g., '10AU/', '100AU/').
#     main_output_dict : dict
#         Nested dictionary mapping radial keys to data sub-structures containing 'abundances' 
#         and 'H_number_density' spatial profiles.
#     itime : int, default -1
#         Time index to slice from the multi-epoch data arrays.
#     verbose : bool, default True
#         If True, outputs warning and missing file notifications to the console.
#     grain_bin : int or str, optional
#         Specific grain size category to filter when analyzing ice layers. If specified, 
#         gas phase contributions are completely ignored for safety.
#     cmap_name : str, default "Set1"
#         Name of the Matplotlib colormap used to dynamically assign unique colors to the 5 phases.
#     """

#     # Dictionary of standard atomic masses used to strictly sort the X-axis
#     atomic_masses = {
#         'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81, 'C': 12.011, 
#         'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305, 
#         'Al': 26.982, 'Si': 28.085, 'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 
#         'K': 39.098, 'Ca': 40.078, 'Fe': 55.845, 'Ni': 58.693
#     }

#     # --- INTERNAL HELPERS ---
#     def parse_species(species_name):
#         if species_name == 'e-':
#             return "gas", None, "e-"
            
#         grain_match = re.match(r'^([JK])(\d+)(.+)', species_name)
        
#         if grain_match:
#             p_code, g_bin, raw_formula = grain_match.groups()
#             sp_phase = "surface" if p_code == 'J' else "mantle"
#             sp_bin = g_bin
#         else:
#             sp_phase = "gas"
#             sp_bin = None
#             raw_formula = species_name
        
#         clean_formula = raw_formula.replace('c-', '').replace('l-', '')
#         return sp_phase, sp_bin, clean_formula

#     def get_chemical_composition(clean_formula):
#         if clean_formula == 'e-': 
#             return {}
        
#         calc_formula = clean_formula
#         if calc_formula.endswith('+') or calc_formula.endswith('-'):
#             calc_formula = calc_formula[:-1]
            
#         calc_formula = calc_formula.replace('-', '')
            
#         pattern = re.compile(r'([A-Z][a-z]?)(\d*)') 
#         composition = {}
#         for atom, n in pattern.findall(calc_formula):
#             if atom == 'X':
#                 continue
#             try:
#                 count = int(n) if n else 1
#             except ValueError:
#                 count = 1
#             composition[atom] = composition.get(atom, 0) + count
#         return composition

#     AU_to_cm = 1.496e13  
#     phases_pool = ["all", "grain", "mantle", "surface", "gas"]

#     # --- 1. RECONSTRUCT RADIAL GRID ---
#     radii_map = {}
#     for original_key in main_output_dict.keys():
#         digits = re.findall(r'\d+', str(original_key))
#         if digits:
#             try:
#                 rad_int = int(digits[0])
#                 radii_map[rad_int] = original_key
#             except ValueError:
#                 continue
                
#     radii = sorted(list(radii_map.keys()))
    
#     if len(radii) == 0:
#         if verbose: print("Main output dictionary is empty or contains no valid numerical keys.")
#         return
        
#     if len(radii) > 1:
#         r_midshifts = 0.5 * np.diff(radii)
#         r_edges = [radii[0] - r_midshifts[0]] + [radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [radii[-1] + r_midshifts[-1]]
#     else:
#         r_edges = [radii[0] * 0.9, radii[0] * 1.1]

#     # --- 2. PRE-PARSING NETWORK SPECIES TO BUILD STATIC MAPS ---
#     first_r_key = radii_map[radii[0]]
#     all_network_species = list(main_output_dict[first_r_key]['abundances'].coords['species'].values)
    
#     species_metadata = []
#     unique_elements = set()
    
#     for idx, sp in enumerate(all_network_species):
#         if "GRAIN" in sp:
#             species_metadata.append((None, None, {}))
#             continue
#         sp_phase, sp_bin, clean_formula = parse_species(sp)
#         elem_map = get_chemical_composition(clean_formula)
#         species_metadata.append((sp_phase, sp_bin, elem_map))
#         for elem in elem_map.keys():
#             unique_elements.add(elem)
            
#     # Sort elements strictly by atomic mass weights
#     unique_elements = sorted(list(unique_elements), key=lambda e: atomic_masses.get(e, 999.0))
#     num_elements = len(unique_elements)

#     # Dictionary to store accumulated disk-wide total absolute atom sums
#     disk_total_budgets = {p: {elem: 0.0 for elem in unique_elements} for p in phases_pool}

#     # --- 3. DATA ACCUMULATION PER CELL SHELL (VECTORIZED) ---
#     for i, r_value in enumerate(radii):
#         folder_name = f"{r_value}AU"
#         file_path = os.path.join(chempath, folder_name, "1D_static.dat")
        
#         if os.path.exists(file_path):
#             try:
#                 z_points = np.loadtxt(file_path, comments='!', usecols=0)
#                 orig_key = radii_map[r_value]
#                 sub_dict = main_output_dict[orig_key]
                
#                 abundance_array = sub_dict['abundances']
#                 nH_profile = sub_dict["H_number_density"][itime,:]
                
#                 if len(z_points) > 1:
#                     z_midshifts = 0.5 * np.diff(z_points)
#                     z_edges = [z_points[0] - z_midshifts[0]] + [z_points[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_points[-1] + z_midshifts[-1])]
#                     dz = np.abs(np.diff(z_edges))
#                 else:
#                     dz = np.array([z_points[0] if z_points[0] > 0 else 1.0])
                
#                 r_left = r_edges[i] * AU_to_cm
#                 r_right = r_edges[i+1] * AU_to_cm
#                 dR = r_right - r_left
#                 R_center = float(r_value) * AU_to_cm
                
#                 cell_volumes = 2 * np.pi * R_center * dR * (dz * AU_to_cm) 
                
#                 y_abundances_2d = abundance_array.isel(time=itime).values
#                 integrated_molecules_per_species = np.sum(y_abundances_2d * nH_profile * cell_volumes, axis=1)
                
#                 coef_matrices = {p: np.zeros((num_elements, len(all_network_species))) for p in phases_pool}
                
#                 for sp_idx, (sp_phase, sp_bin, elem_map) in enumerate(species_metadata):
#                     if sp_phase is None: 
#                         continue
                        
#                     if grain_bin is not None:
#                         if sp_phase == "gas":
#                             continue
#                         if sp_bin != str(grain_bin):
#                             continue
                    
#                     for elem_idx, elem in enumerate(unique_elements):
#                         coef = elem_map.get(elem, 0)
#                         if coef > 0:
#                             coef_matrices[sp_phase][elem_idx, sp_idx] = coef
#                             if sp_phase in ["surface", "mantle"]:
#                                 coef_matrices["grain"][elem_idx, sp_idx] = coef
#                             coef_matrices["all"][elem_idx, sp_idx] = coef

#                 for p in phases_pool:
#                     integrated_atoms_per_element = np.dot(coef_matrices[p], integrated_molecules_per_species)
#                     for elem_idx, elem in enumerate(unique_elements):
#                         disk_total_budgets[p][elem] += integrated_atoms_per_element[elem_idx]

#             except Exception as e:
#                 if verbose: print(f"Error processing cell matrix for R={r_value}: {e}")
#         elif verbose: 
#             print(f"File not found: {file_path}")

#     # --- 4. COLOR ASSIGNMENT VIA COLORMAP ---
#     try:
#         cmap = plt.get_cmap(cmap_name)
#     except ValueError:
#         if verbose: print(f"Colormap '{cmap_name}' not found. Falling back to default 'Set1'.")
#         cmap = plt.get_cmap("Set1")
        
#     phase_colors = {phase_key: cmap(idx) for idx, phase_key in enumerate(phases_pool)}

# # --- 5. GRAPH CONSTRUCTION WITH LOCAL PHASE HYDROGEN NORMALIZATION ---
#     fig, ax = plt.subplots(figsize=(11, 6))
#     x_indexes = np.arange(num_elements)

#     # Ajustement dynamique des phases à tracer selon la présence d'un grain_bin
#     phases_to_plot = ["grain", "mantle", "surface"] if grain_bin is not None else phases_pool

#     for phase_key in phases_to_plot:
#         abundance_ratios_to_local_H = []
        
#         # Pull out the local Hydrogen reference count for THIS specific phase pool
#         local_hydrogen_reference = disk_total_budgets[phase_key].get("H", 0.0)
        
#         for elem in unique_elements:
#             absolute_phase_val = disk_total_budgets[phase_key][elem]
            
#             # Divide the element phase budget by the local phase Hydrogen count
#             if local_hydrogen_reference > 0:
#                 abundance_ratios_to_local_H.append(absolute_phase_val / local_hydrogen_reference)
#             else:
#                 abundance_ratios_to_local_H.append(0.0)
                
#         abundance_ratios_to_local_H = np.array(abundance_ratios_to_local_H)
        
#         # Mask 0 values to protect logarithmic scale bounds rendering
#         valid_mask = abundance_ratios_to_local_H > 0
        
#         if np.any(valid_mask):
#             ax.plot(x_indexes[valid_mask], abundance_ratios_to_local_H[valid_mask], 
#                     label=phase_key.upper(), 
#                     color=phase_colors[phase_key], 
#                     linewidth=2.5, 
#                     marker='o', 
#                     markersize=6, 
#                     alpha=0.9)

#     # Format graph grids and axis layouts
#     ax.set_yscale('log')
#     ax.set_ylabel('Abundance Ratio (X_phase / H_phase)', fontsize=12)
#     ax.set_xlabel('Chemical Elements (Sorted by Atomic Mass Weights)', fontsize=12)
    
#     # Tie the sorted elements text tracking sequence to X ticks
#     ax.set_xticks(x_indexes)
#     ax.set_xticklabels(unique_elements, fontsize=11, fontweight='bold')
    
#     bin_suffix = f" (Grain Size Bin {grain_bin})" if grain_bin else ""
#     ax.set_title(f"Protoplanetary Disk phase-dependant Elemental Abundances (X/H per phase){bin_suffix}", fontsize=13, fontweight='bold', pad=12)
    
#     ax.grid(True, which="both", linestyle="--", alpha=0.4)
#     ax.legend(title="Chemical Phases", title_fontsize='11', loc='upper right', frameon=True)

#     try:
#         first_r = radii[0]
#         time_seconds = main_output_dict[radii_map[first_r]]['abundances'].coords['time'].values[itime]
#         plt.figtext(0.15, -0.01, f"Data integrated over disk total volume at epoch: $t = {time_seconds/3.156e7:.2e}$ years", 
#                     fontsize=10, style='italic')
#     except:
#         pass

#     plt.tight_layout()
#     plt.show()



import os
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def plot_species_evolution_with_grain_size_comparison(PIPE,
                                                     MODEL_NAMES,
                                                     target_radius,
                                                     itime=-1,
                                                     verbose=True,
                                                     spnumber=5,
                                                     cmap_name="tab10"):
    r"""
    Computes and plots the evolution of the globally dominant ice chemical species 
    (Surface + Mantle combined) as a function of grain size across multiple models.

    Dynamically manages structural subplots depending on the number of models in PIPE:
    - 1 Model: Displays standard grid layout for selected radii (max 3 columns).
    - 2 Models: Displays a 3-column layout per radius row (Model 1, Model 2, and Residuals).
    - >2 Models: Displays a multi-column grid where columns represent models and rows represent radii.

    Parameters:
    -----------
    PIPE : list
        Collection of model pipe objects to analyze. Each object must feature 
        `.chemistry` and `.chempath` attributes.
    MODEL_NAMES : list of str
        Display names assigned to the models for legends. If duplicate names 
        are found, a generic fallback nomenclature is generated.
    target_radius : int, float or list
        The specific radius or list of radii (in AU) at which the grain size analysis is performed.
    itime : int, default -1
        Time index to slice from the multi-epoch data arrays.
    verbose : bool, default True
        If True, outputs warning and missing file notifications to the console.
    spnumber : int, default 5
        Number of top global chemical species to display across all bins.
    cmap_name : str, default "tab10"
        Name of the Matplotlib colormap used to dynamically assign unique colors to species.
    res_colormap : str, default "seismic"
        Matplotlib colormap string used to style background shifts or markers in residual modes.
    """
    # Standardize input radius to a list of integers
    if isinstance(target_radius, str):
        target_radius = [target_radius]
    radii_list = [target_radius] if not isinstance(target_radius, list) else target_radius
    radii_list = [int(r) for r in radii_list]

    # Handle duplicate model names via a generic nomenclature fallback
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)):
        MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    # --- INTERNAL HELPERS ---
    def get_grain_size_in_um(file_path, bin_index):
        """Parses 1D_grain_sizes.in to extract physical grain radii converted to micrometers."""
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('!'): continue
                    if '!' in line: line = line.split('!')[0].strip()
                    values = [float(val) for val in line.split()]
                    if not values: continue
                    num_grains = len(values) // 4
                    radii_cm = values[:num_grains]
                    index = int(bin_index) - 1
                    if 0 <= index < num_grains: return radii_cm[index] * 10000.0
            return None
        except FileNotFoundError:
            return None

    def parse_species(species_name):
        """Decodes Nautilus chemical naming strings into phase type, bin index, and formula."""
        if species_name == 'e-': return "gas", None, "e-"
        grain_match = re.match(r'^([JK])(\d+)(.+)', species_name)
        if grain_match:
            p_code, g_bin, raw_formula = grain_match.groups()
            sp_phase = "surface" if p_code == 'J' else "mantle"
            sp_bin = g_bin 
        else:
            sp_phase = "gas"
            sp_bin = None
            raw_formula = species_name
        clean_formula = raw_formula.replace('c-', '').replace('l-', '')
        return sp_phase, sp_bin, clean_formula

    def clean_molec(mol_name):
        """Converts raw chemical strings into standard formatted LaTeX chemical expressions."""
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        return f"${f}$"

    AU_to_cm = 1.496e13  # Astronomical Unit conversion factor

    # --- MULTI-MODEL DATABASE HARVESTING ---
    model_data = {}  
    all_encountered_species = set()
    time_years_string = None

    # Loop over pipelines
    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        main_output_dict = p.chemistry
        chempath = Path(p.chempath)
        
        # Build radial index mapping from keys containing numerical values
        radii_map = {}
        for original_key in main_output_dict.keys():
            digits = re.findall(r'\d+', str(original_key))
            if digits:
                try: radii_map[int(digits[0])] = original_key
                except ValueError: continue
                    
        sorted_all_radii = sorted(list(radii_map.keys()))
        model_data[p_name] = {}

        # Loop over target radii
        for r_value in radii_list:
            if r_value not in radii_map:
                if verbose: print(f"[{p_name}] Radius {r_value} AU missing. Skipping.")
                continue
                
            orig_key = radii_map[r_value]
            sub_dict = main_output_dict[orig_key]
            abundance_array = sub_dict['abundances']
            nH_profile = sub_dict["H_number_density"][itime,:]
            raw_species_list = list(abundance_array.coords['species'].values)

            # Reconstruct 1D radial mid-shifts cell boundaries
            r_idx = sorted_all_radii.index(r_value)
            if len(sorted_all_radii) > 1:
                r_midshifts = 0.5 * np.diff(sorted_all_radii)
                r_edges = [sorted_all_radii[0] - r_midshifts[0]] + [sorted_all_radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [sorted_all_radii[-1] - r_midshifts[-1]]
            else:
                r_edges = [r_value * 0.9, r_value * 1.1]

            # Reconstruct vertical grid cell boundaries from static file configurations
            file_path = os.path.join(chempath, f"{r_value}AU", "1D_static.dat")
            if not os.path.exists(file_path):
                if verbose: print(f"[{p_name}] Static file missing at {r_value}AU. Skipping.")
                continue
                
            z_points = np.loadtxt(file_path, comments='!', usecols=0)
            if len(z_points) > 1:
                z_midshifts = 0.5 * np.diff(z_points)
                z_edges = [z_points[0] - z_midshifts[0]] + [z_points[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_points[-1] + z_midshifts[-1])]
                dz = np.abs(np.diff(z_edges))
            else:
                dz = np.array([z_points[0] if z_points[0] > 0 else 1.0])

            # Calculate 3D cylindrical integration volumes per cell column
            r_left = r_edges[r_idx] * AU_to_cm
            r_right = r_edges[r_idx+1] * AU_to_cm
            cell_volumes = 2 * np.pi * (float(r_value) * AU_to_cm) * (r_right - r_left) * (dz * AU_to_cm)
            
            nH_2d = nH_profile[np.newaxis, :]
            volumes_2d = cell_volumes[np.newaxis, :]

            available_bins_set = set()
            species_metadata = {} 

            # Map valid grain phase targets
            for sp in raw_species_list:
                if "GRAIN" in sp: continue
                sp_phase, sp_bin, clean_formula = parse_species(sp)
                if sp_phase in ["surface", "mantle"] and sp_bin is not None:
                    available_bins_set.add(sp_bin)
                    species_metadata[sp] = (sp_phase, sp_bin, clean_formula)

            try: sorted_bins = sorted(list(available_bins_set), key=lambda x: int(x))
            except ValueError: sorted_bins = sorted(list(available_bins_set))

            if not sorted_bins: continue

            bin_raw_data = {b: {} for b in sorted_bins}
            local_species_scores = {}

            # Integrate absolute particle counts per grain size bin
            for sp in raw_species_list:
                if sp not in species_metadata: continue
                _, sp_bin, clean_formula = species_metadata[sp]

                y_values = abundance_array.isel(time=itime).sel(species=sp).values
                absolute_particles = np.sum(y_values * nH_2d * volumes_2d)
                
                if absolute_particles > 0:
                    bin_raw_data[sp_bin][clean_formula] = bin_raw_data[sp_bin].get(clean_formula, 0.0) + absolute_particles
                    local_species_scores[clean_formula] = local_species_scores.get(clean_formula, 0.0) + absolute_particles

            if not local_species_scores: continue

            # Track top dominating unique carriers globally across bins
            sorted_local = sorted(local_species_scores.items(), key=lambda x: x[1], reverse=True)
            top_local_species = [item[0] for item in sorted_local[:spnumber]]
            
            for sp in top_local_species:
                all_encountered_species.add(sp)

            if time_years_string is None:
                try:
                    t_sec = abundance_array.coords['time'].values[itime]
                    time_years_string = f"{t_sec / 3.156e7:.0f}"
                except: pass

            model_data[p_name][r_value] = {
                'bin_raw_data': bin_raw_data,
                'sorted_bins': sorted_bins,
                'top_species': top_local_species,
                'chempath': chempath
            }

    # --- GEOMETRIC SUBPLOT CONFIGURATION PANEL ---
    num_models = len(PIPE)
    model_names = list(model_data.keys())
    num_radii = len(radii_list)

    if num_models == 1:
        cols = min(3, num_radii)
        rows = (num_radii + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 5.5 * rows), squeeze=False, sharey=True)
        axes = axes.flatten()
    elif num_models == 2:
        cols = 3  # Configuration matrix: [Model 1] [Model 2] [Residuals Map]
        rows = num_radii
        fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows), squeeze=False, sharey=False)
    else:
        cols = num_models
        rows = num_radii
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5 * rows), squeeze=False, sharey=True)

    # --- SYSTEMWIDE COLORMAPPING STABILITY ---
    unique_species_list = sorted(list(all_encountered_species))
    num_unique_species = len(unique_species_list)
    
    cmap = plt.colormaps.get_cmap(cmap_name) if cmap_name in plt.colormaps else plt.colormaps["tab10"]
    if hasattr(cmap, 'colors') and len(cmap.colors) >= num_unique_species:
        species_colors = {sp: cmap(idx) for idx, sp in enumerate(unique_species_list)}
    else:
        color_indices = np.linspace(0, 1, num_unique_species) if num_unique_species > 1 else [0.0]
        species_colors = {sp: cmap(color_indices[idx]) for idx, sp in enumerate(unique_species_list)}

    # --- PLOTTING ENGINE LOOP ---
    for row_idx, r_value in enumerate(radii_list):
        
        parsed_percentages_per_model = {}
        shared_sbins = None
        shared_grain_sizes_um = []

        # 1. Base Profiles Processing Loop
        for col_idx, p_name in enumerate(model_names):
            if r_value not in model_data[p_name]: continue
            
            ax = axes[row_idx, col_idx] if num_models > 1 else axes[row_idx]
            struct = model_data[p_name][row_value] if 'row_value' in locals() else model_data[p_name][r_value]
            bin_raw = struct['bin_raw_data']
            sbins = struct['sorted_bins']
            local_top = struct['top_species']
            cpath = struct['chempath']
            
            # Map labels from bin metadata file references
            if shared_sbins is None:
                shared_sbins = sbins
                grain_sizes_file = cpath / f"{r_value}AU" / "1D_grain_sizes.in"
                for b in sbins:
                    size_um = get_grain_size_in_um(grain_sizes_file, b)
                    shared_grain_sizes_um.append(f"{size_um:.1f}" if size_um is not None else f"B{b}")

            x_positions = np.arange(len(sbins))
            plot_percentages = {sp: [] for sp in unique_species_list}
            
            # Normalize percentage allocations per bin budget
            for b in sbins:
                total_bin_budget = sum(bin_raw[b].values())
                for sp in unique_species_list:
                    if total_bin_budget > 0:
                        pct = (bin_raw[b].get(sp, 0.0) / total_bin_budget) * 100
                        plot_percentages[sp].append(pct)
                    else:
                        plot_percentages[sp].append(0.0)
                        
            parsed_percentages_per_model[p_name] = plot_percentages

            # Generate step line traces for localized dominant carriers
            for sp in local_top:
                latex_label = clean_molec(sp)
                ax.plot(x_positions, plot_percentages[sp], 
                        label=latex_label, 
                        color=species_colors[sp], 
                        linewidth=1.8, 
                        marker='o', 
                        markersize=4)

            ax.set_xlabel('Grain Radius [µm]', fontsize=11)
            ax.set_ylabel('Contribution (%)', fontsize=11)
            ax.set_title(f"{p_name} @ {r_value} AU", fontsize=11, fontweight='bold')
            ax.set_xticks(x_positions)
            ax.set_xticklabels(shared_grain_sizes_um, fontsize=8, rotation=60)
            ax.set_ylim(-2, 105)
            ax.grid(True, linestyle="--", alpha=0.4)
            ax.legend(loc='upper right', ncol=max(1, len(local_top)//2), fontsize=9)

        # 2. Render Residual Map Layout Panels (Triggered only when num_models == 2)
        if num_models == 2:
            ax_res = axes[row_idx, 2]
            p1_name, p2_name = model_names[0], model_names[1]
            
            if p1_name in parsed_percentages_per_model and p2_name in parsed_percentages_per_model:
                pct1 = parsed_percentages_per_model[p1_name]
                pct2 = parsed_percentages_per_model[p2_name]
                
                # Merge unique active tracking components from both models
                combined_top = list(set(model_data[p1_name][r_value]['top_species'] + model_data[p2_name][r_value]['top_species']))
                
                # Compute arithmetic discrepancy delta changes (Model 1 - Model 2)
                for sp in combined_top:
                    diff_vals = np.array(pct1[sp]) - np.array(pct2[sp])
                    latex_label = clean_molec(sp)
                    
                    ax_res.plot(x_positions, diff_vals, 
                                label=latex_label, 
                                color=species_colors[sp], 
                                linewidth=1.5, 
                                linestyle='--', 
                                marker='X', 
                                markersize=5)
                                
                ax_res.set_xlabel('Grain Radius [µm]', fontsize=11)
                ax_res.set_ylabel('Diff Contribution (pts)', fontsize=11)
                ax_res.set_title(f"Residuals ({p1_name} - {p2_name})", fontsize=11, fontweight='bold')
                ax_res.set_xticks(x_positions)
                ax_res.set_xticklabels(shared_grain_sizes_um, fontsize=8, rotation=60)
                
                ax_res.axhline(0, color='black', linestyle=':', alpha=0.7)
                ax_res.grid(True, linestyle="--", alpha=0.4)
                ax_res.legend(loc='best', fontsize=9)

    # Remove extra subplot layout slots in single model configurations
    if num_models == 1 and num_radii > 1:
        for idx in range(num_radii, len(axes)):
            fig.delaxes(axes[idx])

    # Assign global descriptive suptitle labels
    if time_years_string:
        if num_models == 1 and num_radii == 1:
            axes[0].set_title(f"{axes[0].get_title()} \n $t = {time_years_string}$ years")
        else:
            fig.suptitle(f'Top {spnumber} Ice Carriers vs Grain Size Distribution — $t = {time_years_string}$ years', fontsize=15, y=0.99)

    plt.tight_layout()
    plt.show()

def plot_ratio_midplane_gas_vs_grain(chempath,
                                    main_output_dict,
                                    s1="C",
                                    s2="O",
                                    itime=-1,
                                    starratio=None,
                                    verbose=True,
                                    xlim=None,
                                    ylim=None):
    """Plots the midplane (z=0) atomic abundance ratio of two elements for gas vs. grain phases.

    This function extracts multi-species chemical abundances at the disk midplane 
    across all simulated radii. It classifies species into either gas-phase or grain-phase 
    (ice surface + mantle reservoirs), computes the aggregated elemental ratio (s1/s2) 
    for each phase, and renders a comparative 1D radial line plot.

    Args:
        chempath (str): Path to the directory containing radial subfolders.
        main_output_dict (dict): Nested dictionary containing simulation outputs.
        s1 (str, optional): Atomic symbol for the numerator element. Defaults to "C".
        s2 (str, optional): Atomic symbol for the denominator element. Defaults to "O".
        itime (int, optional): Time index to slice from the abundance arrays. Defaults to -1.
        starratio (float, optional): Ratio s1/s2 of the star. Defaults to None.
        verbose (bool, optional): If True, prints status and error messages. Defaults to True.
        xlim (tuple of float, optional): Manual limits for the horizontal Radius axis.
        ylim (tuple of float, optional): Manual limits for the vertical Ratio axis.

    Raises:
        ValueError: If either s1 or s2 is not included in the allowed chemical network.
    """
    
    # Define valid chemical network elements
    elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    if s1 not in elements or s2 not in elements:
        raise ValueError("One of the specified elements does not exist in the chemical network.")

    # --- INTERNAL SPECIES PARSER AND ELEMENT COUNTER ---
    def parse_and_count(species_name, element1, element2):
        """Identifies the chemical phase and counts target atoms in a species formula."""
        # Ignore electrons and generic structural grain notations
        if species_name == 'e-' or 'GRAIN' in species_name: 
            return "ignore", 0, 0
            
        # Match surface (J) and mantle (K) ice species strings
        grain_match = re.match(r'^([JK])\d+(.+)', species_name)
        if grain_match:
            sp_phase = "grain"  
            raw_formula = grain_match.group(2)
        else:
            sp_phase = "gas"
            raw_formula = species_name
            
        # Clean up structural isomer markers, trailing charge states, and internal hyphens
        clean_formula = raw_formula.replace('c-', '').replace('l-', '')
        if clean_formula.endswith('+') or clean_formula.endswith('-'):
            clean_formula = clean_formula[:-1]
        clean_formula = clean_formula.replace('-', '')
        
        # Regex parsing to extract atomic counts
        pattern = re.compile(r'([A-Z][a-z]?)(\d*)')
        composition = {}
        for atom, n in pattern.findall(clean_formula):
            try:
                count = int(n) if n else 1
            except ValueError:
                count = 1
            composition[atom] = composition.get(atom, 0) + count
            
        return sp_phase, composition.get(element1, 0), composition.get(element2, 0)

    # --- DATA ACCUMULATION ---
    radii_list = []
    ratio_gas_list = []
    ratio_grain_list = []

    # Map string keys to numerical radius values for proper spatial sorting
    radii_map = {}
    for original_key in main_output_dict.keys():
        digits = re.findall(r'\d+', str(original_key))
        if digits:
            radii_map[int(digits[0])] = original_key
            
    sorted_radii_ints = sorted(list(radii_map.keys()))

    # Loop over sorted radial bins to retrieve midplane abundances
    for r_int in sorted_radii_ints:
        orig_key = radii_map[r_int]
        sub_dict = main_output_dict[orig_key]
        abundance_array = sub_dict['abundances']
        
        # Extract abundances at the disk midplane (deepest vertical grid cell -> index -1)
        # Expected array slice format: (species,)
        midplane_abundances = abundance_array.isel(time=itime, spatial=-1).values
        species_list = list(abundance_array.coords['species'].values)
        
        # Reset atomic accumulation pools for this radius
        total_s1_gas, total_s2_gas = 0.0, 0.0
        total_s1_grain, total_s2_grain = 0.0, 0.0
        
        # Aggregate physical atomic pools across all network species
        for idx, species in enumerate(species_list):
            phase, c1, c2 = parse_and_count(species, s1, s2)
            abundance = midplane_abundances[idx]
            
            if phase == "gas":
                total_s1_gas += abundance * c1
                total_s2_gas += abundance * c2
            elif phase == "grain":
                total_s1_grain += abundance * c1
                total_s2_grain += abundance * c2

        # Compute ratios with protective fallback division checks against empty reservoirs
        gas_ratio = total_s1_gas / total_s2_gas if total_s2_gas > 0 else 0.0
        grain_ratio = total_s1_grain / total_s2_grain if total_s2_grain > 0 else 0.0
        
        radii_list.append(float(r_int))
        ratio_gas_list.append(gas_ratio)
        ratio_grain_list.append(grain_ratio)

    if not radii_list:
        if verbose: print("No matching physical grid data could be collected.")
        return

    # Convert python collections to structured numpy arrays
    radii_arr = np.array(radii_list)
    gas_arr = np.array(ratio_gas_list)
    grain_arr = np.array(ratio_grain_list)

    # --- PLOTTING ---
    fig, ax = plt.subplots(figsize=(9, 5))
    
    # Apply a logarithmic scale to account for sharp chemical variations
    if max(np.max(gas_arr),np.max(grain_arr))/min(np.min(gas_arr),np.min(grain_arr)) > 10 : ax.set_yscale('log')
    else: ax.set_yscale('linear')
    
    # Render line tracks + markers for each discrete phase
    ax.plot(radii_arr, gas_arr, color="teal", linestyle='-', marker='o', label='Gas', linewidth=1.8)
    ax.plot(radii_arr, grain_arr, color="darkred", linestyle='--', marker='s', label='Grains (Ice)', linewidth=1.8)

    # Label definitions
    ax.set_xlabel('Radius R [AU]', fontsize=11)
    ax.set_ylabel(f'Midplane Atomic Ratio [{s1}/{s2}]', fontsize=11)
    
    # Safely parse coordinate timestamps for dynamic title labeling
    try:
        first_key = radii_map[sorted_radii_ints[0]]
        time_seconds = main_output_dict[first_key]['abundances'].coords['time'].values[itime]
        ax.set_title(f'Atomic Ratio {s1}/{s2} at Midplane ($z=0$) — $t = {time_seconds/3.156e7:.0f}$ yr', fontsize=12, pad=12)
    except:
        ax.set_title(f'Atomic Ratio {s1}/{s2} at Midplane ($z=0$)', fontsize=12, pad=12)

    if starratio is not None:
        ax.axhline(y=starratio,color='black',linestyle='--',label=f'Star {s1}/{s2} ratio')

    # Apply manual axis boundaries if supplied
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)
    
    # Auto-scale vertical bounds while avoiding log(0) clipping exceptions
    if ylim is None:
        all_vals = np.concatenate([gas_arr, grain_arr])
        positive_vals = all_vals[all_vals > 0]
        if len(positive_vals) > 0:
            ax.set_ylim(positive_vals.min() * 0.5, positive_vals.max() * 2)
        else:
            ax.set_ylim(1e-4, 1e2)

    # Render grid lines and legendary context block
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(frameon=True, facecolor='white', edgecolor='gainsboro', loc='best')
    
    plt.tight_layout()
    plt.show()

import os
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from scipy.interpolate import griddata

def plot_grain_properties_midplane_comparison(PIPE,
                                              MODEL_NAMES,
                                              key_list=['CO'],
                                              itime=-1,
                                              fracab=True,
                                              verbose=True,
                                              xlim=None,
                                              ylim=None,
                                              Tmin=None,
                                              Tmax=None,
                                              temp_colormap='hot',
                                              ab_colormap='plasma',
                                              res_colormap='seismic',
                                              common_scale=True):
    r"""
    Plots and compares 2D maps of dust grain temperature (smoothly interpolated) 
    and total ice species properties (Surface J + Mantle K combined inside discrete blocks) 
    strictly at the disk midplane (z = 0) as a function of Disk Radius (R) and Grain Size (r).

    Dynamically manages structural subplots depending on the number of models in PIPE:
    - 1 Model: Displays Temperature followed by selected species inside a standard layout grid.
    - 2 Models: Displays a rows-based framework (Row 1: Temp, Rows 2+: Species) mapped across 
      a 3-column layout (Model 1, Model 2, and Model 1 - Model 2 Residuals).
    - >2 Models: Renders a matrix where columns represent models, Row 1 displays Temperatures, 
      and subsequent rows map the chemical species.

    Parameters:
    -----------
    PIPE : list
        Collection of model pipe objects to analyze. Each object must feature 
        `.chemistry` and `.chempath` attributes.
    MODEL_NAMES : list of str
        Display names assigned to the models for plot headers. If duplicate names 
        are found, a generic fallback nomenclature is generated.
    key_list : list of str, optional
        Target species formulas (without phase prefixes) to isolate and sum. Defaults to ['CO'].
    itime : int, default -1
        Time index to slice from the multi-epoch data arrays.
    fracab : bool, default True
        If True, plots fractional abundances ($n_{\rm ice}/n_H$). If False, plots 
        absolute volume number densities ($cm^{-3}$).
    verbose : bool, default True
        If True, outputs warning and missing file notifications to the console.
    xlim, ylim : tuple of float, optional
        Custom manual boundaries to enforce on the horizontal (Radius) and vertical (Bins) axes.
    Tmin, Tmax : float, optional
        Manual color mapping thresholds imposed over the thermal mesh layers.
    temp_colormap : str, default 'hot'
        Matplotlib colormap string used to render the continuous dust temperatures.
    ab_colormap : str, default 'plasma'
        Matplotlib colormap string used to style the discrete ice species blocks.
    res_colormap : str, default 'seismic'
        Matplotlib colormap string used to render the residual discrepancy panels.
    common_scale : bool, default True
        If True, unifies the colorbar scale ranges across all models for each species row.
    """
    # Enforce unique list elements for key_list
    if isinstance(key_list, str):
        key_list = [key_list]
    key_list = list(dict.fromkeys(key_list))

    # Apply global generic unique names fallback to prevent subplot label collisions
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)):
        MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    # --- INTERNAL HELPERS ---
    def parse_grain_sizes_midplane(file_path):
        """Extracts radii (um) and temperatures (K) for the midplane row from 1D_grain_sizes.in."""
        try:
            valid_lines = []
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('!'): continue
                    if '!' in line: line = line.split('!')[0].strip()
                    values = [float(val) for val in line.split()]
                    if values: valid_lines.append(values)
            if not valid_lines: return None, None
            
            midplane_values = valid_lines[-1]  # Extract last row corresponding to midplane z=0
            num_grains = len(midplane_values) // 4
            radii_cm = midplane_values[:num_grains]
            radii_um = [r * 10000.0 for r in radii_cm]
            temps_k = midplane_values[2 * num_grains : 3 * num_grains]
            return radii_um, temps_k
        except Exception as e:
            if verbose: print(f"Error parsing {file_path}: {e}")
            return None, None

    def clean_molec(mol_name):
        """Converts raw chemical strings into standard formatted LaTeX chemical expressions."""
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        return f"${f}$"

    # --- MULTI-MODEL DATABASE ACQUISITION ---
    model_data = {}
    time_years_string = None
    global_reference_sizes = None

    # Loop over pipelines
    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        main_output_dict = p.chemistry
        chempath = Path(p.chempath)
        
        # Filter directories with valid numerical spatial values
        radii_map = {}
        for original_key in main_output_dict.keys():
            digits = re.findall(r'\d+', str(original_key))
            if digits:
                try: radii_map[int(digits[0])] = original_key
                except ValueError: continue
                    
        extracted_radii = sorted(list(radii_map.keys()))
        
        disk_radii = []
        grain_sizes_matrix = []
        grain_temps_matrix = []
        species_abundance_matrices = {key: [] for key in key_list}

        # Parse file metadata columns
        for r_val in extracted_radii:
            grain_file = chempath / f"{r_val}AU" / "1D_grain_sizes.in"
            radii_um, temps_k = parse_grain_sizes_midplane(grain_file)
            
            if radii_um is None or temps_k is None: continue
                
            orig_key = radii_map[r_val]
            sub_dict = main_output_dict[orig_key]
            abundance_array = sub_dict['abundances']
            MIDPLANE_INDEX = -1
            
            disk_radii.append(r_val)
            grain_sizes_matrix.append(radii_um)
            grain_temps_matrix.append(temps_k)
            
            if global_reference_sizes is None:
                global_reference_sizes = radii_um
            
            if time_years_string is None:
                try:
                    t_sec = abundance_array.coords['time'].values[itime]
                    time_years_string = f"{t_sec / 3.156e7:.0f}"
                except: pass

            # Aggregate Surface (J) and Mantle (K) phases for target chemical keys
            for key in key_list:
                bin_values = []
                for b_idx in range(1, len(radii_um) + 1):
                    v_cell = 0.0
                    surface_name = f"J{b_idx:02d}{key}"
                    mantle_name  = f"K{b_idx:02d}{key}"
                    
                    if surface_name in abundance_array.coords['species'].values:
                        v_cell += float(abundance_array.isel(time=itime).sel(species=surface_name).values[MIDPLANE_INDEX])
                    if mantle_name in abundance_array.coords['species'].values:
                        v_cell += float(abundance_array.isel(time=itime).sel(species=mantle_name).values[MIDPLANE_INDEX])
                        
                    if not fracab:
                        v_cell *= sub_dict["H_number_density"][itime, MIDPLANE_INDEX]
                            
                    bin_values.append(v_cell)
                species_abundance_matrices[key].append(bin_values)

        if not disk_radii: continue

        model_data[p_name] = {
            'disk_radii': np.array(disk_radii),
            'grain_sizes': np.array(grain_sizes_matrix),
            'grain_temps': np.array(grain_temps_matrix),
            'abundance_matrices': species_abundance_matrices
        }

    # --- CANVAS GEOMETRY LAYOUT SETUP ---
    num_models = len(PIPE)
    model_names = list(model_data.keys())
    unit_str = "Abundance [$n_{\\rm ice}/n_H$]" if fracab else "Density [$cm^{-3}$]"

    if num_models == 1:
        num_plots = len(key_list) + 1
        cols = min(3, num_plots)
        rows = (num_plots + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 4.5 * rows), squeeze=False, sharex=True, sharey=True)
        axes = axes.flatten()
    else:
        rows = 1 + len(key_list)
        cols = 3 if num_models == 2 else num_models
        fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.0 * rows), squeeze=False, sharex=True, sharey=True)

    # Establish uniform grid structure layout indices
    p0 = model_names[0]
    num_grain_bins = model_data[p0]['grain_sizes'].shape[1]
    y_centers = np.arange(num_grain_bins)
    y_edges = np.arange(num_grain_bins + 1) - 0.5

    # --- RENDERING ENGINE LOOP ---
    # --- ROW 0 / PANEL 0: DUST GRAIN TEMPERATURES ---
    if common_scale:
        all_temps = np.concatenate([model_data[m]['grain_temps'].flatten() for m in model_names])
        actual_tmin = Tmin if Tmin is not None else all_temps.min()
        actual_tmax = Tmax if Tmax is not None else all_temps.max()
    else:
        actual_tmin, actual_tmax = Tmin, Tmax

    # 1. Plot Smoothly Interpolated Base Thermal Profiles
    for col_idx, p_name in enumerate(model_names):
        ax = axes[0, col_idx] if num_models > 1 else axes[0]
        struct = model_data[p_name]
        radii = struct['disk_radii']
        temps = struct['grain_temps']
        
        grid_R, grid_Y = np.meshgrid(np.linspace(radii.min(), radii.max(), 200), np.linspace(y_edges[0], y_edges[-1], 200))
        
        points_R, points_Y, points_T = [], [], []
        for i, r_val in enumerate(radii):
            for j, y_val in enumerate(y_centers):
                points_R.append(r_val)
                points_Y.append(y_val)
                points_T.append(temps[i, j])
            # Anchor blocks to improve mesh interpolation boundaries
            points_R.extend([r_val, r_val])
            points_Y.extend([y_edges[0], y_edges[-1]])
            points_T.extend([temps[i, 0], temps[i, -1]])
                
        grid_T = griddata((points_R, points_Y), points_T, (grid_R, grid_Y), method='cubic')
        tmin_local = actual_tmin if actual_tmin is not None else temps.min()
        tmax_local = actual_tmax if actual_tmax is not None else temps.max()
        
        cf = ax.contourf(grid_R, grid_Y, grid_T, levels=np.linspace(tmin_local, tmax_local, 50), cmap=temp_colormap)
        fig.colorbar(cf, ax=ax, label=r"$T_{\rm grain}$ [K]")
        ax.set_title(f"{p_name}\n" + r"Grain Temp $T_{\rm grain}$", fontsize=11, fontweight='bold')

    # 2. Temperature Residual Map Panel (Triggered only when exactly 2 models exist)
    if num_models == 2:
        ax_res = axes[0, 2]
        s1, s2 = model_data[model_names[0]], model_data[model_names[1]]
        
        # Interpolate onto a shared layout structure mesh before arithmetic evaluation
        r_min_shared = max(s1['disk_radii'].min(), s2['disk_radii'].min())
        r_max_shared = min(s1['disk_radii'].max(), s2['disk_radii'].max())
        grid_R_res, grid_Y_res = np.meshgrid(np.linspace(r_min_shared, r_max_shared, 200), np.linspace(y_edges[0], y_edges[-1], 200))
        
        g_T1 = griddata(([r for r in s1['disk_radii'] for _ in y_centers], [y for _ in s1['disk_radii'] for y in y_centers]), s1['grain_temps'].flatten(), (grid_R_res, grid_Y_res), method='linear')
        g_T2 = griddata(([r for r in s2['disk_radii'] for _ in y_centers], [y for _ in s2['disk_radii'] for y in y_centers]), s2['grain_temps'].flatten(), (grid_R_res, grid_Y_res), method='linear')
        
        res_temps = g_T1 - g_T2
        max_t_diff = np.nanmax(np.abs(res_temps)) if not np.all(np.isnan(res_temps)) else 1.0
        if max_t_diff == 0: max_t_diff = 1.0
        
        cf_res = ax_res.contourf(grid_R_res, grid_Y_res, res_temps, levels=np.linspace(-max_t_diff, max_t_diff, 50), cmap=res_colormap)
        fig.colorbar(cf_res, ax=ax_res, label=f"Diff $\\Delta T$ [K]")
        ax_res.set_title(f"Residuals Temp\n({model_names[0]} - {model_names[1]})", fontsize=11, fontweight='bold')

    # --- ROWS 1+: CHEMICAL ICE SPECIES DISCRETE PANELS ---
    for row_idx, key in enumerate(key_list):
        current_row = row_idx + 1 if num_models > 1 else row_idx + 1
        
        # Calculate shared scale thresholds across models for unified tracking rows
        if common_scale:
            all_row_vals = []
            for p_name in model_names:
                for r_idx in range(len(model_data[p_name]['disk_radii'])):
                    all_row_vals.extend(model_data[p_name]['abundance_matrices'][key][r_idx])
            all_row_vals = np.array(all_row_vals)
            pos_row_vals = all_row_vals[all_row_vals > 0]
            global_row_vmin = pos_row_vals.min() if len(pos_row_vals) > 0 else 1e-15
            global_row_vmax = all_row_vals.max()
            global_row_is_log = len(pos_row_vals) > 0 and (global_row_vmax / global_row_vmin) > 10.0

        polygons_per_model = {}
        values_per_model = {}

        # 1. Render Base Chemical Blocks
        for col_idx, p_name in enumerate(model_names):
            ax = axes[current_row, col_idx] if num_models > 1 else axes[current_row]
            struct = model_data[p_name]
            radii = struct['disk_radii']
            ab_matrix = struct['abundance_matrices'][key]
            
            if len(radii) > 1:
                r_midshifts = 0.5 * np.diff(radii)
                r_edges = [radii[0] - r_midshifts[0]] + [radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [radii[-1] + r_midshifts[-1]]
            else:
                r_edges = [radii[0] - 0.5, radii[0] + 0.5]
                
            polygons, values = [], []
            for i in range(len(radii)):
                r_left, r_right = r_edges[i], r_edges[i+1]
                for j in range(num_grain_bins):
                    z_bottom, z_top = y_edges[j], y_edges[j+1]
                    poly = [(r_left, z_bottom), (r_right, z_bottom), (r_right, z_top), (r_left, z_top)]
                    polygons.append(poly)
                    values.append(ab_matrix[i][j])
                    
            vals_array = np.array(values)
            polygons_per_model[p_name] = polygons
            values_per_model[p_name] = vals_array

            if common_scale:
                actual_vmin, actual_vmax = global_row_vmin, global_row_vmax
                is_log = global_row_is_log
            else:
                pos_vals = vals_array[vals_array > 0]
                actual_vmin = pos_vals.min() if len(pos_vals) > 0 else 1e-15
                actual_vmax = vals_array.max()
                is_log = len(pos_vals) > 0 and (actual_vmax / actual_vmin) > 10.0
                
            color_norm = plt.cm.colors.LogNorm(vmin=actual_vmin, vmax=actual_vmax) if is_log else plt.cm.colors.Normalize(vmin=actual_vmin, vmax=actual_vmax)
            
            coll = PolyCollection(polygons, array=vals_array, cmap=ab_colormap, norm=color_norm, edgecolors='none')
            ax.add_collection(coll)
            
            sm = plt.cm.ScalarMappable(cmap=ab_colormap, norm=color_norm)
            sm.set_array(vals_array)
            fig.colorbar(sm, ax=ax, label=unit_str)
            ax.set_title(f"{p_name}\n{clean_molec(key)} Ice Profile", fontsize=11)

        # 2. Chemical Species Residual Map Panel (Triggered only when num_models == 2)
        if num_models == 2:
            ax_res = axes[current_row, 2]
            p1_name, p2_name = model_names[0], model_names[1]
            
            # Linear arithmetic deduction on matching block grids
            res_vals = values_per_model[p1_name] - values_per_model[p2_name]
            max_diff = max(abs(res_vals.min()), abs(res_vals.max()))
            if max_diff == 0: max_diff = 1.0
            res_norm = plt.cm.colors.Normalize(vmin=-max_diff, vmax=max_diff)
            
            coll_res = PolyCollection(polygons_per_model[p1_name], array=res_vals, cmap=res_colormap, norm=res_norm, edgecolors='none')
            ax_res.add_collection(coll_res)
            
            sm_res = plt.cm.ScalarMappable(cmap=res_colormap, norm=res_norm)
            sm_res.set_array(res_vals)
            fig.colorbar(sm_res, ax=ax_res, label=f"Diff ({p1_name} - {p2_name})")
            ax_res.set_title(f"Residuals: {clean_molec(key)}\n({p1_name} - {p2_name})", fontsize=11)

    # --- GLOBAL FORMATTING & CLEANUP ---
    if global_reference_sizes is not None:
        tick_labels = [f"{size:.2f}" if size < 10 else f"{size:.1f}" for size in global_reference_sizes]
    else:
        tick_labels = [str(i) for i in y_centers]

    # Map geometric layouts constraints symmetrically across active axes panels
    for ax in axes.flatten():
        ax.set_yticks(y_centers)
        ax.set_yticklabels(tick_labels, fontsize=9)
        ax.set_xlabel('Disk Radius R [AU]', fontsize=10)
        ax.set_ylabel('Grain Radius $r$ [µm]', fontsize=10)
        ax.grid(True, linestyle=":", alpha=0.3)
        
        p_ref = model_data[model_names[0]]
        ax.set_xlim(xlim if xlim is not None else (p_ref['disk_radii'].min(), p_ref['disk_radii'].max()))
        ax.set_ylim(ylim if ylim is not None else (y_edges[0], y_edges[-1]))

    # Prune extra layout slots in single model configurations
    if num_models == 1:
        for idx in range(len(key_list) + 1, len(axes)):
            fig.delaxes(axes[idx])

    # Assign comprehensive descriptive header layout metadata strings
    if time_years_string:
        if num_models == 1 and len(key_list) == 0:
            axes[0].set_title(f"{axes[0].get_title()} \n $t = {time_years_string}$ years")
        else:
            fig.suptitle(f'Midplane ($z=0$) Grain-Gas Properties — $t = {time_years_string}$ years', fontsize=13, fontweight='bold', y=0.99)

    plt.tight_layout()
    plt.show()


import os
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import griddata

def plot_gas_fraction_map(chempath,
                          main_output_dict,
                          molecule="CO",
                          itime=-1,
                          threshold=50.0,
                          overlay_color='black',
                          verbose=True,
                          xlim=None,
                          ylim=None,
                          rmin=None,
                          rmax=None,
                          colormap="RdYlBu_r"):
    """
    Plots a smoothly interpolated 2D vertical cross-section map (R vs z) of the 
    gas-phase fraction of a chemical species relative to its total volumetric 
    budget (Gas + Surface J + Mantle K) at a specific simulation timestep.

    The local percentage is calculated using absolute particle counts integrated 
    over the cylindrical shell volume element (V = 2 * pi * R * dR * dz).
    An overlay is superimposed based on the `threshold` parameter format:
    - Scalar (float/int): Draws a solid line at that exact gas-phase percentage.
    - Sequence (list/tuple of 2 values): Regions a hatched transition zone.

    Parameters:
        chempath (str/Path): Path to the directory containing radial subfolders.
        main_output_dict (dict): Nested dictionary containing simulation outputs.
        molecule (str): Formula of the target chemical species (e.g., 'CO').
        itime (int): Index of the simulation timestep to plot. Defaults to -1.
        threshold (float or list/tuple): Percentage level(s) for the graphic overlay.
        overlay_color (str): Color of the contour line or hatched patterns.
        verbose (bool): If True, prints missing file or parsing diagnostics.
        xlim/ylim (tuple, optional): Custom limits for the R and z axes [AU].
        rmin/rmax (float, optional): Inbound radial filtering constraints [AU].
        colormap (str): Matplotlib colormap identifier. Defaults to 'RdYlBu_r'.
    """
    
    chempath = Path(chempath)
    AU_to_cm = 1.496e13  

    # --- Internal formatting helper ---
    def clean_molec(mol_name):
        f = re.sub(r"(\d+)", r"_{\1}", mol_name)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        return f"${f}$"

    # --- Radial key harvesting and filtering ---
    radii_map = {}
    for original_key in main_output_dict.keys():
        digits = re.findall(r'\d+', str(original_key))
        if digits:
            try:
                radii_map[int(digits[0])] = original_key
            except ValueError:
                continue
                
    extracted_radii = sorted(list(radii_map.keys()))
    
    sorted_radii = []
    for r in extracted_radii:
        if rmin is not None and r < rmin: continue
        if rmax is not None and r > rmax: continue
        sorted_radii.append(r)

    if not sorted_radii:
        if verbose: print(f"No columns found within specified boundaries [rmin={rmin}, rmax={rmax}].")
        return

    disk_radii = []
    max_z_encountered = 0.0

    grid_all_z = []
    scat_R = []
    scat_Z = []
    scat_Perc = []

    # --- Reconstruct radial grid cell edges ---
    if len(sorted_radii) > 1:
        r_midshifts = 0.5 * np.diff(sorted_radii)
        r_edges = [sorted_radii[0] - r_midshifts[0]] + [sorted_radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [sorted_radii[-1] + r_midshifts[-1]]
    else:
        r_edges = [sorted_radii[0] - 0.5, sorted_radii[0] + 0.5]

    # --- Parse physical and chemical profiles ---
    for i, r_val in enumerate(sorted_radii):
        file_path = chempath / f"{r_val}AU" / "1D_static.dat"
        if not os.path.exists(file_path):
            if verbose: print(f"File not found: {file_path}. Skipping radius.")
            continue
            
        orig_key = radii_map[r_val]
        sub_dict = main_output_dict[orig_key]
        abundance_array = sub_dict['abundances']
        nH_profile = sub_dict["H_number_density"][itime, :]  
        available_species = abundance_array.coords['species'].values
        
        try:
            z_points = np.loadtxt(file_path, comments='!', usecols=0)
        except Exception as e:
            if verbose: print(f"Error loading {file_path}: {e}")
            continue
            
        grid_all_z.extend(list(z_points))

        # Reconstruct vertical grid cell edges
        if len(z_points) > 1:
            z_midshifts = 0.5 * np.diff(z_points)
            z_edges = [z_points[0] - z_midshifts[0]] + [z_points[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_points[-1] + z_midshifts[-1])]
            dz_cells = np.abs(np.diff(z_edges))
        else:
            z_edges = [max(0.0, z_points[0] - 0.5), z_points[0] + 0.5]
            dz_cells = np.array([1.0])
            
        max_z_encountered = max(max_z_encountered, max(z_edges))
        disk_radii.append(r_val)
        
        # Determine cell radial width and center in cm
        r_left, r_right = r_edges[i], r_edges[i+1]
        dR = (r_right - r_left) * AU_to_cm
        R_center = float(r_val) * AU_to_cm
        
        # Identify gas and solid phase tokens for all size bins
        gas_name = molecule
        ice_names = []
        for sp in available_species:
            if sp.endswith(molecule) and (sp.startswith('J') or sp.startswith('K')):
                if re.match(r'^[JK]\d+' + re.escape(molecule) + r'$', sp):
                    ice_names.append(sp)

        y_abundances = abundance_array.isel(time=itime)

        # Calculate absolute gas populations integrated over cell volumes
        for j in range(len(z_points)):
            dz = dz_cells[j] * AU_to_cm
            cell_volume = 2 * np.pi * R_center * dR * dz
            nH_local = float(nH_profile[j])
            
            y_gas = float(y_abundances.sel(species=gas_name).values[j]) if gas_name in available_species else 0.0
            
            y_ice = 0.0
            for ice_sp in ice_names:
                y_ice += float(y_abundances.sel(species=ice_sp).values[j])
                
            n_abs_gas = y_gas * nH_local * cell_volume
            n_abs_ice = y_ice * nH_local * cell_volume
            total_particles = n_abs_gas + n_abs_ice
            
            if total_particles > 0.0:
                gas_percentage = (n_abs_gas / total_particles) * 100.0
            else:
                gas_percentage = 100.0

            scat_R.append(r_val)
            scat_Z.append(z_points[j])
            scat_Perc.append(gas_percentage)

    if not scat_R:
        if verbose: print("No spatial nodes successfully generated.")
        return

    scat_R = np.array(scat_R)
    scat_Z = np.array(scat_Z)
    scat_Perc = np.array(scat_Perc)
    disk_radii_arr = np.array(disk_radii)

    # --- Plot setup ---
    fig, ax = plt.subplots(figsize=(8, 4.5))

    # Generate a fine regular grid for 2D interpolation
    grid_R, grid_Z = np.meshgrid(
        np.linspace(disk_radii_arr.min(), disk_radii_arr.max(), 400),
        np.linspace(min(grid_all_z), max(grid_all_z), 400)
    )

    # Use linear interpolation for the background to suppress numerical overshoot artifacts
    grid_Percentage = griddata((scat_R, scat_Z), scat_Perc, (grid_R, grid_Z), method='linear')

    contour_levels = np.linspace(0.0, 100.0, 101)
    color_norm = plt.cm.colors.Normalize(vmin=0.0, vmax=100.0)

    # Render background color fields
    cf = ax.contourf(grid_R, grid_Z, grid_Percentage, levels=contour_levels, cmap=colormap, norm=color_norm, extend='both')
    
    cb = fig.colorbar(cf, ax=ax, label=f"Gas-phase Fraction [ % of Total {clean_molec(molecule)} (Gas+J+K) ]", 
                      ticks=np.linspace(0, 100, 11), extendfrac=0)

    # Use cubic interpolation strictly for smoother contour lines
    grid_Percentage_smooth = griddata((scat_R, scat_Z), scat_Perc, (grid_R, grid_Z), method='cubic')

    # --- Render overlays and legends ---
    try:
        # Case 1: Sequence threshold passed -> Hatched region implementation
        if isinstance(threshold, (list, tuple)) and len(threshold) == 2:
            t_low = float(min(threshold))
            t_high = float(max(threshold))
            
            cs_zone = ax.contourf(grid_R, grid_Z, grid_Percentage_smooth, 
                                  levels=[t_low, t_high], 
                                  colors='none',      
                                  hatches=['////'])   
            
            # Version-agnostic loop to set custom hatch edge color
            if hasattr(cs_zone, 'collections') and cs_zone.collections:
                target_collections = cs_zone.collections
            elif hasattr(cs_zone, 'artists') and cs_zone.artists:
                target_collections = cs_zone.artists
            else:
                target_collections = [ax.collections[-1]]
                
            for col in target_collections:
                col.set_edgecolor(overlay_color)
                col.set_linewidth(0.5) 
            
            hatch_patch = plt.Rectangle((0, 0), 1, 1, fill=False, 
                                        edgecolor=overlay_color, 
                                        hatch='////', linewidth=0.5)
            
            ax.legend(handles=[hatch_patch], 
                      labels=[f'Transition Zone ({t_low:.0f}% - {t_high:.0f}% Gas)'], 
                      loc='upper right', frameon=True, fontsize=10)
            
        # Case 2: Scalar threshold passed -> Solid line implementation
        else:
            t_val = float(threshold)
            cs_line = ax.contour(grid_R, grid_Z, grid_Percentage_smooth, 
                                 levels=[t_val], 
                                 colors=overlay_color, 
                                 linestyles='solid', 
                                 linewidths=1.7)
            
            if hasattr(cs_line, 'legend_elements'):
                artists, _ = cs_line.legend_elements()
                h_item = artists[0]
            elif hasattr(cs_line, 'collections') and cs_line.collections:
                h_item = cs_line.collections[0]
            else:
                h_item = ax.collections[-1]
                
            ax.legend(handles=[h_item], 
                      labels=[f'Snowline ({t_val:.0f}% Gas Fraction)'], 
                      loc='upper right', frameon=True, fontsize=10)
                
    except Exception as e:
        if verbose: print(f"Warning: Could not render threshold overlay structures. {e}")

    # --- Axes properties and bounds ---
    ax.set_xlabel('Disk Radius R [AU]', fontsize=11)
    ax.set_ylabel('Altitude $z$ [AU]', fontsize=11)
    ax.set_title(f"Gas-phase Volumetric Reservoir Distribution Map — {clean_molec(molecule)}", fontsize=12, fontweight='bold', pad=10)

    ax.set_xlim(xlim if xlim is not None else (disk_radii_arr.min(), disk_radii_arr.max()))
    ax.set_ylim(ylim if ylim is not None else (0.0, max_z_encountered))
    ax.grid(True, linestyle=":", alpha=0.3)
    ax.tick_params(labelsize=11)

    try:
        sample_orig_key = radii_map[disk_radii[0]]
        time_seconds = main_output_dict[sample_orig_key]['abundances'].coords['time'].values[itime]
        fig.suptitle(f'Vertical Cross-section Curve — $t = {time_seconds/3.156e7:.0f}$ years', fontsize=11, fontweight='bold', y=0.99)
    except:
        pass

                          
    plt.tight_layout()
    plt.show()


import os 
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm, Normalize 
from pathlib import Path 

def density2D_grid_comparison(PIPE, 
                             MODEL_NAMES, 
                             vmin=1e-30, 
                             vmax=1e-15, 
                             cmap='gnuplot2', 
                             res_colormap='seismic', 
                             dens_type='mass', 
                             xlim=None, 
                             ylim=None, 
                             dust=None, 
                             select_bins=None,  
                             figsize=(14, 16)): 
    r""" 
    Plots and compares 2D poloidal gas/dust density grids (R, Z) across multiple  
    simulation models, dynamically mapping dust size bins and total accumulations. 
    """ 
    # Prevent layout collisions by generating generic fallback names if duplicates exist 
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)): 
        MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))] 

    autocm = 1.495978707e13  # Astronomical Unit converted to centimeters 
    model_structures = {} 

    # --- 1. CORE EXTRACTION ENGINE CROSS PIPELINES --- 
    for p_idx, p in enumerate(PIPE): 
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx]) 
        base_path = str(Path(p.thermalpath)) + '/' 
         
        # Load and parse spatial grid configurations 
        grid_file = base_path + 'amr_grid.inp' 
        if not os.path.exists(grid_file): 
            print(f"[{p_name}] Missing critical amr grid file: {grid_file}. Skipping.") 
            continue 
             
        grid = pd.read_table(grid_file, engine='python', skiprows=5) 
        nr = int(grid.columns[0].split()[0]) 
        nt = int(grid.columns[0].split()[1]) 
        grid_vals = grid[grid.columns[0]].values 

        # Load and parse simulation density profiles 
        dens_file = base_path + 'dust_density.inp' 
        if not os.path.exists(dens_file): 
            print(f"[{p_name}] Missing dust density dataset file: {dens_file}. Skipping.") 
            continue 
             
        dens = pd.read_table(dens_file, engine='python', header=None, skiprows=3) 
        dens_vals = dens[0].values 
        nspecies = int(len(dens_vals) / (nr * nt)) 
        dens_reshaped = np.reshape(dens_vals, (nspecies, nt, nr)) 

        # Reconstruct 2D cylindrical coordinates mapping mesh (Radius, Altitude) 
        r_edge = grid_vals[:nr+1] / autocm 
        theta_edge = grid_vals[nr+1:nr+1+nt+1].copy()  
        theta_edge[-1] = np.pi 
        rr_edge, tt_edge = np.meshgrid(r_edge, theta_edge) 
        R = rr_edge * np.sin(tt_edge) 
        Z = rr_edge * np.cos(tt_edge) 
         
        # Sanitize zero values to protect log-scale normalization structures 
        dens_clean = np.array(dens_reshaped, copy=True) 
        dens_clean[dens_clean <= 1e-100] = 1e-100 

        # Retrieve optional physical grain sizes arrays definitions 
        sizes_file = base_path + 'dust_sizes.inp' 
        grain_mass = None 
        if os.path.isfile(sizes_file): 
            sizes = np.loadtxt(sizes_file) 
            sizes = np.atleast_1d(sizes) 
        elif dust is not None: 
            sizes = dust.sizes()[0] 
            grain_mass = dust.grainmass() 
        else: 
            sizes = None 
            print(f"⚠️ [{p_name}] WARNING: No dust size source found (checked '{sizes_file}' and 'dust' argument).") 

        # --- BIN SELECTION FILTER --- 
        if select_bins is not None: 
            # Guard constraint to prevent out-of-bounds indexing errors 
            valid_bins = [b for b in select_bins if b < nspecies] 
            dens_clean = dens_clean[valid_bins] 
            if sizes is not None: 
                sizes = sizes[valid_bins] 
            if grain_mass is not None: 
                grain_mass = grain_mass[valid_bins] 
            nspecies_to_display = len(valid_bins) 
        else: 
            nspecies_to_display = nspecies 
            valid_bins = list(range(nspecies)) 

        model_structures[p_name] = { 
            'R': R, 'Z': Z, 'dens': dens_clean, 
            'sizes': sizes, 'grain_mass': grain_mass, 
            'nspecies_display': nspecies_to_display,  
            'original_bins': valid_bins, 
            'nr': nr, 'nt': nt 
        } 

    if not model_structures: 
        print("No simulation pipelines successfully parsed.") 
        return 

    # --- 2. CANVAS GEOMETRY CONFIGURATION --- 
    num_models = len(model_structures) 
    model_names = list(model_structures.keys()) 
     
    nspecies_ref = model_structures[model_names[0]]['nspecies_display'] 
    original_bins_ref = model_structures[model_names[0]]['original_bins'] 

    if num_models == 1: 
        npanels = nspecies_ref + 1 
        ncols = min(nspecies_ref, 4) 
        nrows = int(np.ceil(npanels / ncols)) 
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True) 
        axes = np.atleast_2d(axes) 
        cbar_main_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7]) 
    else: 
        nrows = nspecies_ref + 1  
        cols = 3 if num_models == 2 else num_models 
        fig, axes = plt.subplots(nrows, cols, figsize=figsize, sharex=True, sharey=True) 
        axes = np.atleast_2d(axes) 
         
        if num_models == 2: 
            cbar_main_ax = fig.add_axes([0.91, 0.55, 0.015, 0.35]) 
            cbar_res_ax  = fig.add_axes([0.91, 0.15, 0.015, 0.35]) 
        else: 
            cbar_main_ax = fig.add_axes([0.91, 0.15, 0.015, 0.7]) 

    # --- 3. RENDERING ENGINE LOOP --- 
    im, im_res = None, None 
    label_text = "" 

    for row_idx in range(nspecies_ref + 1): 
        for col_idx, p_name in enumerate(model_names): 
            ax = axes[row_idx, col_idx] if num_models > 1 else axes.flat[row_idx] 
            struct = model_structures[p_name] 
            R, Z = struct['R'], struct['Z'] 
            sizes, grain_mass = struct['sizes'], struct['grain_mass'] 
             
            # Map specific grain sizes bins or total calculated matrix sums 
            if row_idx < nspecies_ref: 
                data_matrix = struct['dens'][row_idx] 
                actual_bin_num = original_bins_ref[row_idx] + 1 
                 
                # --- SUBPLOT TITLE GENERATION WITH PHYSICAL SIZE --- 
                if sizes is not None and row_idx < len(sizes): 
                    s = sizes[row_idx] 
                    size_str = f'{s/1e3:.1f} mm' if s >= 1e3 else f'{s:.2f} µm' 
                    title_string = f"bin {actual_bin_num} ({size_str})" 
                else: 
                    title_string = f"bin {actual_bin_num}" 
                 
                if dens_type == 'number' and grain_mass is not None: 
                    plot_data = data_matrix / grain_mass[row_idx] 
                    label_text = r'n$_\mathrm{d}$ [cm$^{-3}$]' 
                elif dens_type == 'surface' and grain_mass is not None and sizes is not None: 
                    plot_data = 4 * np.pi * (sizes[row_idx] * 1e-4)**2 * data_matrix / grain_mass[row_idx] 
                    label_text = r'surfaces [cm$^{-1}$]' 
                else: 
                    plot_data = data_matrix 
                    label_text = r'$\rho_\mathrm{d}$ [g cm$^{-3}$]' 
                    dens_type = 'mass' 
                     
            else: 
                plot_data = struct['dens'].sum(axis=0) 
                title_string = f"Total {dens_type}" 
                label_text = r'$\rho_\mathrm{d}$ [g cm$^{-3}$]' 

            # Render current 2D colormesh panel 
            im = ax.pcolormesh(R, Z, plot_data, cmap=cmap, shading='auto', norm=LogNorm(vmin=vmin, vmax=vmax)) 
            ax.set_title(title_string, fontsize=11) 
             
            # Label model columns headers along the upper row layout border 
            if num_models > 1 and row_idx == 0: 
                ax.text(0.5, 1.25, p_name, transform=ax.transAxes, fontsize=12,  
                        fontweight='bold', ha='center', va='bottom') 

        # 4. Residuals Scenarios (Triggered only when exactly 2 models are present) 
        if num_models == 2: 
            ax_res = axes[row_idx, 2] 
            s1, s2 = model_structures[model_names[0]], model_structures[model_names[1]] 
             
            if row_idx < nspecies_ref: 
                d1 = s1['dens'][row_idx] 
                d2 = s2['dens'][row_idx] 
                actual_bin_num = original_bins_ref[row_idx] + 1 
                 
                # Apply size mappings to the residuals axis subplot titles if available 
                if s1['sizes'] is not None and row_idx < len(s1['sizes']): 
                    s = s1['sizes'][row_idx] 
                    size_str = f'{s/1e3:.1f} mm' if s >= 1e3 else f'{s:.2f} µm' 
                    res_title = f"Residuals bin {actual_bin_num} ({size_str})" 
                else: 
                    res_title = f"Residuals bin {actual_bin_num}" 
            else: 
                d1 = s1['dens'].sum(axis=0) 
                d2 = s2['dens'].sum(axis=0) 
                res_title = "Residuals total" 

            # Compute linear discrepancy profile differences matrix 
            res_data = d1 - d2 
            max_diff = max(abs(res_data.min()), abs(res_data.max())) 
            if max_diff == 0: max_diff = 1.0 
             
            im_res = ax_res.pcolormesh(s1['R'], s1['Z'], res_data, cmap=res_colormap, shading='auto', 
                                       norm=Normalize(vmin=-max_diff, vmax=max_diff)) 
            ax_res.set_title(res_title, fontsize=10, fontweight='bold') 
             
            if row_idx == 0: 
                ax_res.text(0.5, 1.25, f"Diff\n({model_names[0]} - {model_names[1]})", transform=ax_res.transAxes,  
                            fontsize=11, fontweight='bold', ha='center', va='bottom') 

    # --- 4. COLORBARS & SPACING TUNING --- 
    fig.colorbar(im, cax=cbar_main_ax, label=label_text) 
    if num_models == 2 and im_res is not None: 
        fig.colorbar(im_res, cax=cbar_res_ax, label="Linear Discrepancy Map ($M_1 - M_2$)") 

    # Prune unassigned subplot windows in single model configurations 
    if num_models == 1: 
        for idx in range(npanels, nrows * ncols): 
            fig.delaxes(axes.flat[idx]) 

    # Enforce global formatting boundaries limits and structural grids 
    for ax in axes.flat: 
        if not ax.get_visible(): continue 
        if xlim: ax.set_xlim(xlim) 
        if ylim: ax.set_ylim(ylim) 
        ax.grid(True, linestyle=":", alpha=0.3) 
         
    for ax in axes[-1, :]: 
        if ax.get_visible(): ax.set_xlabel('r [au]', fontsize=12) 
    for ax in axes[:, 0]: 
        if ax.get_visible(): ax.set_ylabel('z [au]', fontsize=12) 

    # Adjust vertical hspace to prevent subtitle text overlap collisions across multi-row matrix meshes 
    fig.subplots_adjust(right=0.88, left=0.08, bottom=0.06, top=0.92, hspace=0.40, wspace=0.15) 

    plt.show()


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from pathlib import Path

def temperature2D_grid_comparison(PIPE,
                                 MODEL_NAMES,
                                 vmin=1.0,
                                 vmax=1e3,
                                 cmap='gnuplot2',
                                 res_colormap='seismic',
                                 xlim=None,
                                 ylim=None,
                                 snowline_temp=20.0,
                                 figsize=(14, 12)):
    r"""
    Plots and compares 2D poloidal dust temperature grids (R, Z) across multiple 
    RADMC-3D simulation models, dynamically mapping individual dust size grain profiles.

    This routine extracts the localized thermal structure of circumstellar disks 
    modeled via Monte Carlo radiative transfer. It handles spatial coordinate 
    transformations from native spherical/AMR frameworks into poloidal cartesian 
    mapping slices, applies customizable snowline tracking contours, and builds 
    a multi-panel canvas for cross-model comparative thermodynamics.

    Subplot Layout Engine
    ---------------------
    - 1 Model: Creates a standard optimized panel grid distributing every explicit 
      size bin sequence dynamically (max 4 columns per row block).
    - 2 Models: Forces a strict 3-column matrix per size bin: 
      [Model 1 Absolute T] [Model 2 Absolute T] [Linear Thermodynamic Delta (M1 - M2)]
    - >2 Models: Allocates a multi-column strict canvas where columns isolate 
      individual simulation paths and rows step through dust size populations.

    Parameters
    ----------
    PIPE : list
        Collection of user-defined model pipeline objects. Each object must 
        possess a `.chempath` or `.thermalpath` attribute pointing to its 
        respective RADMC-3D raw file repository.
    MODEL_NAMES : list of str
        Display names assigned to the simulation models for canvas sub-headers. 
        Generic fallbacks ("Model 1", "Model 2") are enforced if duplicates exist.
    vmin, vmax : float, optional
        Logarithmic normalization bounds applied to the absolute temperature profiles [K].
        Defaults to vmin=1.0, vmax=1000.0.
    cmap : str, optional
        Matplotlib sequential colormap code used to style absolute temperature fields. 
        Defaults to 'gnuplot2'.
    res_colormap : str, optional
        Matplotlib divergent colormap code used to style residual absolute delta frames. 
        Defaults to 'seismic'.
    xlim, ylim : tuple of float, optional
        Cylindrical boundaries (min, max) in Astronomical Units [au] used to crop 
        the radial ($R$) and vertical ($Z$) axes windows.
    snowline_temp : float, optional
        Target thermodynamic threshold temperature in Kelvin [K] to track across the 
        disk mesh (e.g., 20 K for CO sublimation, 150 K for H2O sublimation). 
        Set to None to deactivate contour generation. Defaults to 20.0.
    figsize : tuple of float, optional
        Total width and height geometry properties passed to initialize the subplots. 
        Defaults to (14, 12).

    Returns
    -------
    None
        Assembles, formats, and renders an interactive matplotlib canvas window.
    """
    # Force generic names if duplicate entries are found in MODEL_NAMES
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)):
        MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    autocm = 1.495978707e13  # Astronomical Unit conversion factor [cm/au]
    model_structures = {}

    # =========================================================================
    # --- 1. CORE DATA ACQUISITION & COORDINATE TRANSFORMATION ENGINE ---
    # =========================================================================
    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        base_path = str(Path(p.thermalpath)) + '/'
        
        # Verify existence of critical RADMC-3D configuration grid file
        grid_file = base_path + 'amr_grid.inp'
        if not os.path.exists(grid_file):
            print(f"[{p_name}] Critical error: AMR grid file not found at {grid_file}. Skipping pipeline.")
            continue
            
        # Parse spatial spherical grid dimensions
        grid = pd.read_table(grid_file, engine='python', skiprows=5)
        nr = int(grid.columns[0].split()[0])  # Radial cell counts
        nt = int(grid.columns[0].split()[1])  # Poloidal angular cell counts
        grid_vals = grid[grid.columns[0]].values

        # Verify existence of Monte Carlo dust temperature data arrays
        temp_file = base_path + 'dust_temperature.dat'
        if not os.path.exists(temp_file):
            print(f"[{p_name}] Critical error: Dust temperature array missing at {temp_file}. Skipping pipeline.")
            continue
            
        # Extract flat 1D raw thermal data stream
        temp = pd.read_table(temp_file, engine='python', header=None, skiprows=3)
        temp_vals = temp[0].values
        
        # Deduce total explicit dust grain size populations (bins)
        nspecies = int(len(temp_vals) / (nr * nt))
        
        # Reshape data array into 3D tensor: (grain_bin, theta_index, radius_index)
        temp_reshaped = np.reshape(temp_vals, (nspecies, nt, nr))

        # Reconstruct spatial cell boundaries from grid definitions
        r_edge = grid_vals[:nr+1] / autocm  # Convert cm to au
        
        # Isolate angular bounds safely inside a write-enabled copy array
        theta_edge = grid_vals[nr+1:nr+1+nt+1].copy()
        theta_edge[-1] = np.pi  # Force strict boundary limits to eliminate floating noise
        
        # Synthesize 2D boundary grid matrices
        rr_edge, tt_edge = np.meshgrid(r_edge, theta_edge)
        
        # Map spherical boundaries grids directly onto cartesian poloidal slices (R, Z)
        R = rr_edge * np.sin(tt_edge)
        Z = rr_edge * np.cos(tt_edge)

        # Interpolate boundary fields into true cell centers to map unbiased contours
        R_center = 0.25 * (R[:-1, :-1] + R[1:, :-1] + R[:-1, 1:] + R[1:, 1:])
        Z_center = 0.25 * (Z[:-1, :-1] + Z[1:, :-1] + Z[:-1, 1:] + Z[1:, 1:])

        # Retrieve optional physical grain sizes configurations
        sizes_file = base_path + 'dust_sizes.inp'
        if os.path.isfile(sizes_file):
            sizes = np.loadtxt(sizes_file)
            sizes = np.atleast_1d(sizes)
        else:
            sizes = None

        # Archive compiled geometric data blocks per unique pipeline path
        model_structures[p_name] = {
            'R': R, 'Z': Z, 'R_center': R_center, 'Z_center': Z_center,
            'temp': temp_reshaped, 'sizes': sizes, 'nspecies': nspecies
        }

    if not model_structures:
        print("Pipeline parsing sequence terminated. Zero valid models compiled.")
        return

    # =========================================================================
    # --- 2. CANVAS GEOMETRY LAYOUT SETUP ---
    # =========================================================================
    num_models = len(model_structures)
    model_names = list(model_structures.keys())
    nspecies_ref = model_structures[model_names[0]]['nspecies']
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    if num_models == 1:
        ncols = min(nspecies_ref, 4)
        nrows = int(np.ceil(nspecies_ref / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes = np.atleast_2d(axes)
        cbar_main_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])
    else:
        nrows = nspecies_ref
        cols = 3 if num_models == 2 else num_models
        fig, axes = plt.subplots(nrows, cols, figsize=figsize)
        axes = np.atleast_2d(axes)
        
        if num_models == 2:
            cbar_main_ax = fig.add_axes([0.91, 0.55, 0.015, 0.35])
            cbar_res_ax  = fig.add_axes([0.91, 0.15, 0.015, 0.35])
        else:
            cbar_main_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])

    # =========================================================================
    # --- 3. RENDERING ENGINE LOOP ---
    # =========================================================================
    im, im_res = None, None

    # Step sequentially through grain size populations (Rows)
    for row_idx in range(nspecies_ref):
        # Step sequentially through active simulation model streams (Columns)
        for col_idx, p_name in enumerate(model_names):
            ax = axes[row_idx, col_idx] if num_models > 1 else axes.flat[row_idx]
            struct = model_structures[p_name]
            R, Z = struct['R'], struct['Z']
            R_c, Z_c = struct['R_center'], struct['Z_center']
            sizes = struct['sizes']
            
            plot_data = struct['temp'][row_idx]
            title_string = f"bin {row_idx + 1}"

            # Render absolute thermodynamic fields onto boundary meshes
            im = ax.pcolormesh(R, Z, plot_data, cmap=cmap, shading='auto', norm=LogNorm(vmin=vmin, vmax=vmax))
            ax.set_title(title_string, fontsize=11)
            
            # Extract and overlay localized phase condensation contours
            if snowline_temp is not None:
                try:
                    cs = ax.contour(R_c, Z_c, plot_data, levels=[float(snowline_temp)], 
                                    colors='black', linewidths=1.6, linestyles='-')
                    if col_idx == 0:
                        ax.clabel(cs, fmt=f"{snowline_temp:.0f} K", fontsize=8, inline=True)
                except:
                    pass  # Gracefully ignore if isotherm boundaries lie outside current limits
            
            # Map un-cluttered pipeline subheaders onto the upper layer row boundaries
            if num_models > 1 and row_idx == 0:
                ax.text(0.5, 1.20, p_name, transform=ax.transAxes, fontsize=12, 
                        fontweight='bold', ha='center', va='bottom')
            
            # Display grain size annotations inside wheat bounding boxes
            if sizes is not None and row_idx < len(sizes):
                s = sizes[row_idx]
                size_label = f'{s/1e3:.1f} mm' if s >= 1e3 else f'{s:.2f} ' + r'$\mu$m'
                ax.text(0.05, 0.85, size_label, transform=ax.transAxes, fontsize=9, verticalalignment='top', bbox=props)

        # ---------------------------------------------------------------------
        # --- 4. DUAL-MODEL THERMODYNAMIC DELTA RESIDUAL ANALYSIS ---
        # ---------------------------------------------------------------------
        if num_models == 2:
            ax_res = axes[row_idx, 2]
            s1, s2 = model_structures[model_names[0]], model_structures[model_names[1]]
            
            t1 = s1['temp'][row_idx]
            t2 = s2['temp'][row_idx]
            res_title = f"Residuals bin {row_idx+1}"

            # Calculate absolute arithmetic discrepancy delta values (\Delta T = T_1 - T_2)
            res_data = t1 - t2
            max_diff = max(abs(res_data.min()), abs(res_data.max()))
            if max_diff == 0: max_diff = 1.0  
            
            # Render divergent panel using a zero-centered normalization mapping profile
            im_res = ax_res.pcolormesh(s1['R'], s1['Z'], res_data, cmap=res_colormap, shading='auto',
                                       norm=Normalize(vmin=-max_diff, vmax=max_diff))
            ax_res.set_title(res_title, fontsize=10, fontweight='bold')
            
            if row_idx == 0:
                ax_res.text(0.5, 1.20, f"Diff\n({model_names[0]} - {model_names[1]})", transform=ax_res.transAxes, 
                            fontsize=11, fontweight='bold', ha='center', va='bottom')

    # =========================================================================
    # --- 5. COLORBARS & SPATIAL WINDOW POST-PROCESSING ---
    # =========================================================================
    fig.colorbar(im, cax=cbar_main_ax, label=r'Dust Temperature $T_\mathrm{d}$ [K]')
    
    if num_models == 2 and im_res is not None:
        fig.colorbar(im_res, cax=cbar_res_ax, label=r'Thermal Discrepancy $\Delta T$ [K]')

    # Clean unoccupied canvas panels layout in unique standalone model tracking modes
    if num_models == 1:
        for idx in range(nspecies_ref, nrows * ncols):
            fig.delaxes(axes.flat[idx])

    # Enforce crop boundaries limits and share axis labeling context symmetrically
    for r_idx in range(nrows):
        for c_idx in range(cols if num_models > 1 else ncols):
            ax = axes[r_idx, c_idx] if num_models > 1 else axes[r_idx, c_idx]
            if not ax.get_visible(): continue
            
            if xlim: ax.set_xlim(xlim)
            if ylim: ax.set_ylim(ylim)
            ax.grid(True, linestyle=":", alpha=0.3)
            
            # Protect inner panel layout context from vertical/horizontal string labels squeezing
            if r_idx == nrows - 1:
                ax.set_xlabel('r [au]', fontsize=11)
            if c_idx == 0:
                ax.set_ylabel('z [au]', fontsize=11)

    # Tight margin padding allocations block tuning
    fig.subplots_adjust(right=0.88, left=0.08, bottom=0.08, top=0.90, hspace=0.35, wspace=0.22)

    plt.show()
