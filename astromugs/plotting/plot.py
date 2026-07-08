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

import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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
        ax.plot(radii, midtemp[ispec].rolling(window=2, center=True).mean(), linewidth=2, linestyle='-', label='bin: {}'.format(ispec+1))
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
