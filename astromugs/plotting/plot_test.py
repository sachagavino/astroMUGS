#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
File name: plot_test.py
Author: Sacha Gavino (Updated)
Language: PYTHON 3
Description: Comprehensive plotting suite for disk thermal models and chemical evolution structures (RADMC-3D & Nautilus).
"""

import glob
import sys
import os
import re
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata

# Fallback global constant if astromugs is not installed locally
try:
    from astromugs.constants.constants import autocm
except ImportError:
    autocm = 1.495978707e13  # 1 AU in cm


def density2D_grid(path='thermal/', vmin=1e-30, vmax=1e-15, cmap='gnuplot2', dens_type='mass',
                    xlim=None, ylim=None, dust=None, figsize=(10, 14)):
    """
    Plots 2D poloidal views of all dust species density distributions from a single model run.
    """
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
    dens = np.array(dens, copy=True)
    dens[dens <= 1e-100] = 1e-100

    sizes_file = path + 'dust_sizes.inp'
    if os.path.isfile(sizes_file):
        sizes = np.loadtxt(sizes_file)
        sizes = np.atleast_1d(sizes)
    elif dust is not None:
        rho_m = dust.rho_m
        sizes = dust.sizes()[0]
        grain_mass = dust.grainmass()
    else:
        sizes = None

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
                size_label = f'{s/1e3:.1f} mm' if s >= 1e3 else f'{s:.2f} ' + r'$\mu$m'
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

    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel('r [au]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel('z [au]', fontsize=14)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    plt.show()


def density1D_midplane(path='thermal/', vmin=1e-30, vmax=1e-15, dens_type='mass',
                       xlim=None, dust=None, figsize=(12, 8)):
    """
    Plots the 1D dust density profile in the midplane as a function of radius for each dust species.
    """
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = np.array(grid[grid.columns[0]].values, copy=True)
    
    dens = pd.read_table(path + 'dust_density.inp', engine='python', header=None, skiprows=3)
    dens = dens[0].values
    
    nspecies = int(len(dens) / (nr * nt))
    dens = np.reshape(dens, (nspecies, nt, nr))

    r_edge = grid[:nr+1] / autocm
    r_center = 0.5 * (r_edge[:-1] + r_edge[1:])
    idx_midplane = nt // 2 

    sizes_file = path + 'dust_sizes.inp'
    if os.path.isfile(sizes_file):
        sizes = np.loadtxt(sizes_file)
        sizes = np.atleast_1d(sizes)
    elif dust is not None:
        rho_m = dust.rho_m
        sizes = dust.sizes()[0]
        grain_mass = dust.grainmass()
    else:
        sizes = None

    npanels = nspecies + 1 
    ncols = min(npanels, 3)
    nrows = int(np.ceil(npanels / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)

    if dens_type == 'number':
        ylabel = r'$n_\mathrm{d}$ [cm$^{-3}$]'
    elif dens_type == 'surface':
        ylabel = r'Surfaces [cm$^{-1}$]'
    else:
        ylabel = r'$\rho_\mathrm{d}$ [g cm$^{-3}$]'

    for idx in range(nrows * ncols):
        ax = axes.flat[idx]
        
        if idx < nspecies:
            profile = dens[idx, idx_midplane, :]
            
            if dens_type == 'number':
                y_data = profile / grain_mass[idx]
            elif dens_type == 'surface':
                y_data = 4 * np.pi * (sizes[idx] * 1e-4) * profile / grain_mass[idx]
            elif dens_type == 'mass':
                y_data = profile
                
            ax.plot(r_center, y_data, color='darkblue', lw=2)
            ax.set_title(f'Bin {idx+1}', fontsize=12)
            
            if sizes is not None and idx < len(sizes):
                s = sizes[idx]
                size_label = f'{s/1e3:.1f} mm' if s >= 1e3 else f'{s:.2f} ' + r'$\mu$m'
                ax.text(0.05, 0.95, size_label, transform=ax.transAxes,
                        fontsize=12, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        elif idx == nspecies:
            if dens_type == 'mass':
                total_profile = dens[:, idx_midplane, :].sum(axis=0)
                ax.plot(r_center, total_profile, color='black', lw=2.5, linestyle='--')
                ax.set_title('Total Mass', fontsize=12)
            else:
                ax.axis('off')
        else:
            ax.axis('off')

        ax.set_yscale('log')
        ax.set_ylim(vmin, vmax)
        if xlim:
            ax.set_xlim(xlim)
        else:
            ax.set_xscale('log')

    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel('r [au]', fontsize=12)
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel, fontsize=12)

    fig.tight_layout()
    plt.show()


def density2D_grid_interactive(path='thermal/', vmin=1e-30, vmax=1e-15, cmap='gnuplot2', dens_type='mass',
                                xlim=None, ylim=None, dust=None, figsize=(10, 14)):
    """
    Interactive version of density2D_grid with live sliders for vmin/vmax limits optimization.
    Requires %matplotlib widget in the notebook.
    """
    import ipywidgets as widgets
    from IPython.display import display

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

    plot_data = []
    for idx in range(nspecies):
        if dens_type == 'number' and dust is not None:
            plot_data.append(4 * np.pi * sizes[idx] * 1e-4 * dens[idx] / grain_mass[idx])
        else:
            plot_data.append(dens[idx])
    plot_data.append(dens.sum(axis=0))

    npanels = nspecies + 1
    ncols = min(nspecies, 4)
    nrows = int(np.ceil(npanels / ncols))

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
                size_label = f'{s/1e3:.1f} mm' if s >= 1e3 else f'{s:.2f} ' + r'$\mu$m'
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
        if xlim: ax.set_xlim(xlim)
        if ylim: ax.set_ylim(ylim)

    for ax in axes[-1, :]:
        if ax.get_visible(): ax.set_xlabel('r [au]', fontsize=14)
    for ax in axes[:, 0]: ax.set_ylabel('z [au]', fontsize=14)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(meshes[-1], cax=cbar_ax, label=r'$\rho_\mathrm{d}$ [g cm$^{-3}$]')

    log_vmin, log_vmax = np.log10(vmin), np.log10(vmax)
    vmin_slider = widgets.FloatSlider(value=log_vmin, min=-50, max=0, step=0.5, description='log(vmin)', continuous_update=False, layout=widgets.Layout(width='400px'))
    vmax_slider = widgets.FloatSlider(value=log_vmax, min=-50, max=0, step=0.5, description='log(vmax)', continuous_update=False, layout=widgets.Layout(width='400px'))

    def update_clim(change):
        new_vmin, new_vmax = 10**vmin_slider.value, 10**vmax_slider.value
        if new_vmin >= new_vmax: return
        new_norm = LogNorm(vmin=new_vmin, vmax=new_vmax)
        for m in meshes: m.set_norm(new_norm)
        cbar.update_normal(meshes[-1])
        fig.canvas.draw_idle()

    vmin_slider.observe(update_clim, names='value')
    vmax_slider.observe(update_clim, names='value')
    display(widgets.HBox([vmin_slider, vmax_slider]))
    plt.show()


def temperature2D_grid(path='thermal/', vmin=1e0, vmax=1e3, cmap='gnuplot2',
                       xlim=None, ylim=None, figsize=(10, 14)):
    """
    Plots 2D poloidal dust temperature maps from a single model run.
    """
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[grid.columns[0]].values

    temp = pd.read_table(path + 'dust_temperature.dat', engine='python', header=None, skiprows=3)
    temp = temp[0].values
    nspecies = int(len(temp) / (nr * nt))
    temp = np.reshape(temp, (nspecies, nt, nr))
    grid = np.array(grid, copy=True)
    r_edge = grid[:nr+1] / autocm
    theta_edge = grid[nr+1:nr+1+nt+1]
    theta_edge[-1] = np.pi
    rr_edge, tt_edge = np.meshgrid(r_edge, theta_edge)
    R = rr_edge * np.sin(tt_edge)
    Z = rr_edge * np.cos(tt_edge)

    R_center = 0.25 * (R[:-1, :-1] + R[1:, :-1] + R[:-1, 1:] + R[1:, 1:])
    Z_center = 0.25 * (Z[:-1, :-1] + Z[1:, :-1] + Z[:-1, 1:] + Z[1:, 1:])

    sizes_file = path + 'dust_sizes.inp'
    if os.path.isfile(sizes_file):
        sizes = np.loadtxt(sizes_file)
        sizes = np.atleast_1d(sizes)
    else:
        sizes = None

    npanels = nspecies + 1
    ncols = min(nspecies, 4)
    nrows = int(np.ceil(npanels / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    for idx in range(nrows * ncols):
        ax = axes.flat[idx]
        if idx < nspecies:
            im = ax.pcolormesh(R, Z, temp[idx], cmap=cmap, shading='auto',
                               norm=LogNorm(vmin=vmin, vmax=vmax))
            ax.contour(R_center, Z_center, temp[idx], levels=[20], colors='black', linewidths=2)
            ax.set_title(f'bin {idx+1}', fontsize=12)
            if sizes is not None and idx < len(sizes):
                s = sizes[idx]
                size_label = f'{s/1e3:.1f} mm' if s >= 1e3 else f'{s:.2f} ' + r'$\mu$m'
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

        if xlim: ax.set_xlim(xlim)
        if ylim: ax.set_ylim(ylim)

    for ax in axes[-1, :]:
        if ax.get_visible(): ax.set_xlabel('r [au]', fontsize=14)
    for ax in axes[:, 0]: ax.set_ylabel('z [au]', fontsize=14)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label=r'T [K]')
    plt.show()


def midplane_temp(path='thermal/', xlim=None, ylim=None):
    """
    Plots midplane temperature profiles for all dust species.
    """
    grid = pd.read_table(path+'amr_grid.inp', engine='python', skiprows=5)
    head = grid.columns
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[head[0]].values
    try:
        temp = pd.read_table(path+'dust_temperature.dat', engine='python', header=None, skiprows=3)
    except IOError:
        print('plot.midplane_temp: dust_temperature.dat missing. Run thermal simulation first.')
        sys.exit(1)
    temp = temp[0].values
    nbspecies = int(len(temp)/(nr*nt))
    temp = np.reshape(temp, (nbspecies, nt, nr))
    grid = np.array(grid, copy=True)
    dist = grid[:nr+1]/autocm
    theta = grid[nr+1:nr+1+nt+1]
    theta[-1] = np.pi
    dist, tt = np.meshgrid(dist, theta)
    rr = dist*np.sin(tt)
    midtemp = temp[:, 90, :]
    radii = 0.5*(rr[90][0:rr[90].size-1] + rr[90][1:rr[90].size])

    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111)
    midtemp = pd.DataFrame(data=midtemp.transpose())
    for ispec in range(nbspecies):
        ax.plot(radii, midtemp[ispec].rolling(window=6, center=True).mean(), linewidth=2, label='bin: {}'.format(ispec+1))
        if xlim: ax.set_xlim(xlim)
        if ylim: ax.set_ylim(ylim)
    ax.set_xlabel(r'r [au]', fontsize=20)
    ax.set_ylabel(r'T [K]', fontsize=20)
    ax.legend(fontsize=15)
    ax.tick_params(labelsize=18)
    plt.show()


def vertical_temp(thermpath='thermal/', chempath='chemistry/', r=100):
    """
    Plots the vertical temperature profile at a specified disk radius coordinate.
    """
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
            print('plot.vertical_temp: radius {} does not exist.'.format(r))
            sys.exit(1)
        fig = plt.figure(figsize=(9.6, 8.2))
        ax = fig.add_subplot(111)
        ax.plot(temp[5], temp[0], linewidth=2, label='{} AU'.format(r))
        ax.set_ylabel(r'z [au]', fontsize=20)
        ax.set_xlabel(r'T$_\mathrm{d}$ [K]', fontsize=20)
        ax.legend(fontsize=15)
        ax.tick_params(labelsize=18)
        plt.show()
    elif nbspecies > 1:
        try:
            static = pd.read_table(chempath+str(r)+'AU/1D_static.dat', sep=r"\s+", engine='python', header=None, comment='!')
            temp = pd.read_table(chempath+str(r)+'AU/temperatures.dat', sep=r"\s+", engine='python', header=None)
        except IOError:
            print('plot.vertical_temp: radius = {} au does not exist.'.format(r))
            sys.exit(1)
        fig = plt.figure(figsize=(9.6, 8.2))
        ax = fig.add_subplot(111)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.91, 0.05, '{} AU'.format(r), transform=ax.transAxes, fontsize=16, bbox=props)
        for ai in range(nbspecies):
            ax.plot(temp[ai], static[0], linewidth=2, label='bin: {}'.format(ai+1))
        ax.set_ylabel(r'z [au]', fontsize=20)
        ax.set_xlabel(r'T$_\mathrm{d}$ [K]', fontsize=20)
        ax.legend(fontsize=15)
        ax.tick_params(labelsize=18)
        plt.show()


def avz(chempath='thermal/', r=100):
    """
    Plots the vertical visual extinction profile (Av) as a function of disk altitude.
    """
    static = pd.read_table(chempath+str(r)+'AU/1D_static.dat', sep=r"\s+", engine='python', header=None, comment='!', skiprows=1)
    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.91, 0.05, '{} AU'.format(r), transform=ax.transAxes, fontsize=16, bbox=props)
    ax.plot(static[3], static[0], linewidth=2, label='vertical Av')
    ax.set_xlabel(r'z [au]', fontsize=20)
    ax.set_ylabel(r'A$_\mathrm{\nu}$ [mag]', fontsize=20)
    ax.legend(fontsize=15)
    ax.tick_params(labelsize=18)
    plt.show()


def opacity(path='thermal/'):
    """
    Plots dust grain opacities across wavelengths for absorption and scattering components.
    """
    opaclist = sorted(glob.glob(path+'dustkap*'))

    # Absorption components
    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111) 
    ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize=18)
    ax.set_ylabel(r'$\kappa_\mathrm{abs}$ [cm$^2$/g]', fontsize=18)
    ax.set_xlim(1e-1, 1e4)
    ax.set_ylim(1e-2, 1e5)
    for opac in opaclist:
        name = opac.split("_")[1].split(".")[0]
        kappa = pd.read_table(opac, sep=r"\s+", comment='#', header=None, skiprows=10)
        ax.loglog(kappa[0], kappa[1], linewidth=2, label=name)
    ax.tick_params(labelsize=22)
    ax.legend(fontsize=15)
    plt.show()

    # Scattering components tracks
    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111) 
    ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize=18)
    ax.set_ylabel(r'$\kappa_\mathrm{scat}$ [cm$^2$/g]', fontsize=18)
    ax.set_xlim(1e-1, 1e4)
    ax.set_ylim(1e-2, 1e5)
    for opac in opaclist:
        name = opac.split("_")[1].split(".")[0]
        kappa = pd.read_table(opac, sep=r"\s+", comment='#', header=None, skiprows=10)
        ax.loglog(kappa[0], kappa[2], linewidth=2, label=name)
    ax.tick_params(labelsize=22)
    ax.legend(fontsize=15)
    plt.show()

    # Scattering angles cosine tracks
    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111) 
    ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize=18)
    ax.set_ylabel(r'<cos($\theta$)>', fontsize=18)
    ax.set_xlim(1e-1, 1e4)
    for opac in opaclist:
        name = opac.split("_")[1].split(".")[0]
        kappa = pd.read_table(opac, sep=r"\s+", comment='#', header=None, skiprows=10)
        ax.loglog(kappa[0], kappa[3], linewidth=2, label=name)
    ax.tick_params(labelsize=22)
    ax.legend(fontsize=15)
    plt.show()


def localflux(path='thermal/'):
    """
    Plots the midplane local mean radiation field intensity profile across wavelengths.
    """
    flux = pd.read_table(path+'mean_intensity.out', sep=r"\s+", comment='#', header=None, skiprows=4)[0].values
    grid = pd.read_table(path+'amr_grid.inp', engine='python', skiprows=5)
    
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
    nlam = int(len(flux)/(nr*nt))
    flux = np.reshape(flux, (nlam, nt, nr))
    midflux = flux[:, 90, :]

    fig = plt.figure(figsize=(9.6, 8.2))
    ax = fig.add_subplot(111)
    midflux_df = pd.DataFrame(data=midflux.transpose())
    for ilam in range(0, nlam, 2):
        ax.semilogy(radii, midflux_df[ilam].rolling(window=5, center=True).mean(), linewidth=1)
    ax.set_xlim(1, 200)
    ax.set_ylim(1e-30, 1e-10)
    ax.set_xlabel(r'r [au]', fontsize=20)
    ax.set_ylabel(r'Flux', fontsize=20)
    ax.grid()
    ax.tick_params(labelsize=22)
    plt.show()


def image(pathfile='thermal/', distance=100, vmin=1e-10, vmax=1e3, cmap='gnuplot2', labels=None):
    """
    Plots multi-wavelength RADMC-3D continuum synthetic spatial intensity maps (Stokes I).
    """
    with open(pathfile, 'r') as f:
        iformat = int(f.readline())                                 
        npix_x, npix_y = [int(x) for x in f.readline().split()]   
        nlam = int(f.readline())                                    
        pix_cm, _ = [float(x) for x in f.readline().split()]       
        wavelengths = [float(f.readline()) for _ in range(nlam)]   

    pix_au = pix_cm / autocm
    box_au = npix_x * pix_au
    half_box = box_au / 2.0

    distance_cm = distance * 3.086e18
    omega_pix = (pix_cm / distance_cm) ** 2
    to_jy = 1e23 * omega_pix

    data = np.loadtxt(pathfile, skiprows=4 + nlam + 1)
    nstokes = 4 if iformat == 3 else 1
    data = np.reshape(data, (nlam, npix_y, npix_x, nstokes))
    extent = [-half_box, half_box, -half_box, half_box]
    imgs = [data[i, :, :, 0] * to_jy for i in range(nlam)]

    if labels is None:
        def _wave_label(lam):
            if lam >= 1000: return f'Stokes I - {lam/1000:.2g} mm'
            elif lam >= 1: return f'Stokes I - {lam:.4g} μm'
            return f'Stokes I - {lam*1000:.4g} nm'
        labels = [_wave_label(w) for w in wavelengths]

    fig, axes = plt.subplots(1, nlam, figsize=(5*nlam, 5), sharex=True, sharey=True)
    if nlam == 1: axes = [axes]

    for i, ax in enumerate(axes):
        im = ax.imshow(imgs[i], origin='lower', extent=extent, cmap=cmap,
                       norm=LogNorm(vmin=vmin, vmax=vmax), interpolation='nearest')
        ax.tick_params(labelsize=17)
        ax.set_xlabel(r'x [au]', fontsize=17)
        if i == 0: ax.set_ylabel(r'y [au]', fontsize=17)
        ax.text(0.05, 0.95, labels[i], color='red', transform=ax.transAxes, fontsize=16, fontweight='bold', va='top')

    cbar = fig.colorbar(im, ax=axes, location='right', fraction=0.025, pad=0.02)
    cbar.set_label(r'$I_\nu$ [Jy pixel$^{-1}$]', fontsize=17)
    plt.show()


def image_vertical_cut(pathfile='thermal/', distance=100, xlim=None, ylim=None, labels=None, figsize=(9.6, 8.2)):
    """
    Plots 1D vertical cuts (along y at x=0) of Stokes I for each wavelength.
    """
    with open(pathfile, 'r') as f:
        iformat = int(f.readline())
        npix_x, npix_y = [int(x) for x in f.readline().split()]
        nlam = int(f.readline())
        pix_cm, _ = [float(x) for x in f.readline().split()]
        wavelengths = [float(f.readline()) for _ in range(nlam)]

    pix_au = pix_cm / autocm
    half_box = (npix_y * pix_au) / 2.0
    to_jy = 1e23 * ((pix_cm / (distance * 3.086e18)) ** 2)

    data = np.loadtxt(pathfile, skiprows=4 + nlam + 1)
    data = np.reshape(data, (nlam, npix_y, npix_x, 4 if iformat == 3 else 1))
    y_au = np.linspace(-half_box, half_box, npix_y)
    ix0 = npix_x // 2

    if labels is None:
        labels = [f'{w/1000:.2g} mm' if w >= 1000 else (f'{w:.4g} μm' if w >= 1 else f'{w*1000:.4g} nm') for w in wavelengths]

    fig, ax = plt.subplots(figsize=figsize)
    for i in range(nlam):
        ax.semilogy(y_au, data[i, :, ix0, 0] * to_jy, linewidth=2, label=labels[i])

    ax.set_xlabel(r'y [au]', fontsize=20)
    ax.set_ylabel(r'$I_\nu$ [Jy pixel$^{-1}$]', fontsize=20)
    ax.legend(fontsize=15)
    ax.tick_params(labelsize=18)
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    plt.show()


def numberdens(species='CO', path='thermal/', vmin=1e0, vmax=1e8, cmap='gnuplot2',
               ncols=3, xlim=None, ylim=None, figsize=None, save=False, savename='numberdens.pdf'):
    """
    Plots 2D maps of molecular number densities from files.
    """
    species_list = [species] if isinstance(species, str) else list(species)
    nspecies = len(species_list)

    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = np.array(grid[grid.columns[0]].values, copy=True)

    r_edge = grid[:nr+1] / autocm
    theta_edge = grid[nr+1:nr+1+nt+1]
    theta_edge[-1] = np.pi
    r_cen = 0.5 * (r_edge[:-1] + r_edge[1:])
    theta_cen = 0.5 * (theta_edge[:-1] + theta_edge[1:])
    rr, tt = np.meshgrid(r_cen, theta_cen)
    R = rr * np.sin(tt)
    Z = rr * np.cos(tt)

    ncols = min(ncols, nspecies)
    nrows = int(np.ceil(nspecies / ncols))
    if figsize is None: figsize = (5 * ncols + 1, 4 * nrows)

    norm = LogNorm(vmin=vmin, vmax=vmax)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()

    im = None
    for idx, sp in enumerate(species_list):
        ax = axes[idx]
        nd = pd.read_table(path + f'numberdens_{sp}.inp', engine='python', header=None, skiprows=2)
        nd = nd[0].values.reshape(nt, nr)
        nd = np.where(nd <= 0, 1e-100, nd)

        im = ax.pcolormesh(R, Z, nd, cmap=cmap, shading='gouraud', norm=norm, rasterized=True)
        ax.set_title(sp, fontsize=13)
        if xlim: ax.set_xlim(xlim)
        if ylim: ax.set_ylim(ylim)

    for idx in range(nspecies, len(axes)): axes[idx].set_visible(False)
    for ax in axes[(nrows - 1) * ncols:]: ax.set_xlabel('R [au]', fontsize=13)
    for i in range(nrows): axes[i * ncols].set_ylabel('Z [au]', fontsize=13)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label=r'$n$ [cm$^{-3}$]')
    if save: fig.savefig(savename, bbox_inches='tight')
    plt.show()


def plot_velocity_and_temperature(path='thermal/', vmin=0.0, vmax=10.0, logscale=False, cmap_v='viridis',
                                  Tmin=10.0, Tmax=1000.0, logscale_T=True, cmap_T='inferno',
                                  xlim=None, ylim=None, figsize=None, save=False, savename='gas_properties.pdf'):
    """
    Plots a 2D side-by-side comparison of gas velocity (v_phi) and gas temperature.
    """
    grid = pd.read_table(path + 'amr_grid.inp', engine='python', skiprows=5)
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = np.array(grid[grid.columns[0]].values, copy=True)

    r_edge = grid[:nr+1] / autocm
    theta_edge = grid[nr+1:nr+1+nt+1]
    theta_edge[-1] = np.pi
    r_cen = 0.5 * (r_edge[:-1] + r_edge[1:])
    theta_cen = 0.5 * (theta_edge[:-1] + theta_edge[1:])
    rr, tt = np.meshgrid(r_cen, theta_cen)
    R, Z = rr * np.sin(tt), rr * np.cos(tt)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize if figsize else (12, 5), sharex=True, sharey=True)
    norm_v = LogNorm(vmin=max(1e-2, vmin), vmax=vmax) if logscale else Normalize(vmin=vmin, vmax=vmax)

    vphi_all = pd.read_table(path + 'gas_velocity.inp', engine='python', skiprows=2, sep=r'\s+')[2].values / 1e5
    nphi = len(vphi_all) // (nt * nr)
    im1 = ax1.pcolormesh(R, Z, vphi_all.reshape(nphi, nt, nr)[0, :, :], cmap=cmap_v, shading='gouraud', norm=norm_v, rasterized=True)
    ax1.set_title('Gas Velocity' + (' (Log)' if logscale else ''), fontsize=13)
    ax1.set_xlabel('R [au]')
    ax1.set_ylabel('Z [au]')

    t_all = pd.read_table(path + 'gas_temperature.inp', engine='python', skiprows=2)[0].values
    t_2d = t_all.reshape(nphi, nt, nr)[0, :, :]
    norm_T = LogNorm(vmin=max(1e-1, Tmin), vmax=Tmax) if logscale_T else Normalize(vmin=Tmin, vmax=Tmax)
    
    im2 = ax2.pcolormesh(R, Z, t_2d, cmap=cmap_T, shading='gouraud', norm=norm_T, rasterized=True)
    ax2.set_title('Gas Temperature', fontsize=13)
    ax2.set_xlabel('R [au]')

    if xlim: ax1.set_xlim(xlim)
    if ylim: ax1.set_ylim(ylim)

    fig.colorbar(im1, cax=make_axes_locatable(ax1).append_axes("right", size="5%", pad=0.07), label=r'$v_\phi$ [km s$^{-1}$]')
    fig.colorbar(im2, cax=make_axes_locatable(ax2).append_axes("right", size="5%", pad=0.07), label=r'$T_{\rm gas}$ [K]')
    fig.tight_layout()
    if save: fig.savefig(savename, bbox_inches='tight')
    plt.show()


def static(chempath='chemistry/', column='nH', vmin=1, vmax=50, iso=None, cmap='gnuplot2',
           xlim=None, ylim=None, figsize=(6, 6), save=False, savename='filename.pdf'):
    """
    Plots a 2D map of a chosen physical parameter column from the 1D_static.dat files.
    """
    columns = ['z', 'nH', 'Tg', 'Av', 'diff', 'Td', 'inv_ab', 'conv_factor', 'a', 'uv']
    folders = [d for d in os.listdir(chempath) if os.path.isdir(os.path.join(chempath, d)) and re.match(r'^\d+AU$', d)]
    rchem = sorted([int(d.replace('AU', '')) for d in folders])

    all_data = []
    for r in rchem:
        df = pd.read_table(os.path.join(chempath, f'{r}AU', '1D_static.dat'), sep=r"\s+", comment='!', header=None, engine='python')
        df.columns = columns
        all_data.append(df)

    nbz_max = max(len(d) for d in all_data)
    static_map = np.full((nbz_max, len(rchem)), np.nan)
    temp_map = np.full((nbz_max, len(rchem)), np.nan)
    zz = np.zeros((nbz_max, len(rchem)))

    for idx, data in enumerate(all_data):
        nbz_r = len(data)
        start = nbz_max - nbz_r
        static_map[start:, idx] = data[column].values
        temp_map[start:, idx] = data['Td'].values
        z_col = data['z'].values
        zz[start:, idx] = z_col
        if start > 0:
            dz = (z_col[0] - z_col[1]) if nbz_r > 1 else z_col[0] * 0.1
            zz[:start, idx] = z_col[0] + np.arange(start, 0, -1) * dz

    rr, _ = np.meshgrid(rchem, np.arange(nbz_max))
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal', adjustable='box')
    t = ax.pcolormesh(rr, zz, static_map, cmap=cmap, shading='auto', norm=LogNorm(vmin=vmin, vmax=vmax), rasterized=True)
    fig.colorbar(t, pad=0.01, label=column)

    if iso is not None:
        ax.contour(rr, zz, temp_map, [iso] if not isinstance(iso, (list, tuple)) else iso, colors='black', linewidths=2.5)

    ax.set_xlabel(r'r [au]', fontsize=20)
    ax.set_ylabel(r'z [au]', fontsize=20)
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    if save: fig.savefig(savename, bbox_inches='tight')
    plt.show()


def nmgc_grainsizes(chempath='chemistry/', quantity='Td', vmin=None, vmax=None, cmap='gnuplot2',
                    xlim=None, ylim=None, figsize=(14, 10), save=False, savename='grain_sizes.pdf'):
    """
    Plots localized 2D maps broken down across each unique chemical model grain size bin.
    """
    static_columns = ['z', 'nH', 'Tg', 'Av', 'diff', 'Td', 'inv_ab', 'conv_factor', 'a', 'uv']
    folders = [d for d in os.listdir(chempath) if os.path.isdir(os.path.join(chempath, d)) and re.match(r'^\d+AU$', d)]
    rchem = sorted([int(d.replace('AU', '')) for d in folders])

    with open(os.path.join(chempath, f'{rchem[0]}AU', '1D_grain_sizes.in'), 'r') as f:
        for line in f:
            stripped = line.split('!')[0].strip()
            if not stripped: continue
            vals = stripped.split()
            grain_radii_um = np.array([float(v) for v in vals[:len(vals) // 4]]) * 1e4
            break
    ngrains = len(grain_radii_um)

    all_static = []
    for r in rchem:
        sd = pd.read_table(os.path.join(chempath, f'{r}AU', '1D_static.dat'), sep=r"\s+", comment='!', header=None, engine='python')
        sd.columns = static_columns
        all_static.append(sd)

    nbz_max = max(len(sd) for sd in all_static)
    data_map = np.full((ngrains, nbz_max, len(rchem)), np.nan)
    zz = np.zeros((nbz_max, len(rchem)))

    for idx, r in enumerate(rchem):
        static_data = all_static[idx]
        nbz_r = len(static_data)
        start = nbz_max - nbz_r
        zz[start:, idx] = static_data['z'].values
        if start > 0:
            dz = (zz[start, idx] - zz[start+1, idx]) if nbz_r > 1 else zz[start, idx] * 0.1
            zz[:start, idx] = zz[start, idx] + np.arange(start, 0, -1) * dz

        gs_lines = []
        with open(os.path.join(chempath, f'{r}AU', '1D_grain_sizes.in'), 'r') as f:
            for line in f:
                stripped = line.split('!')[0].strip()
                if not stripped: continue
                gs_lines.append([float(v) for v in stripped.split()])
        gs_array = np.array(gs_lines)

        inv_ab = gs_array[:, ngrains:2*ngrains]
        Td = gs_array[:, 2*ngrains:3*ngrains]

        for ig in range(ngrains):
            if quantity == 'Td':
                data_map[ig, start:, idx] = Td[:, ig]
            elif quantity == 'nd':
                data_map[ig, start:, idx] = static_data['nH'].values / inv_ab[:, ig]

    rr, _ = np.meshgrid(rchem, np.arange(nbz_max))
    ncols_plot = min(ngrains, 4)
    nrows_plot = int(np.ceil(ngrains / ncols_plot))
    fig, axes = plt.subplots(nrows_plot, ncols_plot, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)

    for ig in range(nrows_plot * ncols_plot):
        ax = axes.flat[ig]
        if ig < ngrains:
            ax.set_aspect('equal', adjustable='box')
            im = ax.pcolormesh(rr, zz, data_map[ig], cmap=cmap, shading='gouraud', norm=LogNorm(vmin=vmin, vmax=vmax), rasterized=True)
            s = grain_radii_um[ig]
            ax.text(0.05, 0.95, f'{s/1e3:.1f} mm' if s >= 1e3 else f'{s:.2f} µm', transform=ax.transAxes, fontsize=13, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5), va='top')
            if xlim: ax.set_xlim(xlim)
            if ylim: ax.set_ylim(ylim)
        else:
            ax.set_visible(False)

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    fig.colorbar(im, cax=fig.add_axes([0.90, 0.15, 0.02, 0.7]), label=r'T$_\mathrm{d}$ [K]' if quantity=='Td' else r'n$_\mathrm{d}$ [cm$^{-3}$]')
    if save: fig.savefig(savename, bbox_inches='tight')
    plt.show()


def plot_outputs_nautilus(PIPE, MODEL_NAMES, itime=-1, MODE='chemistry', key_list=['CO'],
                          fracab=True, verbose=True, xlim=None, ylim=None, colormap="plasma",
                          vmin=None, vmax=None, common_scale=True):
    """
    Plots a 2D vertical poloidal cross-section comparison across models.
    Each model maps strictly into its own dedicated column layout.
    """
    if len(MODEL_NAMES) != len(PIPE):
        raise ValueError("MODEL_NAMES and PIPE must have the same length")

    if isinstance(key_list, str):
        key_list = [key_list]

    if len(MODEL_NAMES) != len(set(MODEL_NAMES)):
        MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    def title_mol(mol_name, frac, path, verbose):
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
        return f"${f}$ [$n_{{{f}}}/n_H$]" if frac else f"${f}$ [$n_{{{f}}}$] [cm$^{{-3}}$]"

    def title_phys(variable):
        name = variable.replace("_", " ").title()
        if "temperature" in name.lower(): name += " [K]"
        elif "extinction" in name.lower(): name += " [mag]"
        elif "density" in name.lower(): name += " [$cm^{-3}$]"
        return name

    def get_grain_size_in_um(file_path, bin_index, verbose=False):
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('!'): continue
                    if '!' in line: line = line.split('!')[0].strip()
                    values = [float(val) for val in line.split()]
                    if not values: continue
                    num_grains = len(values) // 4
                    return values[:num_grains][int(bin_index) - 1] * 10000.0
            return None
        except FileNotFoundError:
            return None
    
    model_data = {}
    all_global_values = []
    
    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', f"{MODEL_NAMES[p_idx]}")
        main_output_dict = p.chemistry
        chempath = Path(p.chempath)
        species_data = {k: [] for k in key_list}
        
        for r_value in main_output_dict.keys():
            file_path = os.path.join(chempath, f"{r_value}AU", "1D_static.dat")
            if os.path.exists(file_path):
                try:
                    z_points = np.loadtxt(file_path, comments='!', usecols=0)
                    sub_dict = main_output_dict[r_value]
                    
                    for key in key_list:
                        if MODE == 'physical':
                            v_points = sub_dict[key][itime, :].copy()
                        elif MODE == 'chemistry':
                            abundance_array = sub_dict['abundances']
                            v_points = abundance_array.isel(time=itime).sel(species=key).values.copy()
                            if not fracab:
                                v_points *= sub_dict["H_number_density"][itime, :]
                        
                        if len(z_points) == len(v_points):
                            species_data[key].append({'R': float(r_value), 'z': np.array(z_points), 'v': np.array(v_points)})
                except Exception as e:
                    if verbose: print(f"Error processing {key} for R={r_value}: {e}")

        plot_structures = {}
        for key in key_list:
            columns_data = sorted(species_data[key], key=lambda x: x['R'])
            if not columns_data: continue
            polygons, values = [], []
            radii = [col['R'] for col in columns_data]
            
            if len(radii) > 1:
                r_midshifts = 0.5 * np.diff(radii)
                r_edges = [radii[0] - r_midshifts[0]] + [radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [radii[-1] + r_midshifts[-1]]
            else:
                r_edges = [radii[0] - 0.5, radii[0] + 0.5]
                
            for i, col in enumerate(columns_data):
                r_left, r_right = r_edges[i], r_edges[i+1]
                z_pts, v_pts = col['z'], col['v']
                z_midshifts = 0.5 * np.diff(z_pts)
                z_edges = [z_pts[0] - z_midshifts[0]] + [z_pts[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_pts[-1] + z_midshifts[-1])]
                
                for j in range(len(v_pts)):
                    poly = [(r_left, z_edges[j]), (r_right, z_edges[j]), (r_right, z_edges[j+1]), (r_left, z_edges[j+1])]
                    polygons.append(poly)
                    values.append(v_pts[j])
                    
            plot_structures[key] = {
                'polygons': polygons, 'values': np.array(values), 'radii': radii,
                'all_z': np.concatenate([c['z'] for c in columns_data])
            }
            all_global_values.extend(values)
        model_data[p_name] = plot_structures
    
    num_models = len(PIPE)
    model_names = list(model_data.keys())
    
    if num_models == 1:
        cols = min(3, len(key_list))
        rows = (len(key_list) + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4.5 * rows), squeeze=False)
        axes = axes.flatten()
    else:
        cols = num_models
        rows = len(key_list)
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4.5 * rows), squeeze=False)

    has_log_anywhere = MODE == 'chemistry' or any("density" in k.lower() or "extinction" in k.lower() for k in key_list)
    if common_scale and all_global_values:
        all_global_values = np.array(all_global_values)
        if has_log_anywhere:
            pos_vals = all_global_values[all_global_values > 0]
            global_vmin = vmin if vmin is not None else (max(1e-15, pos_vals.min()) if len(pos_vals) > 0 else 1e-15)
        else:
            global_vmin = vmin if vmin is not None else all_global_values.min()
        global_vmax = vmax if vmax is not None else all_global_values.max()
        
    for row_idx, key in enumerate(key_list):
        is_log = MODE == 'chemistry' or "density" in key.lower() or "extinction" in key.lower()

        for col_idx, p_name in enumerate(model_names):
            if key not in model_data[p_name]: continue
            ax = axes[row_idx, col_idx] if num_models > 1 else axes[row_idx]
            struct = model_data[p_name][key]
            vals = struct['values']
            
            if common_scale:
                actual_vmin, actual_vmax = global_vmin, global_vmax
            else:
                if is_log:
                    pos_vals = vals[vals > 0]
                    actual_vmin = vmin if vmin is not None else (max(1e-15, pos_vals.min()) if len(pos_vals) > 0 else 1e-15)
                else:
                    actual_vmin = vmin if vmin is not None else vals.min()
                actual_vmax = vmax if vmax is not None else vals.max()
            
            color_norm = LogNorm(vmin=actual_vmin, vmax=actual_vmax) if is_log else Normalize(vmin=actual_vmin, vmax=actual_vmax)
            coll = PolyCollection(struct['polygons'], array=vals, cmap=colormap, norm=color_norm, edgecolors='none')
            ax.add_collection(coll)
            
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


def plot_midplane_nautilus_multi_comparison(PIPE, MODEL_NAMES, itime=-1, MODE='chemistry', key_list=['CO'],
                                           fracab=True, verbose=True, xlim=None, ylim=None, colormap="turbo",
                                           vmin=None, vmax=None, figsize=None):
    """
    Plots and compares 1D midplane radial tracking profiles across models.
    Each curve utilizes its own independent spatial grids to ensure correct dimensions mapping.
    """
    if len(MODEL_NAMES) != len(PIPE):
        raise ValueError("MODEL_NAMES and PIPE must have the same length")
        
    if isinstance(key_list, str):
        key_list = [key_list]

    if len(MODEL_NAMES) != len(set(MODEL_NAMES)):
        MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    def clean_molec(mol_name):
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        return f"${f}$"

    def title_phys(variable):
        name = variable.replace("_", " ").title()
        if "temperature" in name.lower(): name += " [K]"
        elif "extinction" in name.lower(): name += " [mag]"
        elif "density" in name.lower(): name += " [$cm^{-3}$]"
        return name

    fig, ax = plt.subplots(figsize=figsize if figsize else (10, 6))
    colors = [plt.colormaps[colormap](0.5)] if len(key_list) == 1 else plt.colormaps[colormap](np.linspace(0, 0.9, len(key_list)))

    marker_pool = ['o', 's', '^', 'D', 'v', 'p', '*', 'X', 'h']
    linestyle_pool = ['-', '--', ':', '-.']

    all_valid_mins, all_valid_maxs = [], []
    time_years_string = None

    for k_idx, key in enumerate(key_list):
        for p_idx, p in enumerate(PIPE):
            p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
            main_output_dict = p.chemistry
            
            radii_list, values_list = [], []

            for r_value in main_output_dict.keys():
                try:
                    sub_dict = main_output_dict[r_value]
                    if MODE == 'physical':
                        v_midplane = sub_dict[key][itime, -1]
                    elif MODE == 'chemistry':
                        abundance_array = sub_dict['abundances']
                        v_midplane = float(abundance_array.isel(time=itime).sel(species=key).values[-1])
                        if not fracab:
                            v_midplane *= sub_dict["H_number_density"][itime, -1]
                    
                    radii_list.append(float(r_value))
                    values_list.append(v_midplane)
                    
                    if time_years_string is None:
                        try:
                            t_sec = sub_dict['abundances'].coords['time'].values[itime]
                            time_years_string = f"{t_sec / 3.156e7:.0f}"
                        except: pass
                except Exception as e:
                    if verbose: print(f"[{p_name}] Error processing R={r_value}, Key={key}: {e}")

            if not radii_list: continue

            sort_indices = np.argsort(radii_list)
            radii_arr = np.array(radii_list)[sort_indices]
            values_arr = np.array(values_list)[sort_indices]

            pos_values = values_arr[values_arr > 0]
            if len(pos_values) > 0: all_valid_mins.append(pos_values.min())
            all_valid_maxs.append(values_arr.max())

            clean_label = title_phys(key) if MODE == 'physical' else clean_molec(key)
            full_label = f"{p_name}" if len(key_list) == 1 else f"{clean_label} ({p_name})"

            ax.plot(radii_arr, values_arr, color=colors[k_idx], 
                    linestyle=linestyle_pool[p_idx % len(linestyle_pool)], 
                    marker=marker_pool[p_idx % len(marker_pool)], markersize=5, lw=1.5, label=full_label)

    is_log = MODE == 'chemistry' or any("density" in k.lower() or "extinction" in k.lower() for k in key_list)
    if is_log:
        ax.set_yscale('log')
        ax.set_ylim(vmin if vmin is not None else (max(1e-15, min(all_valid_mins)) if all_valid_mins else 1e-15),
                    vmax if vmax is not None else (max(all_valid_maxs)*2 if all_valid_maxs else 1.0))
    elif vmin is not None or vmax is not None:
        ax.set_ylim(vmin, vmax)

    ax.set_xlabel('Radius R [AU]')
    ax.set_ylabel(f"{clean_label}" if len(key_list) == 1 else "Values (See Legend)")
    ax.set_title(f'Midplane ($z = 0$) Profile Comparison — $t = {time_years_string}$ years' if time_years_string else 'Midplane ($z = 0$) Profile Comparison')
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.legend(loc='best', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.5)
    plt.tight_layout()
    plt.show()


def plot_midplane_nautilus_multi_CoverO(PIPE, MODEL_NAMES, itime=-1, MODE='chemistry',
                                        species_list=['CO', 'CO2', 'CH4', 'H2O'], fracab=True,
                                        verbose=True, xlim=None, ylim=None, colormap="turbo",
                                        vmin=None, vmax=None, figsize=None):
    """
    Plots individual molecular profiles alongside integrated local physical C/O ratio signatures.
    Each profile traces distinct coordinates to guarantee dimensions parity.
    """
    if len(MODEL_NAMES) != len(PIPE):
        raise ValueError("MODEL_NAMES and PIPE must have the same length")
        
    if isinstance(species_list, str):
        species_list = [species_list]

    if len(MODEL_NAMES) != len(set(MODEL_NAMES)):
        MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    def clean_molec(mol_name):
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        return f"${f}$"

    def count_atoms(mol_name, element):
        clean_name = re.sub(r"^[JK]\d+_*", "", mol_name).split('+')[0].split('-')[0]
        matches = re.findall(r'([A-Z][a-z]*)(\d*)', clean_name)
        count = 0
        for el, num in matches:
            if el == element: count += int(num) if num else 1
        return count

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize if figsize else (15, 6))
    colors = [plt.colormaps[colormap](0.5)] if len(species_list) == 1 else plt.colormaps[colormap](np.linspace(0, 0.9, len(species_list)))

    marker_pool = ['o', 's', '^', 'D', 'v', 'p', '*', 'X', 'h']
    linestyle_pool = ['-', '--', ':', '-.']

    all_valid_mins, all_valid_maxs = [], []
    time_years_string = None

    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        main_output_dict = p.chemistry
        
        marker_style = marker_pool[p_idx % len(marker_pool)]
        line_style = linestyle_pool[p_idx % len(linestyle_pool)]
        co_ratio_radii, co_ratio_values = [], []

        radii_sorted = sorted([float(r) for r in main_output_dict.keys()])
        
        for r_val in radii_sorted:
            r_str = str(r_val) if str(r_val) in main_output_dict else r_val
            if r_str not in main_output_dict: continue
            
            sub_dict = main_output_dict[r_str]
            total_C_density, total_O_density = 0.0, 0.0
            valid_cell = False

            for key in species_list:
                try:
                    if MODE == 'physical':
                        v_midplane = sub_dict[key][itime, -1]
                    elif MODE == 'chemistry':
                        abundance_array = sub_dict['abundances']
                        v_midplane = float(abundance_array.isel(time=itime).sel(species=key).values[-1])
                        if not fracab: v_midplane *= sub_dict["H_number_density"][itime, -1]
                    
                    valid_cell = True
                    if MODE == 'chemistry':
                        total_C_density += v_midplane * count_atoms(key, 'C')
                        total_O_density += v_midplane * count_atoms(key, 'O')
                except Exception as e:
                    if verbose and r_val == radii_sorted[0]: print(f"[{p_name}] Warning for Key={key}: {e}")
            
            if valid_cell and MODE == 'chemistry':
                co_ratio_radii.append(r_val)
                co_ratio_values.append(total_C_density / total_O_density if total_O_density > 0 else np.nan)

            if time_years_string is None:
                try:
                    t_sec = sub_dict['abundances'].coords['time'].values[itime]
                    time_years_string = f"{t_sec / 3.156e7:.0f}"
                except: pass

        if MODE == 'chemistry' and co_ratio_radii:
            ax2.plot(co_ratio_radii, co_ratio_values, color=plt.cm.tab10(p_idx % 10), 
                     linestyle=line_style, marker=marker_style, markersize=5, lw=1.5, label=f"C/O ({p_name})")

    for k_idx, key in enumerate(species_list):
        clean_label = clean_molec(key)
        for p_idx, p in enumerate(PIPE):
            p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
            radii_list, values_list = [], []
            
            for r_val in sorted([float(r) for r in p.chemistry.keys()]):
                try:
                    sub_dict = p.chemistry[str(r_val) if str(r_val) in p.chemistry else r_val]
                    if MODE == 'physical': v = sub_dict[key][itime, -1]
                    else:
                        v = float(sub_dict['abundances'].isel(time=itime).sel(species=key).values[-1])
                        if not fracab: v *= sub_dict["H_number_density"][itime, -1]
                    radii_list.append(r_val)
                    values_list.append(v)
                except: continue
            
            if radii_list:
                pos_values = np.array(values_list)[np.array(values_list) > 0]
                if len(pos_values) > 0: all_valid_mins.append(pos_values.min())
                all_valid_maxs.append(max(values_list))
                
                full_label = f"{p_name}" if len(species_list) == 1 else f"{clean_label} ({p_name})"
                ax1.plot(radii_list, values_list, color=colors[k_idx], 
                         linestyle=linestyle_pool[p_idx % len(linestyle_pool)], 
                         marker=marker_pool[p_idx % len(marker_pool)], markersize=5, lw=1.5, label=full_label)

    is_log = MODE == 'chemistry' or any("density" in k.lower() or "extinction" in k.lower() for k in species_list)
    if is_log:
        ax1.set_yscale('log')
        ax1.set_ylim(vmin if vmin is not None else (max(1e-15, min(all_valid_mins)) if all_valid_mins else 1e-15), 
                     vmax if vmax is not None else max(all_valid_maxs) * 2)
    else:
        if vmin is not None or vmax is not None: ax1.set_ylim(vmin, vmax)

    ax1.set_ylabel("Abundance" if fracab else "Density [cm$^{-3}$]")
    ax1.set_xlabel('Radius R [AU]')
    ax1.legend(loc='best', frameon=True, fontsize='small')
    ax1.grid(True, linestyle=':', alpha=0.5)

    if MODE == 'chemistry':
        ax2.set_ylabel("Gas+Ice C/O Ratio")
        ax2.legend(loc='best', frameon=True, fontsize='small')
        ax2.grid(True, linestyle=':', alpha=0.5)
    ax2.set_xlabel('Radius R [AU]')

    if xlim:
        ax1.set_xlim(xlim)
        ax2.set_xlim(xlim)
    if ylim: ax2.set_ylim(ylim)

    title_suffix = f' — $t = {time_years_string}$ years' if time_years_string else ''
    fig.suptitle(f'Midplane ($z = 0$) Radial Profile Comparison{title_suffix}', y=0.98, fontsize=12)
    plt.tight_layout()
    plt.show()


def plot_grain_surface_midplane_comparison(PIPE, MODEL_NAMES, itime=-1, verbose=True, xlim=None, ylim=None, colormap="viridis"):
    """
    Computes total integrated midplane grain surface area across separate pipelines, using independent spatial coordinate tracks.
    """
    if len(MODEL_NAMES) != len(PIPE): raise ValueError("MODEL_NAMES and PIPE must have the same length")
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)): MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.set_yscale('log')
    colors = plt.colormaps[colormap].resampled(len(PIPE))(np.linspace(0, 0.95, len(PIPE)))
    
    global_max_val, global_min_pos_val = -np.inf, np.inf
    time_years_string = None

    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        main_output_dict = p.chemistry
        chempath = Path(p.chempath)
        radii_list, surface_list = [], []

        for r_value in main_output_dict.keys():
            file_path = os.path.join(chempath, f"{r_value}AU", "1D_grain_sizes.in")
            if os.path.exists(file_path):
                try:
                    grain_data = np.loadtxt(file_path, comments='!')
                    midplane_row = grain_data[-1, :]
                    N = int(len(midplane_row) / 4)
                    
                    a_array = midplane_row[0:N]
                    gtodn_array = midplane_row[N:2*N]
                    nH_midplane = main_output_dict[r_value]["H_number_density"][itime, -1]
                    
                    total_surface = 4 * np.pi * nH_midplane * np.sum((a_array**2) / gtodn_array)
                    radii_list.append(float(r_value))
                    surface_list.append(total_surface)
                    
                    if time_years_string is None:
                        try:
                            t_sec = main_output_dict[r_value]['abundances'].coords['time'].values[itime]
                            time_years_string = f"{t_sec / 3.156e7:.0f}"
                        except: pass
                except Exception as e:
                    if verbose: print(f"[{p_name}] Error for R={r_value}: {e}")

        if radii_list:
            sort_indices = np.argsort(radii_list)
            radii_arr = np.array(radii_list)[sort_indices]
            surface_arr = np.array(surface_list)[sort_indices]
            
            global_max_val = max(global_max_val, surface_arr.max())
            pos_vals = surface_arr[surface_arr > 0]
            if len(pos_vals) > 0: global_min_pos_val = min(global_min_pos_val, pos_vals.min())

            ax.plot(radii_arr, surface_arr, color=colors[p_idx], linestyle='-', marker='s', markersize=4, lw=1.5, label=p_name)

    ax.set_xlabel('Radius R [AU]')
    ax.set_ylabel(r'Total Grain Surface Area [$\text{cm}^{2}/\text{cm}^{3}$]')
    ax.set_title(f'Total Grain Surface Area at Midplane ($z=0$) - $t = {time_years_string}$ yr' if time_years_string else 'Total Grain Surface Area at Midplane ($z=0$)')
    
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    else: ax.set_ylim(max(1e-25, global_min_pos_val), global_max_val * 2 if global_max_val != -np.inf else 1e-10)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.legend(loc='best', frameon=True)
    plt.tight_layout()
    plt.show()


def plot_vertical_cut_nautilus_comparison(PIPE, MODEL_NAMES, R, species='CO', itime=-1, fracab=True,
                                         colormap="turbo", xlim=None, ylim=None, xscale="linear", yscale="linear", verbose=True):
    """
    Plots and compares 1D vertical cuts across multiple tracking nodes without grid mismatch collisions.
    """
    if len(MODEL_NAMES) != len(PIPE): raise ValueError("MODEL_NAMES and PIPE must have the same length")
    r_list = [R] if not isinstance(R, list) else R
    species_list = [species] if not isinstance(species, list) else species
    r_list = [int(r) for r in r_list]

    if len(r_list) > 1 and len(species_list) > 1:
        raise ValueError("Cannot supply multiple values for both 'R' and 'species' simultaneously.")

    if len(MODEL_NAMES) != len(set(MODEL_NAMES)): MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    fig, ax = plt.subplots(figsize=(8, 6))
    total_distinct_targets = max(len(r_list), len(species_list))
    colors = plt.colormaps[colormap](np.linspace(0, 0.9, total_distinct_targets))

    marker_pool = ['o', 's', '^', 'D', 'v', 'p', '*', 'X']
    linestyle_pool = ['-', '--', ':', '-.']
    time_years_string = None

    target_idx = 0
    for r_val in r_list:
        for spec_val in species_list:
            for p_idx, p in enumerate(PIPE):
                p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
                chempath = Path(p.chempath)
                if r_val not in p.chemistry: continue

                try:
                    static_file = chempath / f"{r_val}AU" / "1D_static.dat"
                    static = pd.read_table(static_file, sep=r'\s+', comment='!', header=None, engine='python')
                    z = static[0].values  
                    
                    ab = p.chemistry[r_val]['abundances']
                    sp_arr = ab.isel(time=itime).sel(species=spec_val).values
                    n_plot = sp_arr if fracab else sp_arr * p.chemistry[r_val]["H_number_density"][itime, :]

                    if time_years_string is None:
                        try: time_years_string = f"{ab.coords['time'].values[itime] / 3.156e7:.0f}"
                        except: pass

                    marker_style = marker_pool[p_idx % len(marker_pool)]
                    line_style = linestyle_pool[p_idx % len(linestyle_pool)]
                    
                    ax.scatter(n_plot, z, color=colors[target_idx], s=15, marker=marker_style)
                    ax.plot(n_plot, z, color=colors[target_idx], linestyle=line_style, lw=1.5,
                            label=f"{spec_val} @ {r_val} AU ({p_name})" if total_distinct_targets > 1 else f"{p_name}")
                except Exception as e:
                    if verbose: print(f"[{p_name}] Error parsing vertical cut for R={r_val}, species={spec_val}: {e}")
            if total_distinct_targets > 1: target_idx += 1

    ax.set_ylabel("z [AU]")
    ax.set_xlabel("Fractional Abundance" if fracab else "Number Density [cm$^{-3}$]")
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.legend(loc='best', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.5)
    plt.tight_layout()
    plt.show()


def plot_atom_ratio_nautilus_midplane(PIPE, MODEL_NAMES, ratio_list=['C/O'], itime=-1,
                                      verbose=True, xlim=None, ylim=None, common_scale=True):
    """
    Plots midplane elemental ratios. Each model maps cleanly to its own independent column tracking framework.
    """
    if len(MODEL_NAMES) != len(PIPE): raise ValueError("MODEL_NAMES and PIPE must have the same length")
    if isinstance(ratio_list, str): ratio_list = [ratio_list]

    elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    parsed_ratios = []
    for item in ratio_list:
        s1, s2 = item.split('/')
        parsed_ratios.append((s1, s2, item))

    if len(MODEL_NAMES) != len(set(MODEL_NAMES)): MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    def count_species_elements(species_name, element1, element2):
        if species_name == 'e-': return {element1: 0, element2: 0}
        formula = species_name.replace('c-', '').replace('l-', '').split('+')[0].split('-')[0]
        pattern = re.compile(r'([A-Z][a-z]?)(-?\d*)')
        composition = {}
        for atom, n in pattern.findall(formula):
            composition[atom] = composition.get(atom, 0) + (int(n) if n else 1)
        return {element1: composition.get(element1, 0), element2: composition.get(element2, 0)}

    def keep_gas_species_only(species):
        motif = re.compile(r'^[JK]\d{2}|^GRAIN')
        return [e for e in species if not motif.match(e)]

    model_data = {}
    atom_cache = {}
    time_years_string = None

    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        ratio_database = {token: [] for _, _, token in parsed_ratios}

        for r_value in p.chemistry.keys():
            file_path = os.path.join(p.chempath, f"{r_value}AU", "1D_static.dat")
            if os.path.exists(file_path):
                try:
                    z_points = np.atleast_1d(np.loadtxt(file_path, comments='!', usecols=0))
                    abundance_array = p.chemistry[r_value]['abundances']
                    
                    if time_years_string is None:
                        try: time_years_string = f"{float(abundance_array.coords['time'].values[itime]) / 3.156e7:.0f}"
                        except: pass

                    local_species_list = keep_gas_species_only(list(abundance_array.coords['species'].values))
                    time_slice = abundance_array.isel(time=itime).sel(species=local_species_list)
                    
                    for s1, s2, token in parsed_ratios:
                        s1_coeffs, s2_coeffs = [], []
                        for species in local_species_list:
                            cache_key = f"{species}_{s1}_{s2}"
                            if cache_key not in atom_cache:
                                atom_cache[cache_key] = count_species_elements(species, s1, s2)
                            s1_coeffs.append(atom_cache[cache_key][s1])
                            s2_coeffs.append(atom_cache[cache_key][s2])
                        
                        da_s1 = xr.DataArray(s1_coeffs, coords=[local_species_list], dims=['species'])
                        da_s2 = xr.DataArray(s2_coeffs, coords=[local_species_list], dims=['species'])
                        
                        total_s1 = (time_slice * da_s1).sum(dim='species')
                        total_s2 = (time_slice * da_s2).sum(dim='species')

                        with np.errstate(divide='ignore', invalid='ignore'):
                            ratio_values = xr.where(total_s2 > 0, total_s1 / total_s2, 0.0).values
                        
                        if len(z_points) == len(ratio_values):
                            ratio_database[token].append({'R': float(r_value), 'v_midplane': ratio_values[-1]})
                except Exception as e:
                    if verbose: print(f"[{p_name}] Error for R={r_value} AU: {e}")

        plot_structures = {}
        for s1, s2, token in parsed_ratios:
            columns_data = sorted(ratio_database[token], key=lambda x: x['R'])
            if not columns_data: continue
            plot_structures[token] = {'radii': np.array([c['R'] for c in columns_data]), 'values': np.array([c['v_midplane'] for c in columns_data])}
        model_data[p_name] = plot_structures

    num_models = len(PIPE)
    model_names = list(model_data.keys())
    num_ratios = len(parsed_ratios)

    if num_models == 1:
        cols = min(3, num_ratios)
        rows = (num_ratios + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), squeeze=False)
        axes = axes.flatten()
    else:
        cols = num_models
        rows = num_ratios
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), squeeze=False)

    for row_idx, (_, _, token) in enumerate(parsed_ratios):
        global_is_log, global_ymin, global_ymax = False, 0.0, 1.0
        if common_scale:
            all_row_vals = []
            for p_name in model_names:
                if token in model_data[p_name]: all_row_vals.extend(model_data[p_name][token]['values'])
            all_row_vals = np.array(all_row_vals)
            if len(all_row_vals) > 0:
                pos_vals = all_row_vals[all_row_vals > 0]
                global_ymin = pos_vals.min() if len(pos_vals) > 0 else 1e-15
                global_ymax = all_row_vals.max()
                if len(pos_vals) > 0 and (global_ymax / global_ymin) > 100.0: global_is_log = True

        for col_idx, p_name in enumerate(model_names):
            if token not in model_data[p_name]: continue
            ax = axes[row_idx, col_idx] if num_models > 1 else axes[row_idx]
            struct = model_data[p_name][token]
            vals = struct['values']
            
            is_log = global_is_log if common_scale else (len(vals[vals>0]) > 0 and (vals.max()/vals[vals>0].min()) > 100.0)
            if is_log: ax.set_yscale('log')
            
            ax.plot(struct['radii'], vals, color='black', lw=2, marker='o', ms=4)
            ax.set_title(f"{p_name}\nRatio {token} (z=0)")
            ax.set_xlabel('R [AU]')
            ax.set_ylabel(f'Ratio {token}')
            ax.grid(True, which="both", linestyle='--', alpha=0.5)
            if xlim: ax.set_xlim(xlim)
            if ylim: ax.set_ylim(ylim)
            else: ax.set_ylim(global_ymin if common_scale else vals[vals>0].min(), global_ymax if common_scale else vals.max())

    plt.tight_layout()
    plt.show()


def plot_atom_ratio_nautilus_comparison(PIPE, MODEL_NAMES, ratio_list=['C/O'], itime=-1,
                                       verbose=True, xlim=None, ylim=None, colormap="gnuplot", common_scale=True):
    """
    Plots 2D vertical atomic cuts mapping models column-by-column without processing differences or residuals.
    """
    if len(MODEL_NAMES) != len(PIPE): raise ValueError("MODEL_NAMES and PIPE must have the same length")
    if isinstance(ratio_list, str): ratio_list = [ratio_list]

    elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    parsed_ratios = []
    for item in ratio_list:
        s1, s2 = item.split('/')
        parsed_ratios.append((s1, s2, item))

    if len(MODEL_NAMES) != len(set(MODEL_NAMES)): MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    def count_species_elements(species_name, element1, element2):
        if species_name == 'e-': return {element1: 0, element2: 0}
        formula = species_name.replace('c-', '').replace('l-', '').split('+')[0].split('-')[0]
        pattern = re.compile(r'([A-Z][a-z]?)(-?\d*)')
        composition = {}
        for atom, n in pattern.findall(formula):
            composition[atom] = composition.get(atom, 0) + (int(n) if n else 1)
        return {element1: composition.get(element1, 0), element2: composition.get(element2, 0)}

    def keep_gas_species_only(species):
        motif = re.compile(r'^[JK]\d{2}|^GRAIN')
        return [e for e in species if not motif.match(e)]

    model_data = {}
    atom_cache = {}
    time_years_string = None

    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        ratio_database = {token: [] for _, _, token in parsed_ratios}

        for r_value in p.chemistry.keys():
            file_path = os.path.join(p.chempath, f"{r_value}AU", "1D_static.dat")
            if os.path.exists(file_path):
                try:
                    z_points = np.loadtxt(file_path, comments='!', usecols=0)
                    abundance_array = p.chemistry[r_value]['abundances']
                    local_species_list = keep_gas_species_only(list(abundance_array.coords['species'].values))
                    sliced_abundances = abundance_array.isel(time=itime).sel(species=local_species_list).values
                    
                    if time_years_string is None:
                        try: time_years_string = f"{abundance_array.coords['time'].values[itime] / 3.156e7:.0f}"
                        except: pass

                    for s1, s2, token in parsed_ratios:
                        s1_coeffs = [count_species_elements(sp, s1, s2)[s1] for sp in local_species_list]
                        s2_coeffs = [count_species_elements(sp, s1, s2)[s2] for sp in local_species_list]
                        
                        s1_coeffs = np.array(s1_coeffs)[:, np.newaxis]
                        s2_coeffs = np.array(s2_coeffs)[:, np.newaxis]
                        
                        total_s1 = np.sum(sliced_abundances * s1_coeffs, axis=0)
                        total_s2 = np.sum(sliced_abundances * s2_coeffs, axis=0)

                        with np.errstate(divide='ignore', invalid='ignore'):
                            v_points = np.where(total_s2 > 0, total_s1 / total_s2, 0.0)
                        
                        if len(z_points) == len(v_points):
                            ratio_database[token].append({'R': float(r_value), 'z': np.array(z_points), 'v': np.array(v_points)})
                except Exception as e:
                    if verbose: print(f"[{p_name}] Error compiling ratios for R={r_value}: {e}")

        plot_structures = {}
        for s1, s2, token in parsed_ratios:
            columns_data = sorted(ratio_database[token], key=lambda x: x['R'])
            if not columns_data: continue
            radii = [col['R'] for col in columns_data]
            
            r_edges = [radii[0]-0.5, radii[0]+0.5] if len(radii) == 1 else [radii[0]-0.5*np.diff(radii)[0]] + [radii[i]+0.5*np.diff(radii)[i] for i in range(len(np.diff(radii)))] + [radii[-1]+0.5*np.diff(radii)[-1]]
            
            polygons, values = [], []
            for i, col in enumerate(columns_data):
                r_left, r_right = r_edges[i], r_edges[i+1]
                z_pts, v_pts = col['z'], col['v']
                z_edges = [z_pts[0]-0.5*np.diff(z_pts)[0]] + [z_pts[j]+0.5*np.diff(z_pts)[j] for j in range(len(np.diff(z_pts)))] + [max(0.0, z_pts[-1]+0.5*np.diff(z_pts)[-1])]
                for j in range(len(v_pts)):
                    poly = [(r_left, z_edges[j]), (r_right, z_edges[j]), (r_right, z_edges[j+1]), (r_left, z_edges[j+1])]
                    polygons.append(poly)
                    values.append(v_pts[j])

            plot_structures[token] = {
                'polygons': polygons, 'values': np.array(values), 'radii': radii,
                'all_z': np.concatenate([c['z'] for c in columns_data])
            }
        model_data[p_name] = plot_structures

    num_models = len(PIPE)
    model_names = list(model_data.keys())
    num_ratios = len(parsed_ratios)

    cols = num_models if num_models > 1 else min(3, num_ratios)
    rows = num_ratios if num_models > 1 else (num_ratios + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.5 * rows), squeeze=False)

    for row_idx, (_, _, token) in enumerate(parsed_ratios):
        if common_scale:
            all_row_vals = np.concatenate([model_data[m][token]['values'] for m in model_names if token in model_data[m]])
            global_row_vmin = all_row_vals[all_row_vals>0].min() if len(all_row_vals[all_row_vals>0]) > 0 else 1e-15
            global_row_vmax = all_row_vals.max()
            global_row_is_log = (global_row_vmax / global_row_vmin) > 10.0

        for col_idx, p_name in enumerate(model_names):
            if token not in model_data[p_name]: continue
            ax = axes[row_idx, col_idx]
            struct = model_data[p_name][token]
            vals = struct['values']
            
            vmin_loc = global_row_vmin if common_scale else (vals[vals>0].min() if len(vals[vals>0])>0 else 1e-15)
            vmax_loc = global_row_vmax if common_scale else vals.max()
            is_log = global_row_is_log if common_scale else (vmax_loc / vmin_loc > 10.0)
            
            color_norm = LogNorm(vmin=vmin_loc, vmax=vmax_loc) if is_log else Normalize(vmin=vmin_loc, vmax=vmax_loc)
            coll = PolyCollection(struct['polygons'], array=vals, cmap=colormap, norm=color_norm, edgecolors='none')
            ax.add_collection(coll)
            
            fig.colorbar(plt.cm.ScalarMappable(norm=color_norm, cmap=colormap), ax=ax, label=f"Ratio [{token}]")
            ax.set_title(f"{p_name}\nAtomic Ratio: {token}")
            ax.set_xlabel('R [AU]')
            ax.set_ylabel('z [AU]')
            ax.set_xlim(xlim if xlim else (0, max(struct['radii']) * 1.02))
            ax.set_ylim(ylim if ylim else (0, max(struct['all_z']) * 1.07))

    plt.tight_layout()
    plt.show()


def plot_top_contributing_species(chempath, main_output_dict, target_atom="C", itime=-1,
                                  verbose=True, spnumber=5, color="darkred", phase="gas", grain_bin=None):
    """
    Plots the top volume-integrated chemical carrier species reservoir contributions for a target element.
    """
    allowed_elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    if target_atom not in allowed_elements: raise ValueError(f"Atom '{target_atom}' unrecognized.")
    if grain_bin is not None and phase in ["gas", "all"]: raise ValueError("grain_bin and gas phase are mutually exclusive.")

    chempath = Path(chempath)

    def parse_species(species_name):
        if species_name == 'e-': return "gas", None, "e-"
        m = re.match(r'^([JK])(\d+)(.+)', species_name)
        if m: return ("surface" if m.group(1) == 'J' else "mantle"), m.group(2), m.group(3).replace('c-', '').replace('l-', '')
        return "gas", None, species_name.replace('c-', '').replace('l-', '')

    def count_target_atom(clean_formula, target):
        if clean_formula == 'e-': return 0
        f = clean_formula.split('+')[0].split('-')[0]
        return sum(int(n) if n else 1 for atom, n in re.findall(r'([A-Z][a-z]*)(\d*)', f) if atom == target)

    global_species_contributions = {}
    radii_map = {int(re.findall(r'\d+', str(k))[0]): k for k in main_output_dict.keys() if re.findall(r'\d+', str(k))}
    radii = sorted(list(radii_map.keys()))

    r_edges = [radii[0]*0.9, radii[0]*1.1] if len(radii)==1 else [radii[0]-0.5*np.diff(radii)[0]] + [radii[i]+0.5*np.diff(radii)[i] for i in range(len(np.diff(radii)))] + [radii[-1]+0.5*np.diff(radii)[-1]]

    for i, r_val in enumerate(radii):
        file_path = chempath / f"{r_val}AU" / "1D_static.dat"
        if os.path.exists(file_path):
            try:
                z_points = np.loadtxt(file_path, comments='!', usecols=0)
                sub_dict = main_output_dict[radii_map[r_val]]
                ab_arr = sub_dict['abundances']
                nH = sub_dict["H_number_density"][itime, :]

                dz = np.abs(np.diff([z_points[0]-0.5*np.diff(z_points)[0]] + [z_points[j]+0.5*np.diff(z_points)[j] for j in range(len(np.diff(z_points)))] + [max(0.0, z_points[-1]+0.5*np.diff(z_points)[-1])])) if len(z_points)>1 else np.array([z_points[0]])
                vol = 2 * np.pi * (r_val * 1.496e13) * ((r_edges[i+1]-r_edges[i])*1.496e13) * (dz * 1.496e13)

                for sp in ab_arr.coords['species'].values:
                    if "GRAIN" in sp: continue
                    p, b, f = parse_species(sp)
                    if phase != "all" and (phase == "grain" wholesaler and p not in ["surface", "mantle"] or phase != "grain" and p != phase): continue
                    if grain_bin is not None and b != str(grain_bin): continue

                    coef = count_target_atom(f, target_atom)
                    if coef > 0:
                        tot_atoms = np.sum(ab_arr.isel(time=itime).sel(species=sp).values * nH * vol) * coef
                        if tot_atoms > 0: global_species_contributions[f] = global_species_contributions.get(f, 0.0) + tot_atoms
            except Exception as e:
                if verbose: print(f"Error processing data for R={r_val}: {e}")

    tot_sum = sum(global_species_contributions.values())
    if tot_sum == 0: return

    sorted_sp = sorted({sp: (v/tot_sum)*100 for sp, v in global_species_contributions.items()}.items(), key=lambda x: x[1], reverse=True)[:spnumber]
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar([f"${re.sub(r'(\d+)', r'_{\1}', s[0])}$" for s in sorted_sp], [s[1] for s in sorted_sp], color=color, edgecolor='grey', alpha=0.85)
    for bar in bars: ax.annotate(f'{bar.get_height():.2f}%', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()), xytext=(0, 4), textcoords="offset points", ha='center', va='bottom', fontweight='bold')
    ax.set_ylabel('Global Contribution (%)')
    plt.show()


def plot_top_species_per_radius(chempath, main_output_dict, target_atom="C", itime=-1,
                                verbose=True, spnumber=5, phase="gas", grain_bin=None, cmap_name="tab10", rmin=None, rmax=None):
    """
    Plots the horizontal breakdown of the top N contributing carriers across radial shells.
    """
    allowed_elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    if target_atom not in allowed_elements: raise ValueError(f"Atom '{target_atom}' unrecognized.")

    chempath = Path(chempath)

    def parse_species(species_name):
        if species_name == 'e-': return "gas", None, "e-"
        m = re.match(r'^([JK])(\d+)(.+)', species_name)
        if m: return ("surface" if m.group(1) == 'J' else "mantle"), m.group(2), m.group(3).replace('c-', '').replace('l-', '')
        return "gas", None, species_name.replace('c-', '').replace('l-', '')

    def count_target_atom(clean_formula, target):
        if clean_formula == 'e-': return 0
        f = clean_formula.split('+')[0].split('-')[0]
        return sum(int(n) if n else 1 for atom, n in re.findall(r'([A-Z][a-z]*)(\d*)', f) if atom == target)

    radii_map = {int(re.findall(r'\d+', str(k))[0]): k for k in main_output_dict.keys() if re.findall(r'\d+', str(k))}
    radii = sorted([r for r in radii_map.keys() if (rmin is None or r >= rmin) and (rmax is None or r <= rmax)])

    r_edges = [radii[0]*0.9, radii[0]*1.1] if len(radii)==1 else [radii[0]-0.5*np.diff(radii)[0]] + [radii[i]+0.5*np.diff(radii)[i] for i in range(len(np.diff(radii)))] + [radii[-1]+0.5*np.diff(radii)[-1]]
    
    radial_plot_data = {}
    all_encountered = set()

    for i, r_val in enumerate(radii):
        file_path = chempath / f"{r_val}AU" / "1D_static.dat"
        if os.path.exists(file_path):
            try:
                z_points = np.loadtxt(file_path, comments='!', usecols=0)
                sub_dict = main_output_dict[radii_map[r_val]]
                ab_arr = sub_dict['abundances']
                nH = sub_dict["H_number_density"][itime, :]

                dz = np.abs(np.diff([z_points[0]-0.5*np.diff(z_points)[0]] + [z_points[j]+0.5*np.diff(z_points)[j] for j in range(len(np.diff(z_points)))] + [max(0.0, z_points[-1]+0.5*np.diff(z_points)[-1])])) if len(z_points)>1 else np.array([z_points[0]])
                vol = 2 * np.pi * (r_val * 1.496e13) * ((r_edges[i+1]-r_edges[i])*1.496e13) * (dz * 1.496e13)

                local_cont = {}
                for sp in ab_arr.coords['species'].values:
                    if "GRAIN" in sp: continue
                    p, b, f = parse_species(sp)
                    if phase != "all" and (phase == "grain" and p not in ["surface", "mantle"] or phase != "grain" and p != phase): continue
                    if grain_bin is not None and b != str(grain_bin): continue

                    coef = count_target_atom(f, target_atom)
                    if coef > 0:
                        tot_atoms = np.sum(ab_arr.isel(time=itime).sel(species=sp).values * nH * vol) * coef
                        if tot_atoms > 0: local_cont[f] = local_cont.get(f, 0.0) + tot_atoms

                r_sum = sum(local_cont.values())
                if r_sum > 0:
                    top_entries = sorted({sp: (v/r_sum)*100 for sp, v in local_cont.items()}.items(), key=lambda x: x[1], reverse=True)[:spnumber]
                    radial_plot_data[r_val] = top_entries
                    for sp, _ in top_entries: all_encountered.add(sp)
            except Exception as e:
                if verbose: print(f"Error processing R={r_val}: {e}")

    if not radial_plot_data: return
    unique_sp = sorted(list(all_encountered))
    cmap = plt.get_cmap(cmap_name)
    species_colors = {sp: cmap(idx/max(1, len(unique_sp)-1)) for idx, sp in enumerate(unique_sp)}

    y_labels, x_pct, bar_colors = [], [], []
    for r_val in sorted(list(radial_plot_data.keys()), reverse=True):
        for sp, pct in reversed(radial_plot_data[r_val]):
            y_labels.append(f"{r_val} AU — ${re.sub(r'(\d+)', r'_{\1}', sp)}$")
            x_pct.append(pct)
            bar_colors.append(species_colors[sp])

    fig, ax = plt.subplots(figsize=(11, max(4, len(y_labels)*0.45)))
    ax.barh(y_labels, x_pct, color=bar_colors, edgecolor='grey', alpha=0.85, height=0.7)
    ax.set_xlabel('Local Budget Contribution (%)')
    plt.tight_layout()
    plt.show()


def plot_top_species_per_radius_midplane(chempath, main_output_dict, target_atom="C", itime=-1,
                                         verbose=True, spnumber=5, phase="gas", grain_bin=None, cmap_name="tab10", rmin=None, rmax=None):
    """
    Computes and plots the horizontal tracking carrier breakdown limited strictly to the midplane (z=0) slice row elements.
    """
    allowed_elements = ['H', 'He', 'C', 'N', 'O', 'Si', 'S', 'Fe', 'Na', 'Mg', 'Cl', 'P', 'F']
    if target_atom not in allowed_elements: raise ValueError(f"Atom '{target_atom}' unrecognized.")

    def parse_species(species_name):
        if species_name == 'e-': return "gas", None, "e-"
        m = re.match(r'^([JK])(\d+)(.+)', species_name)
        if m: return ("surface" if m.group(1) == 'J' else "mantle"), m.group(2), m.group(3).replace('c-', '').replace('l-', '')
        return "gas", None, species_name.replace('c-', '').replace('l-', '')

    def count_target_atom(clean_formula, target):
        if clean_formula == 'e-': return 0
        f = clean_formula.split('+')[0].split('-')[0]
        return sum(int(n) if n else 1 for atom, n in re.findall(r'([A-Z][a-z]*)(\d*)', f) if atom == target)

    radii_map = {int(re.findall(r'\d+', str(k))[0]): k for k in main_output_dict.keys() if re.findall(r'\d+', str(k))}
    radii = sorted([r for r in radii_map.keys() if (rmin is None or r >= rmin) and (rmax is None or r <= rmax)])

    radial_plot_data = {}
    all_encountered = set()

    for r_val in radii:
        sub_dict = main_output_dict[radii_map[r_val]]
        try:
            ab_arr = sub_dict['abundances']
            nH_mid = sub_dict["H_number_density"][itime, -1]
            y_ab = ab_arr.isel(time=itime, spatial=-1).values if 'spatial' in ab_arr.dims else ab_arr.isel(time=itime)[..., -1].values

            local_cont = {}
            for idx, sp in enumerate(ab_arr.coords['species'].values):
                if "GRAIN" in sp: continue
                p, b, f = parse_species(sp)
                if phase != "all" and (phase == "grain" and p not in ["surface", "mantle"] or phase != "grain" and p != phase): continue
                if grain_bin is not None and b != str(grain_bin): continue

                coef = count_target_atom(f, target_atom)
                if coef > 0:
                    tot_atoms = y_ab[idx] * nH_mid * coef
                    if tot_atoms > 0: local_cont[f] = local_cont.get(f, 0.0) + tot_atoms

            r_sum = sum(local_cont.values())
            if r_sum > 0:
                top_entries = sorted({sp: (v/r_sum)*100 for sp, v in local_cont.items()}.items(), key=lambda x: x[1], reverse=True)[:spnumber]
                radial_plot_data[r_val] = top_entries
                for sp, _ in top_entries: all_encountered.add(sp)
        except Exception as e:
            if verbose: print(f"Error processing midplane at R={r_val}: {e}")

    if not radial_plot_data: return
    unique_sp = sorted(list(all_encountered))
    cmap = plt.get_cmap(cmap_name)
    species_colors = {sp: cmap(idx/max(1, len(unique_sp)-1)) for idx, sp in enumerate(unique_sp)}

    y_labels, x_pct, bar_colors = [], [], []
    for r_val in sorted(list(radial_plot_data.keys()), reverse=True):
        for sp, pct in reversed(radial_plot_data[r_val]):
            y_labels.append(f"{r_val} AU — ${re.sub(r'(\d+)', r'_{\1}', sp)}$")
            x_pct.append(pct)
            bar_colors.append(species_colors[sp])

    fig, ax = plt.subplots(figsize=(11, max(4, len(y_labels)*0.45)))
    ax.barh(y_labels, x_pct, color=bar_colors, edgecolor='grey', alpha=0.85, height=0.7)
    ax.set_xlabel('Disk Midplane Budget Contribution (%)')
    plt.tight_layout()
    plt.show()


def plot_species_evolution_with_grain_size_comparison(PIPE, MODEL_NAMES, target_radius, itime=-1,
                                                      verbose=True, spnumber=5, cmap_name="tab10"):
    """
    Plots the structural size distribution tracking breakdown of key carrier reservoirs without residual mapping columns.
    """
    radii_list = [target_radius] if not isinstance(target_radius, list) else target_radius
    radii_list = [int(r) for r in radii_list]

    if len(MODEL_NAMES) != len(set(MODEL_NAMES)): MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    def get_grain_size_in_um(file_path, bin_index):
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('!'): continue
                    if '!' in line: line = line.split('!')[0].strip()
                    values = [float(val) for val in line.split()]
                    if not values: continue
                    num_grains = len(values) // 4
                    return values[:num_grains][int(bin_index) - 1] * 10000.0
            return None
        except: return None

    def parse_species(species_name):
        if species_name == 'e-': return "gas", None, "e-"
        grain_match = re.match(r'^([JK])(\d+)(.+)', species_name)
        if grain_match:
            p_code, g_bin, raw_formula = grain_match.groups()
            return ("surface" if p_code == 'J' else "mantle"), g_bin, raw_formula.replace('c-', '').replace('l-', '')
        return "gas", None, species_name.replace('c-', '').replace('l-', '')

    def clean_molec(mol_name):
        raw = re.sub(r"^[JK]\d+", "", mol_name)
        f = re.sub(r"(\d+)", r"_{\1}", raw)
        f = re.sub(r"([+-]+)$", r"^{\1}", f)
        return f"${f}$"

    AU_to_cm = 1.496e13
    model_data = {}
    all_encountered_species = set()
    time_years_string = None

    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        main_output_dict = p.chemistry
        chempath = Path(p.chempath)
        
        radii_map = {int(re.findall(r'\d+', str(k))[0]): k for k in main_output_dict.keys() if re.findall(r'\d+', str(k))}
        sorted_all_radii = sorted(list(radii_map.keys()))
        model_data[p_name] = {}

        for r_value in radii_list:
            if r_value not in radii_map: continue
            sub_dict = main_output_dict[radii_map[r_value]]
            abundance_array = sub_dict['abundances']
            nH_profile = sub_dict["H_number_density"][itime,:]
            raw_species_list = list(abundance_array.coords['species'].values)

            r_idx = sorted_all_radii.index(r_value)
            if len(sorted_all_radii) > 1:
                r_midshifts = 0.5 * np.diff(sorted_all_radii)
                r_edges = [sorted_all_radii[0] - r_midshifts[0]] + [sorted_all_radii[i] + r_midshifts[i] for i in range(len(r_midshifts))] + [sorted_all_radii[-1] - r_midshifts[-1]]
            else:
                r_edges = [r_value * 0.9, r_value * 1.1]

            file_path = os.path.join(chempath, f"{r_value}AU", "1D_static.dat")
            z_points = np.loadtxt(file_path, comments='!', usecols=0)
            if len(z_points) > 1:
                z_midshifts = 0.5 * np.diff(z_points)
                z_edges = [z_points[0] - z_midshifts[0]] + [z_points[j] + z_midshifts[j] for j in range(len(z_midshifts))] + [max(0.0, z_points[-1] + z_midshifts[-1])]
                dz = np.abs(np.diff(z_edges))
            else:
                dz = np.array([z_points[0] if z_points[0] > 0 else 1.0])

            cell_volumes = 2 * np.pi * (float(r_value) * AU_to_cm) * ((r_edges[r_idx+1] - r_edges[r_idx]) * AU_to_cm) * (dz * AU_to_cm)
            nH_2d, volumes_2d = nH_profile[np.newaxis, :], cell_volumes[np.newaxis, :]

            available_bins_set = set()
            species_metadata = {}
            for sp in raw_species_list:
                if "GRAIN" in sp: continue
                sp_phase, sp_bin, clean_formula = parse_species(sp)
                if sp_phase in ["surface", "mantle"] and sp_bin is not None:
                    available_bins_set.add(sp_bin)
                    species_metadata[sp] = (sp_phase, sp_bin, clean_formula)

            sorted_bins = sorted(list(available_bins_set), key=lambda x: int(x)) if available_bins_set else []
            if not sorted_bins: continue

            bin_raw_data = {b: {} for b in sorted_bins}
            local_species_scores = {}

            for sp in raw_species_list:
                if sp not in species_metadata: continue
                _, sp_bin, clean_formula = species_metadata[sp]
                absolute_particles = np.sum(abundance_array.sel(species=sp).values * nH_2d * volumes_2d)
                if absolute_particles > 0:
                    bin_raw_data[sp_bin][clean_formula] = bin_raw_data[sp_bin].get(clean_formula, 0.0) + absolute_particles
                    local_species_scores[clean_formula] = local_species_scores.get(clean_formula, 0.0) + absolute_particles

            sorted_local = sorted(local_species_scores.items(), key=lambda x: x[1], reverse=True)
            top_local_species = [item[0] for item in sorted_local[:spnumber]]
            for sp in top_local_species: all_encountered_species.add(sp)

            if time_years_string is None:
                try: time_years_string = f"{abundance_array.coords['time'].values[itime] / 3.156e7:.0f}"
                except: pass

            model_data[p_name][r_value] = {'bin_raw_data': bin_raw_data, 'sorted_bins': sorted_bins, 'top_species': top_local_species, 'chempath': chempath}

    num_models = len(PIPE)
    model_names = list(model_data.keys())
    num_radii = len(radii_list)

    cols = num_models
    rows = num_radii
    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 5 * rows), squeeze=False, sharey=True)

    unique_species_list = sorted(list(all_encountered_species))
    cmap = plt.colormaps.get_cmap(cmap_name)
    species_colors = {sp: cmap(idx / max(1, len(unique_species_list)-1)) for idx, sp in enumerate(unique_species_list)}

    for row_idx, r_value in enumerate(radii_list):
        for col_idx, p_name in enumerate(model_names):
            if r_value not in model_data[p_name]: continue
            ax = axes[row_idx, col_idx]
            struct = model_data[p_name][r_value]
            bin_raw, sbins, local_top, cpath = struct['bin_raw_data'], struct['sorted_bins'], struct['top_species'], struct['chempath']

            shared_grain_sizes_um = []
            for b in sbins:
                size_um = get_grain_size_in_um(cpath / f"{r_value}AU" / "1D_grain_sizes.in", b)
                shared_grain_sizes_um.append(f"{size_um:.1f}" if size_um is not None else f"B{b}")

            x_positions = np.arange(len(sbins))
            for sp in local_top:
                plot_pct = []
                for b in sbins:
                    tot = sum(bin_raw[b].values())
                    plot_pct.append((bin_raw[b].get(sp, 0.0) / tot * 100) if tot > 0 else 0.0)

                ax.plot(x_positions, plot_pct, label=clean_molec(sp), color=species_colors[sp], lw=1.8, marker='o', ms=4)

            ax.set_xlabel('Grain Radius [µm]', fontsize=11)
            ax.set_ylabel('Contribution (%)', fontsize=11)
            ax.set_title(f"{p_name} @ {r_value} AU", fontsize=11, fontweight='bold')
            ax.set_xticks(x_positions)
            ax.set_xticklabels(shared_grain_sizes_um, fontsize=8, rotation=60)
            ax.set_ylim(-2, 105)
            ax.grid(True, linestyle="--", alpha=0.4)
            ax.legend(loc='upper right', ncol=2, fontsize=9)

    fig.suptitle(f'Top Carriers vs Grain Size Distribution — $t = {time_years_string}$ years' if time_years_string else 'Top Ice Carriers vs Grain Size', fontsize=15, y=0.99)
    plt.tight_layout()
    plt.show()


def plot_ratio_midplane_gas_vs_grain(chempath, main_output_dict, s1="C", s2="O", itime=-1, starratio=None, verbose=True, xlim=None, ylim=None):
    """
    Plots the comparative midplane abundance ratio of two specific elements splitting gas vs aggregate solid grain phases.
    """
    if s1 not in ['H','He','C','N','O','Si','S','Fe','Na','Mg','Cl','P','F'] or s2 not in ['H','He','C','N','O','Si','S','Fe','Na','Mg','Cl','P','F']:
        raise ValueError("Elements not supported.")

    def parse_and_count(sp_name, element1, element2):
        if sp_name == 'e-' or 'GRAIN' in sp_name: return "ignore", 0, 0
        m = re.match(r'^([JK])\d+(.+)', sp_name)
        f = m.group(2).replace('c-', '').replace('l-', '').split('+')[0].split('-')[0] if m else sp_name.replace('c-', '').replace('l-', '').split('+')[0].split('-')[0]
        composition = {atom: int(n) if n else 1 for atom, n in re.findall(r'([A-Z][a-z]*)(\d*)', f)}
        return ("grain" if m else "gas"), composition.get(element1, 0), composition.get(element2, 0)

    radii_map = {int(re.findall(r'\d+', str(k))[0]): k for k in main_output_dict.keys() if re.findall(r'\d+', str(k))}
    radii = sorted(list(radii_map.keys()))

    r_list, gas_r, grain_r = [], [], []
    for r in radii:
        ab = main_output_dict[radii_map[r]]['abundances'].isel(time=itime, spatial=-1).values
        species = list(main_output_dict[radii_map[r]]['abundances'].coords['species'].values)
        
        g1, g2, s1_tot, s2_tot = 0.0, 0.0, 0.0, 0.0
        for idx, sp in enumerate(species):
            phase, c1, c2 = parse_and_count(sp, s1, s2)
            if phase == "gas": g1 += ab[idx]*c1; g2 += ab[idx]*c2
            elif phase == "grain": s1_tot += ab[idx]*c1; s2_tot += ab[idx]*c2

        r_list.append(float(r))
        gas_r.append(g1/g2 if g2>0 else 0.0)
        grain_r.append(s1_tot/s2_tot if s2_tot>0 else 0.0)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(r_list, gas_r, color="teal", marker='o', label='Gas')
    ax.plot(r_list, grain_r, color="darkred", linestyle='--', marker='s', label='Grains (Ice)')
    ax.set_xlabel('Radius R [AU]')
    ax.set_ylabel(f'Midplane Atomic Ratio [{s1}/{s2}]')
    if starratio: ax.axhline(starratio, color='black', linestyle=':')
    ax.grid(True, linestyle=':')
    ax.legend()
    plt.show()


def plot_grain_properties_midplane_comparison(PIPE, MODEL_NAMES, key_list=['CO'], itime=-1, fracab=True,
                                             verbose=True, xlim=None, ylim=None, Tmin=None, Tmax=None, vmin=None, vmax=None,
                                             temp_colormap='hot', ab_colormap='plasma', common_scale=True, species_scale_common=False):
    """
    Renders 2D grain profile snapshots matching each individual model dynamically to a dedicated column space.
    """
    if isinstance(key_list, str): key_list = [key_list]
    key_list = list(dict.fromkeys(key_list))
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)): MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]

    def parse_grain_sizes_midplane(file_path):
        try:
            valid_lines = []
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('!'): continue
                    if '!' in line: line = line.split('!')[0].strip()
                    values = [float(val) for val in line.split()]
                    if values: valid_lines.append(values)
            midplane_values = valid_lines[-1]
            num_grains = len(midplane_values) // 4
            return [r * 10000.0 for r in midplane_values[:num_grains]], midplane_values[2 * num_grains : 3 * num_grains]
        except: return None, None

    def clean_molec(mol_name):
        return f"${re.sub(r'([+-]+)$', r'^{\1}', re.sub(r'(\d+)', r'_{\1}', mol_name))}$"

    model_data = {}
    time_years_string = None

    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        main_output_dict = p.chemistry
        chempath = Path(p.chempath)
        
        radii_map = {int(re.findall(r'\d+', str(k))[0]): k for k in main_output_dict.keys() if re.findall(r'\d+', str(k))}
        extracted_radii = sorted(list(radii_map.keys()))
        
        disk_radii, grain_sizes_matrix, grain_temps_matrix = [], [], []
        species_abundance_matrices = {key: [] for key in key_list}

        for r_val in extracted_radii:
            radii_um, temps_k = parse_grain_sizes_midplane(chempath / f"{r_val}AU" / "1D_grain_sizes.in")
            if radii_um is None: continue
            
            disk_radii.append(r_val)
            grain_sizes_matrix.append(radii_um)
            grain_temps_matrix.append(temps_k)
            
            if time_years_string is None:
                try: time_years_string = f"{p.chemistry[radii_map[r_val]]['abundances'].coords['time'].values[itime] / 3.156e7:.0f}"
                except: pass

            for key in key_list:
                bin_values = []
                ab_array = main_output_dict[radii_map[r_val]]['abundances']
                for b_idx in range(1, len(radii_um) + 1):
                    v_cell = 0.0
                    s_name, m_name = f"J{b_idx:02d}{key}", f"K{b_idx:02d}{key}"
                    if s_name in ab_array.coords['species'].values: v_cell += float(ab_array.isel(time=itime).sel(species=s_name).values[-1])
                    if m_name in ab_array.coords['species'].values: v_cell += float(ab_array.isel(time=itime).sel(species=m_name).values[-1])
                    if not fracab: v_cell *= main_output_dict[radii_map[r_val]]["H_number_density"][itime, -1]
                    bin_values.append(v_cell)
                species_abundance_matrices[key].append(bin_values)

        model_data[p_name] = {
            'disk_radii': np.array(disk_radii), 'grain_sizes': np.array(grain_sizes_matrix),
            'grain_temps': np.array(grain_temps_matrix), 'abundance_matrices': species_abundance_matrices
        }

    num_models = len(PIPE)
    model_names = list(model_data.keys())
    rows = 1 + len(key_list)
    cols = num_models
    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.0 * rows), squeeze=False, sharex=True, sharey=True)

    num_grain_bins = model_data[model_names[0]]['grain_sizes'].shape[1]
    y_centers = np.arange(num_grain_bins)
    y_edges = np.arange(num_grain_bins + 1) - 0.5

    # Render Dust Temperatures
    if common_scale:
        all_temps = np.concatenate([model_data[m]['grain_temps'].flatten() for m in model_names])
        actual_tmin, actual_tmax = Tmin if Tmin else all_temps.min(), Tmax if Tmax else all_temps.max()

    for col_idx, p_name in enumerate(model_names):
        ax = axes[0, col_idx]
        struct = model_data[p_name]
        radii = struct['disk_radii']
        temps = struct['grain_temps']
        
        grid_R, grid_Y = np.meshgrid(np.linspace(radii.min(), radii.max(), 200), np.linspace(y_edges[0], y_edges[-1], 200))
        points_R = [r for r in radii for _ in y_centers]
        points_Y = [y for _ in radii for y in y_centers]
        grid_T = griddata((points_R, points_Y), temps.flatten(), (grid_R, grid_Y), method='cubic')
        
        cf = ax.contourf(grid_R, grid_Y, grid_T, levels=np.linspace(actual_tmin, actual_tmax, 50), cmap=temp_colormap)
        fig.colorbar(cf, ax=ax, label=r"$T_{\rm grain}$ [K]")
        ax.set_title(f"{p_name}\nGrain Temp $T_{\rm grain}$", fontsize=11, fontweight='bold')

    # Render Chemical Species Profiles
    for row_idx, key in enumerate(key_list):
        current_row = row_idx + 1
        for col_idx, p_name in enumerate(model_names):
            ax = axes[current_row, col_idx]
            struct = model_data[p_name]
            radii = struct['disk_radii']
            ab_matrix = struct['abundance_matrices'][key]
            
            r_edges = [radii[0]-0.5, radii[0]+0.5] if len(radii)==1 else [radii[0]-0.5*np.diff(radii)[0]] + [radii[i]+0.5*np.diff(radii)[i] for i in range(len(np.diff(radii)))] + [radii[-1]+0.5*np.diff(radii)[-1]]
            polygons, values = [], []
            for i in range(len(radii)):
                for j in range(num_grain_bins):
                    polygons.append([(r_edges[i], y_edges[j]), (r_edges[i+1], y_edges[j]), (r_edges[i+1], y_edges[j+1]), (r_edges[i], y_edges[j+1])])
                    values.append(ab_matrix[i][j])
            
            vals_array = np.array(values)
            v_min_loc = vals_array[vals_array>0].min() if len(vals_array[vals_array>0])>0 else 1e-15
            v_max_loc = vals_array.max()
            norm = LogNorm(vmin=v_min_loc, vmax=v_max_loc) if v_max_loc/v_min_loc > 10 else Normalize(vmin=v_min_loc, vmax=v_max_loc)
            
            coll = PolyCollection(polygons, array=vals_array, cmap=ab_colormap, norm=norm, edgecolors='none')
            ax.add_collection(coll)
            fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=ab_colormap), ax=ax, label="Ice Abundance")
            ax.set_title(f"{p_name}\n{clean_molec(key)} Ice Profile", fontsize=11)

    for ax in axes.flatten():
        ax.set_yticks(y_centers)
        ax.set_yticklabels([f"{s:.1f}" for size_row in model_data[model_names[0]]['grain_sizes'] for s in size_row][:num_grain_bins])
        ax.set_xlim(xlim if xlim else (model_data[model_names[0]]['disk_radii'].min(), model_data[model_names[0]]['disk_radii'].max()))
        ax.set_ylim(ylim if ylim else (y_edges[0], y_edges[-1]))
        ax.set_xlabel('Radius R [AU]')
        ax.set_ylabel('Grain Size [µm]')

    fig.suptitle(f'Midplane Grain Properties — $t = {time_years_string}$ years' if time_years_string else 'Midplane Grain Properties', fontsize=13, fontweight='bold', y=0.99)
    plt.tight_layout()
    plt.show()


def density2D_grid_comparison(PIPE, MODEL_NAMES, vmin=1e-30, vmax=1e-15, cmap='gnuplot2',
                              dens_type='mass', xlim=None, ylim=None, dust=None, select_bins=None, figsize=(14, 16)):
    """
    Compares 2D density structures across multiple runs, arranging each model run into an isolated, sequential column slice layout.
    """
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)): MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]
    model_structures = {}

    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        base_path = str(Path(p.thermalpath)) + '/'
        
        grid_vals = pd.read_table(base_path + 'amr_grid.inp', engine='python', skiprows=5)
        nr = int(grid_vals.columns[0].split()[0])
        nt = int(grid_vals.columns[0].split()[1])
        g_vals = grid_vals[grid_vals.columns[0]].values

        dens_vals = pd.read_table(base_path + 'dust_density.inp', engine='python', header=None, skiprows=3)[0].values
        nspecies = int(len(dens_vals) / (nr * nt))
        dens_tensor = np.reshape(dens_vals, (nspecies, nt, nr))

        r_edge = g_vals[:nr+1] / autocm
        theta_edge = g_vals[nr+1:nr+1+nt+1].copy()
        theta_edge[-1] = np.pi
        rr_edge, tt_edge = np.meshgrid(r_edge, theta_edge)
        
        dens_tensor[dens_tensor <= 1e-100] = 1e-100
        sizes = np.atleast_1d(np.loadtxt(base_path + 'dust_sizes.inp')) if os.path.isfile(base_path + 'dust_sizes.inp') else None

        valid_bins = select_bins if select_bins else list(range(nspecies))
        model_structures[p_name] = {
            'R': rr_edge * np.sin(tt_edge), 'Z': rr_edge * np.cos(tt_edge),
            'dens': dens_tensor[valid_bins], 'sizes': sizes[valid_bins] if sizes is not None else None,
            'nspecies_display': len(valid_bins), 'original_bins': valid_bins
        }

    num_models = len(model_structures)
    model_names = list(model_structures.keys())
    nspecies_ref = model_structures[model_names[0]]['nspecies_display']

    nrows = nspecies_ref + 1
    cols = num_models
    fig, axes = plt.subplots(nrows, cols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)

    cbar_main_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])

    for row_idx in range(nspecies_ref + 1):
        for col_idx, p_name in enumerate(model_names):
            ax = axes[row_idx, col_idx]
            struct = model_structures[p_name]
            
            if row_idx < nspecies_ref:
                plot_data = struct['dens'][row_idx]
                sz = struct['sizes'][row_idx] if struct['sizes'] is not None else None
                title_str = f"bin {struct['original_bins'][row_idx]+1} ({f'{sz/1e3:.1f} mm' if sz>=1e3 else f'{sz:.2f} µm'})" if sz else f"bin {struct['original_bins'][row_idx]+1}"
            else:
                plot_data = struct['dens'].sum(axis=0)
                title_str = f"Total {dens_type.title()}"

            im = ax.pcolormesh(struct['R'], struct['Z'], plot_data, cmap=cmap, shading='auto', norm=LogNorm(vmin=vmin, vmax=vmax))
            ax.set_title(title_str, fontsize=11)
            ax.grid(True, linestyle=":", alpha=0.3)
            
            if row_idx == 0:
                ax.text(0.5, 1.2, p_name, transform=ax.transAxes, fontsize=12, fontweight='bold', ha='center')

    fig.colorbar(im, cax=cbar_main_ax, label=f"Density Field [{dens_type}]")
    for ax in axes[-1, :]: ax.set_xlabel('r [au]', fontsize=12)
    for ax in axes[:, 0]: ax.set_ylabel('z [au]', fontsize=12)
    fig.subplots_adjust(right=0.90, left=0.08, bottom=0.06, top=0.90, hspace=0.35, wspace=0.15)
    plt.show()


def temperature2D_grid_comparison(PIPE, MODEL_NAMES, vmin=1.0, vmax=1e3, cmap='gnuplot2',
                                  xlim=None, ylim=None, snowline_temp=20.0, select_bins=None, figsize=(14, 12), dust=None):
    """
    Compares 2D poloidal dust temperature cross-sections. Each model occupies its own individual subplot column mapping index.
    """
    if len(MODEL_NAMES) != len(set(MODEL_NAMES)): MODEL_NAMES = [f"Model {i+1}" for i in range(len(PIPE))]
    model_structures = {}

    for p_idx, p in enumerate(PIPE):
        p_name = getattr(p, 'name', MODEL_NAMES[p_idx])
        base_path = str(Path(p.thermalpath)) + '/'
        
        grid_vals = pd.read_table(base_path + 'amr_grid.inp', engine='python', skiprows=5)
        nr = int(grid_vals.columns[0].split()[0])
        nt = int(grid_vals.columns[0].split()[1])
        g_vals = grid_vals[grid_vals.columns[0]].values

        temp_vals = pd.read_table(base_path + 'dust_temperature.dat', engine='python', header=None, skiprows=3)[0].values
        nspecies = int(len(temp_vals) / (nr * nt))
        temp_tensor = np.reshape(temp_vals, (nspecies, nt, nr))

        r_edge = g_vals[:nr+1] / autocm
        theta_edge = g_vals[nr+1:nr+1+nt+1].copy()
        theta_edge[-1] = np.pi
        rr_edge, tt_edge = np.meshgrid(r_edge, theta_edge)
        
        R = rr_edge * np.sin(tt_edge)
        Z = rr_edge * np.cos(tt_edge)
        R_center = 0.25 * (R[:-1, :-1] + R[1:, :-1] + R[:-1, 1:] + R[1:, 1:])
        Z_center = 0.25 * (Z[:-1, :-1] + Z[1:, :-1] + Z[:-1, 1:] + Z[1:, 1:])

        sizes = np.atleast_1d(np.loadtxt(base_path + 'dust_sizes.inp')) if os.path.isfile(base_path + 'dust_sizes.inp') else None
        valid_bins = select_bins if select_bins else list(range(nspecies))

        model_structures[p_name] = {
            'R': R, 'Z': Z, 'R_center': R_center, 'Z_center': Z_center,
            'temp': temp_tensor[valid_bins], 'sizes': sizes[valid_bins] if sizes is not None else None,
            'nspecies_display': len(valid_bins), 'original_bins': valid_bins
        }

    num_models = len(model_structures)
    model_names = list(model_structures.keys())
    nspecies_ref = model_structures[model_names[0]]['nspecies_display']

    fig, axes = plt.subplots(nspecies_ref, num_models, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    cbar_main_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])

    for row_idx in range(nspecies_ref):
        for col_idx, p_name in enumerate(model_names):
            ax = axes[row_idx, col_idx]
            struct = model_structures[p_name]
            plot_data = struct['temp'][row_idx]
            
            sz = struct['sizes'][row_idx] if struct['sizes'] is not None else None
            title_str = f"bin {struct['original_bins'][row_idx]+1} ({f'{sz/1e3:.1f} mm' if sz>=1e3 else f'{sz:.2f} µm'})" if sz else f"bin {struct['original_bins'][row_idx]+1}"

            im = ax.pcolormesh(struct['R'], struct['Z'], plot_data, cmap=cmap, shading='auto', norm=LogNorm(vmin=vmin, vmax=vmax))
            ax.set_title(title_str, fontsize=11)
            ax.grid(True, linestyle=":", alpha=0.3)

            if snowline_temp:
                try: ax.contour(struct['R_center'], struct['Z_center'], plot_data, levels=[float(snowline_temp)], colors='black', linewidths=1.5)
                except: pass

            if row_idx == 0:
                ax.text(0.5, 1.2, p_name, transform=ax.transAxes, fontsize=12, fontweight='bold', ha='center')

    fig.colorbar(im, cax=cbar_main_ax, label=r'Dust Temperature $T_\mathrm{d}$ [K]')
    for ax in axes[-1, :]: ax.set_xlabel('r [au]', fontsize=11)
    for ax in axes[:, 0]: ax.set_ylabel('z [au]', fontsize=11)
    fig.subplots_adjust(right=0.90, left=0.08, bottom=0.08, top=0.90, hspace=0.35, wspace=0.22)
    plt.show()
