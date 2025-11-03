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

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from astromugs.constants.constants import autocm


def density2D(mass1, mass2=None, overlap=False):
    grid = pd.read_table('thermal/amr_grid.inp', engine='python', skiprows=5)
    head = grid.columns
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[head[0]].values
    dens = pd.read_table('thermal/dust_density.inp', engine='python', header=None, skiprows=3)
    dens = dens[0].values
    nbspecies = int(len(dens)/(nr*nt))
    dens = np.reshape(dens, (nbspecies, nt, nr))
    dist = grid[:nr+1]/autocm
    theta = grid[nr+1:nr+1+nt+1]
    theta[-1] = np.pi
    dist, tt = np.meshgrid(dist, theta)
    rr = dist*np.sin(tt)
    zz = dist*np.cos(tt)

    dens[dens<=1e-100] = 1e-100
    if overlap == False:
        # #--PLOT FIGURE--
        if nbspecies == 1:
            fig = plt.figure(figsize=(10, 8.))
            ax = fig.add_subplot(111)
            plt.xlabel(r'r [au]', fontsize = 17)
            plt.ylabel(r'z [au]', fontsize = 17, labelpad=-7.4)
            
            numdens = dens[0]#/mass1[0]
            t = plt.pcolor(rr, zz, numdens, cmap='gnuplot2', shading='auto', norm=LogNorm(vmin=1e-80, vmax=1e-1))
            clr = plt.colorbar(t)
            clr.set_label(r'$n_\mathrm{d}$ [cm${-3}$]', labelpad=-33, y=1.06, rotation=0, fontsize = 16)
            #plt.xlim(1, 500)
            #plt.ylim(-300, 300)
            ax.tick_params(labelsize=17)
            clr.ax.tick_params(labelsize=16) 
            plt.show()

        else:
            for ispec in range(nbspecies):
                fig = plt.figure(figsize=(8, 8.))
                ax = fig.add_subplot(111)
                plt.xlabel(r'r [au]', fontsize = 17)
                plt.ylabel(r'z [au]', fontsize = 17, labelpad=-7.4)
                numdens = dens[ispec]#/mass1[ispec]
                t = plt.pcolor(rr, zz, numdens, cmap='gnuplot2', shading='auto', norm=LogNorm(vmin=1e-30, vmax=1e-17))
                clr = plt.colorbar(t)
                clr.set_label(r'$\rho_\mathrm{d}$ [g.cm${-3}$]', labelpad=-33, y=1.06, rotation=0, fontsize = 16)
                plt.xlim(1, 1000)
                plt.ylim(-1000, 1000)
                ax.tick_params(labelsize=17)
                clr.ax.tick_params(labelsize=16) 
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
                ax.text(0.90, 0.95, 'bin: {}'.format(ispec+1), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=16, bbox=props)
                plt.show()

    if overlap == True:
        density = np.zeros((nt, nr))
        fig = plt.figure(figsize=(10, 10.))
        ax = fig.add_subplot(111)
        plt.xlabel(r'r [au]', fontsize = 17)
        plt.ylabel(r'z [au]', fontsize = 17, labelpad=-7.4)
        numdens = dens#/mass1[ispec]
        for ispec in range(0, nbspecies):
            density += numdens[ispec]
        t = plt.pcolor(rr, zz, density, cmap='gnuplot2', shading='auto', norm=LogNorm(vmin=1e-25, vmax=1e-17))
        clr = plt.colorbar(t)
        clr.set_label(r'$\rho_\mathrm{d}$ [g.cm${-3}$]', labelpad=-33, y=1.06, rotation=0, fontsize = 16)
        #plt.xlim(1, 500)
        #plt.ylim(-300, 300)
        ax.tick_params(labelsize=17)
        clr.ax.tick_params(labelsize=16) 
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.90, 0.95, 'env+disk', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=16, bbox=props)
        plt.show()


def temperature2D():
    grid = pd.read_table('thermal/amr_grid.inp', engine='python', skiprows=5)
    head = grid.columns
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[head[0]].values
    temp = pd.read_table('thermal/dust_temperature.dat', engine='python', header=None, skiprows=3)
    temp = temp[0].values
    nbspecies = int(len(temp)/(nr*nt))
    temp = np.reshape(temp, (nbspecies, nt, nr))
    dist = grid[:nr+1]/autocm
    dist = 0.5*(dist[0:dist.size-1] + dist[1:dist.size])
    theta = grid[nr+1:nr+1+nt+1]
    theta[-1] = np.pi
    theta = 0.5*(theta[0:theta.size-1] + theta[1:theta.size])
    dist, tt = np.meshgrid(dist, theta)
    rr = dist*np.sin(tt)
    zz = dist*np.cos(tt)

    #--PLOT FIGURE--
    if nbspecies == 1:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        plt.xlabel(r'r [au]', fontsize = 17)
        plt.ylabel(r'z [au]', fontsize = 17)
        t = plt.pcolormesh(rr, zz, temp[0], cmap='hot', shading='gouraud', vmin=18, vmax=70)
        clr = plt.colorbar(t)
        clr.set_label(r'$T_\mathrm{d}$ [K]', labelpad=-33, y=1.06, rotation=0, fontsize = 16)
        ax.tick_params(labelsize=17)
        clr.ax.tick_params(labelsize=16) 
        plt.xlim(1, 200)
        plt.ylim(-100, 100)
        plt.show()
    else:
        for ispec in range(nbspecies):
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111)
            plt.xlabel(r'r [au]', fontsize = 17)
            plt.ylabel(r'z [au]', fontsize = 17)
            t = plt.pcolormesh(rr, zz, temp[ispec], cmap='hot', shading='gouraud', vmin=18, vmax=70)
            clr = plt.colorbar(t)
            clr.set_label(r'$T_\mathrm{d}$ [K]', labelpad=-33, y=1.06, rotation=0, fontsize = 16)
            ax.tick_params(labelsize=17)
            clr.ax.tick_params(labelsize=16)
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.90, 0.95, 'bin: {}'.format(ispec+1), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=16, bbox=props)
            plt.xlim(1, 200)
            plt.ylim(-100, 100)
            plt.show()

def midplane_temp():
    grid = pd.read_table('thermal/amr_grid.inp', engine='python', skiprows=5)
    head = grid.columns
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[head[0]].values
    try:
        temp = pd.read_table('thermal/dust_temperature.dat', engine='python', header=None, skiprows=3)
    except IOError:
        print('plot.midplane_temp: the file thermal/dust_temperature.dat is not present. Run a dust thermal simulation first.')
        sys.exit(1)
    temp = temp[0].values
    nbspecies = int(len(temp)/(nr*nt))
    temp = np.reshape(temp, (nbspecies, nt, nr))
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
    ax.set_ylim(0,60)
    #ax.set_xlim(1,350)
    ax.set_xlabel(r'r [au]', fontsize = 20)
    ax.set_ylabel(r'T [K]', fontsize = 20)
    ax.legend(fontsize=15)
    ax.tick_params(labelsize=18)
    plt.show()

def vertical_temp(r=100):
    grid = pd.read_table('thermal/amr_grid.inp', engine='python', skiprows=5)
    head = grid.columns
    nr = int(grid.columns[0].split("  ")[0])
    nt = int(grid.columns[0].split("  ")[1])
    grid = grid[head[0]].values
    temp = pd.read_table('thermal/dust_temperature.dat', engine='python', header=None, skiprows=3)
    temp = temp[0].values
    nbspecies = int(len(temp)/(nr*nt))

    if nbspecies == 1:
        try:
            temp = pd.read_table('chemistry/'+str(r)+'AU/1D_static.dat', sep=r"\s+", engine='python', header=None, comment='!')
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
            static = pd.read_table('chemistry/'+str(r)+'AU/1D_static.dat', sep=r"\s+", engine='python', header=None, comment='!')
            temp = pd.read_table('chemistry/'+str(r)+'AU/temperatures.dat', sep=r"\s+", engine='python', header=None)
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

def avz(r=100):
    static = pd.read_table('chemistry/'+str(r)+'AU/1D_static.dat', sep=r"\s+", engine='python', header=None, comment='!', skiprows=1)
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

def opacity():
    opaclist = sorted(glob.glob('thermal/dustkap*'))

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


def localflux():
    #---1/ Get grid shape and reshape the local flux array accordingly
    flux = pd.read_table('thermal/mean_intensity.out', sep=r"\s+", comment='#', header=None, skiprows=4)
    grid = pd.read_table('thermal/amr_grid.inp', engine='python', skiprows=5)
    lam = pd.read_table('thermal/mcmono_wavelength_micron.inp', engine='python', header=None, skiprows=1)
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
