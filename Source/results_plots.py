#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 10:47:46 2022

results plot functions for intercomparrison

@author: tjor
"""

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.dates as mdates

import dataframe_image as dfi


import numpy as np
import pandas as pd
import datetime 
import os

from scipy.interpolate import interp1d
import read_IP_data as rd_IP

def scatter_subplot(spec_type, system, plot_index, ylab, xlab, limits, ticks, df_sys, df_R,bands):
    ''' suplot routine for scatter plot'''
   
    colors = cm.rainbow(np.linspace(0,1,10))# color mask to match rrs with time series   
    plt.subplot(2,4,plot_index)    
    X = np.arange(0,2000,1)
    plt.plot(X,X,color='gray')
    plt.title(system)
    
    # fill up
    if spec_type == 'nLw':
        for i in range(1,8,1):
            plt.scatter(df_R[str(bands[i])], df_sys[str(bands[i])], color=colors[i], facecolors='none',label=str(bands[i]) + ' nm')  
    else :
        for i in range(10):
            plt.scatter(df_R[str(bands[i])], df_sys[str(bands[i])], color=colors[i], facecolors='none',label=str(bands[i]) + ' nm')
    
    if plot_index == 1: 
        plt.legend(fontsize=10)
    
    # ables
    plt.gca().set_aspect('equal')
    plt.xlim(limits)
    plt.ylim(limits)
    plt.xticks(ticks)
    plt.yticks(ticks)
    
    if plot_index < 4:
        plt.xticks() 

    if spec_type == 'Rrs':    
        plt.xticks(ticks, rotation=45)
        plt.yticks(ticks)
    
    return

def scatter_subplot_CP(spec_type, system, plot_index, ylab, xlab, limits, ticks, df_sys, df_R,bands):
    ''' suplot routine for scatter plot - CP version'''
   
    colors = cm.rainbow(np.linspace(0,1,10))# color mask to match rrs with time series   
    plt.subplot(1,4,plot_index)    
    X = np.arange(0,2000,1)
    plt.plot(X,X,color='gray')
    plt.title(system)
    
    # fill up
    if spec_type == 'nLw':
        for i in range(1,8,1):
            plt.scatter(df_R[str(bands[i])], df_sys[str(bands[i])], color=colors[i], facecolors='none',label=str(bands[i]) + ' nm')  
    else :
        for i in range(10):
            plt.scatter(df_R[str(bands[i])], df_sys[str(bands[i])], color=colors[i], facecolors='none',label=str(bands[i]) + ' nm')
    
    if plot_index == 1: 
        plt.legend(fontsize=10)
    
    # ables
    plt.gca().set_aspect('equal')
    plt.xlim(limits)
    plt.ylim(limits)
    plt.xticks(ticks)
    plt.yticks(ticks)
    
    #if plot_index < 4:
     #   plt.xticks() 

    if spec_type == 'Rrs':    
        plt.xticks(ticks, rotation=45)
        plt.yticks(ticks)
    
    return

def plot_scatter(spec_type,df_R, df_PML, df_NASA, df_TARTU, df_HEREON, df_RBINS, df_CNR, df_NOAA, bands, path_output,  Q_mask, Qtype = 'AOC_3'):
    ''' scatter plot'''
  
    # qc filter
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
    df_RBINS = df_RBINS[Q_mask[Qtype]==1]
    df_CNR = df_CNR[Q_mask[Qtype]==1]
    df_NOAA = df_NOAA[Q_mask[Qtype]==1]

    df_R = df_R[Q_mask[Qtype]==1]
    
    # scatter plot figure
    fig = plt.figure(figsize=(24,12))
    plt.rc('font', size=18)      
    if spec_type == 'Ed':
        xlab ='Reference: $E_{d}^{r}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$]'
        ylab = '$E_{d}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$]'
        limits = [800, 1650]
        ticks = [800, 1000, 1200,1400, 1600]
        plt.suptitle('System inter-comparison for downwelling irradiance: $E_{d}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Lsky':
        xlab ='Reference: $L_{sky}^{r}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = '$L_{sky}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 130]
        ticks = [0, 30, 60, 90, 120]
        plt.suptitle('System inter-comparison for sky radiance: $L_{sky}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Lt':
        xlab ='Reference: $L_{t}^{r}$(0$^{+})$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = '$L_{t}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 25]
        ticks = [0, 5, 10, 15, 20 , 25]
        plt.suptitle('System inter-comparison for upwelling radiance: $L_{t}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Rrs':
        xlab ='Reference: $R_{rs}^{r}$($\lambda)$ [sr$^{-1}$]'
        ylab = '$R_{rs}$($\lambda)$ [sr$^{-1}$]'
        limits = [0, 0.016]
        ticks = [0, 0.004, 0.008, 0.012, 0.016]
        plt.suptitle('System inter-comparison for remote-sensing reflectance: $R_{rs}$($\lambda$)')
    elif spec_type == 'nLw':
        xlab ='Reference: $L_{wn}^{r}$($\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = '$L_{nw}$($\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 25]
        ticks = [0, 5, 10, 15, 20 , 25]
        plt.suptitle('System inter-comparison for normalized water-leaving radiance: $L_{wn}$($\lambda$)')
       
    fig.supxlabel(xlab)
    fig.supylabel(ylab)
        
    # subplots
    subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML['400'])))
    index = 1
    scatter_subplot(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_PML, df_R ,bands)
    
    subtitle = 'HEREON: N = '  + str(np.sum(~np.isnan(df_HEREON['400'])))
    index = 2
    scatter_subplot(spec_type,subtitle, index,  ylab, xlab, limits, ticks, df_HEREON, df_R,bands)
    
    subtitle = 'TARTU: N = '  + str(np.sum(~np.isnan(df_TARTU['400'])))
    index = 3
    scatter_subplot(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_TARTU, df_R,bands)
    
    subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
    index = 4
    scatter_subplot(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_NASA, df_R,bands)

    subtitle = 'RBINS: N = ' + str(np.sum(~np.isnan(df_RBINS['400'])))
    index = 5
    scatter_subplot(spec_type,subtitle,index, ylab, xlab, limits, ticks, df_RBINS, df_R,bands)
    
    subtitle = 'CNR: N = ' + str(np.sum(~np.isnan(df_CNR['400'])))
    index = 6
    scatter_subplot(spec_type,subtitle,index, ylab, xlab, limits, ticks, df_CNR, df_R, bands)
    
    if spec_type =='Ed' or spec_type == 'nLw':
        subtitle = 'NOAA: N = ' + str(np.sum(~np.isnan(df_NOAA['400'])))
        index = 7
        scatter_subplot(spec_type,subtitle,index, ylab, xlab, limits, ticks, df_NOAA, df_R, bands)

    plt.tight_layout(pad=1.8)

    
    filename  =  path_output +  '/' + spec_type + '_scattterplot.png'
    plt.savefig(filename)
    
    return

def plot_scatter_CP(spec_type,df_R, df_PML, df_NASA, df_TARTU, df_HEREON, bands, path_output,  Q_mask, Qtype = 'AOC_3'):
    '''  scatter plot - CP version'''
  
    # qc filter
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]

    df_R = df_R[Q_mask[Qtype]==1]
    
    # scatter plot figure
    fig = plt.figure(figsize=(20,6))    
    plt.rc('font', size=18)      
    if spec_type == 'Ed':
        xlab ='Reference: $E_{d}^{r}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$]'
        ylab = '$E_{d}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$]'
        limits = [800, 1650]
        ticks = [800, 1000, 1200,1400, 1600]
        plt.suptitle('System inter-comparison for downwelling irradiance: $E_{d}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Lsky':
        xlab ='Reference: $L_{sky}^{r}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = '$L_{sky}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 130]
        ticks = [0, 30, 60, 90, 120]
        plt.suptitle('System inter-comparison for sky radiance: $L_{sky}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Lt':
        xlab ='Reference: $L_{t}^{r}$(0$^{+})$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = '$L_{t}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 25]
        ticks = [0, 5, 10, 15, 20 , 25]
        plt.suptitle('System inter-comparison for upwelling radiance: $L_{t}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Rrs':
        xlab ='Reference: $R_{rs}^{r}$($\lambda)$ [sr$^{-1}$]'
        ylab = '$R_{rs}$($\lambda)$ [sr$^{-1}$]'
        limits = [0, 0.016]
        ticks = [0, 0.004, 0.008, 0.012, 0.016]
        plt.suptitle('System inter-comparison for remote-sensing reflectance: $R_{rs}$($\lambda$)')
    elif spec_type == 'nLw':
        xlab ='Reference: $L_{wn}^{r}$($\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = '$L_{nw}$($\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 25]
        ticks = [0, 5, 10, 15, 20 , 25]
        plt.suptitle('System inter-comparison for normalized water-leaving radiance: $L_{wn}$($\lambda$)')
       
    fig.supxlabel(xlab)
    fig.supylabel(ylab)
        
    # subplots
    subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML['400'])))
    index = 1
    scatter_subplot_CP(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_PML, df_R ,bands)
    
    subtitle = 'HEREON: N = ' + str(np.sum(~np.isnan(df_HEREON['400'])))
    index = 2
    scatter_subplot_CP(spec_type,subtitle, index,  ylab, xlab, limits, ticks, df_HEREON, df_R,bands)
    
    subtitle = 'TARTU: N = ' + str(np.sum(~np.isnan(df_TARTU['400'])))
    index = 3
    scatter_subplot_CP(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_TARTU, df_R,bands)
    
    subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
    index = 4
    scatter_subplot_CP(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_NASA, df_R,bands)
    
        
    plt.tight_layout(pad=1.8) 
    
    filename  =  path_output +  '/' + spec_type + '_scattterplot_CP.png'
    plt.savefig(filename)
    
    return


def plot_scatter_IP(spec_type,df_R, df_PML, df_NASA, df_TARTU, df_HEREON, bands, path_output,  Q_mask, Qtype = 'AOC_3'):
    '''  scatter plot - IP version with 4 teams'''
  
    # qc filter
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]

    df_R = df_R[Q_mask[Qtype]==1]
    
    # scatter plot figure
    fig = plt.figure(figsize=(20,6))    
    plt.rc('font', size=18)      
    if spec_type == 'Ed':
        xlab ='Reference: $E_{d}^{r}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$]'
        ylab = '$E_{d}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$]'
        limits = [800, 1650]
        ticks = [800, 1000, 1200,1400, 1600]
        plt.suptitle('System inter-comparison for downwelling irradiance: $E_{d}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Lsky':
        xlab ='Reference: $L_{sky}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = '$L_{sky}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 130]
        ticks = [0, 30, 60, 90, 120]
        plt.suptitle('System inter-comparison for sky radiance: $L_{sky}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Lt':
        xlab ='Reference: $L_{t}$(0$^{+})$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = '$L_{t}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 25]
        ticks = [0, 5, 10, 15, 20 , 25]
        plt.suptitle('System inter-comparison for upwelling radiance: $L_{t}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Rrs':
        xlab ='Reference: $R_{rs}$($\lambda)$ [sr$^{-1}$]'
        ylab = '$R_{rs}$($\lambda)$ [sr$^{-1}$]'
        limits = [0, 0.016]
        ticks = [0, 0.004, 0.008, 0.012, 0.016]
        plt.suptitle('System inter-comparison for remote-sensing reflectance: $R_{rs}$($\lambda$)')
    elif spec_type == 'nLw':
        xlab ='Reference: $L_{wn}$($\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = '$L_{nw}$($\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 25]
        ticks = [0, 5, 10, 15, 20 , 25]
        plt.suptitle('System inter-comparison for normalized water-leaving radiance: $L_{wn}$($\lambda$)')
       
    fig.supxlabel(xlab)
    fig.supylabel(ylab)
        
    # subplots
    subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML['400'])))
    index = 1
    scatter_subplot_CP(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_PML, df_R ,bands)
    
    subtitle = 'HEREON: N = ' + str(np.sum(~np.isnan(df_HEREON['400'])))
    index = 2
    scatter_subplot_CP(spec_type,subtitle, index,  ylab, xlab, limits, ticks, df_HEREON, df_R,bands)
    
    subtitle = 'TARTU: N = ' + str(np.sum(~np.isnan(df_TARTU['400'])))
    index = 3
    scatter_subplot_CP(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_TARTU, df_R,bands)
    
    subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
    index = 4
    scatter_subplot_CP(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_NASA, df_R,bands)
    
        
    plt.tight_layout(pad=1.8) 
    
    filename  =  path_output +  '/' + spec_type + '_scattterplot_IP.png'
    plt.savefig(filename)
    
    return


def plot_scatter_CP_vs_IP(spec_type, df_PML_IP, df_NASA_IP, df_TARTU_IP, df_HEREON_IP, df_PML_CP, df_NASA_CP, df_TARTU_CP, df_HEREON_CP, bands, path_output,  Q_mask, Qtype = 'AOC_3'):
    '''  scatter plot - variant that plots IP versus CP processed data'''
  
    # qc filter
    df_PML_IP = df_PML_IP[Q_mask[Qtype]==1]
    df_NASA_IP = df_NASA_IP[Q_mask[Qtype]==1]
    df_TARTU_IP = df_TARTU_IP[Q_mask[Qtype]==1]
    df_HEREON_IP = df_HEREON_IP[Q_mask[Qtype]==1]
    
    df_PML_CP = df_PML_CP[Q_mask[Qtype]==1]
    df_NASA_CP = df_NASA_CP[Q_mask[Qtype]==1]
    df_TARTU_CP = df_TARTU_CP[Q_mask[Qtype]==1]
    df_HEREON_CP = df_HEREON_CP[Q_mask[Qtype]==1]
    
    # scatter plot figure
    fig = plt.figure(figsize=(22,6))
    plt.rc('font', size=16)      
    if spec_type == 'Ed':
        xlab ='CP: $E_{d}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$]'
        ylab = 'IP: $E_{d}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$]'
        limits = [800, 1650]
        ticks = [800, 1000, 1200,1400, 1600]
        plt.suptitle('Individual versus Community Processor: $E_{d}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Lsky':
        xlab ='CP: $L_{sky}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = 'IP: $L_{sky}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 130]
        ticks = [0, 30, 60, 90, 120]
        plt.suptitle('Individual versus Community Processor: $L_{sky}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Lt':
        xlab ='CP: $L_{t}$(0$^{+})$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = 'IP: $L_{t}$(0$^{+}$, $\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 25]
        ticks = [0, 5, 10, 15, 20 , 25]
        plt.suptitle('Individual versus Community Processor: $L_{t}$(0$^{+}$,$\lambda$)')
    elif spec_type == 'Rrs':
        xlab ='CP: $R_{rs}$($\lambda)$ [sr$^{-1}$]'
        ylab = 'IP: $R_{rs}$($\lambda)$ [sr$^{-1}$]'
        limits = [0, 0.016]
        ticks = [0, 0.004, 0.008, 0.012, 0.016]
        plt.suptitle('Individual versus Community Processor: $R_{rs}$(0$^{+}$, $\lambda$)')
    elif spec_type == 'nLw':
        xlab ='CP: $L_{n}$($\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        ylab = 'IP: $L_{nw}$($\lambda)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
        limits = [0, 25]
        ticks = [0, 5, 10, 15, 20 , 25]
        plt.suptitle('Individual versus Community Processor: $L_{wn}$($\lambda$)')
       
    fig.supxlabel(xlab)
    fig.supylabel(ylab)
        
    # subplots
    subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML_CP['400'])))
    index = 1
    scatter_subplot_CP(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_PML_IP, df_PML_CP,bands)
    
    subtitle = 'HEREON: N = ' + str(np.sum(~np.isnan(df_HEREON_CP['400'])))
    index = 2
    scatter_subplot_CP(spec_type,subtitle, index,  ylab, xlab, limits, ticks, df_HEREON_IP, df_HEREON_CP,bands)
    
    subtitle = 'TARTU: N = ' + str(np.sum(~np.isnan(df_TARTU_CP['400'])))
    index = 3
    scatter_subplot_CP(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_TARTU_IP, df_TARTU_CP,bands)
    
    subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA_CP['400'])))
    index = 4
    scatter_subplot_CP(spec_type,subtitle, index, ylab, xlab, limits, ticks, df_NASA_IP, df_NASA_CP,bands)
    
    filename  =  path_output +  '/' + spec_type + '_scattterplot_CP_vs_IP.png'
    plt.savefig(filename)
    
    plt.tight_layout(pad=2.2) 
    
    return



def _resid_subplot(spec_type,system, plot_index, ylab, percent_limits, df_sys, df_R ,bands):
    ''' suplot routine for residuals'''  
    colors = cm.rainbow(np.linspace(0,1,10)) 
    
    if spec_type == 'nLw':
        resid = [] # res
        resid.append([])
        for i in range(1,8,1):# residual disrtibution in each band
            resid_i = 100*np.array((df_sys[str(bands[i])] -  df_R[str(bands[i])])/df_R[str(bands[i])])
            resid_i = resid_i[~np.isnan(resid_i)]
            resid.append(resid_i)
        resid.append([])
        resid.append([])
    else:
        resid = [] # residual distribution in each band
        for i in range(10):
            resid_i = 100*np.array((df_sys[str(bands[i])] -  df_R[str(bands[i])])/df_R[str(bands[i])])
            resid_i =  resid_i[~np.isnan(resid_i)]
            resid.append(resid_i)

    #
    plt.subplot(2,4,plot_index) 
    plt.title(system)
    plt.plot(np.arange(0,12,1), np.zeros(12),linestyle ='dashed',color='gray')
    bp = plt.boxplot(resid ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90]) 
    plt.xlim(0.5,10.5)
    plt.ylim(-percent_limits, percent_limits)
    for i in range(10):
        bp['boxes'][i].set_facecolor(colors[i])
   
    plt.xticks([1,2,3,4,5,6,7,8,9,10], bands[0:10])
    plt.xticks(rotation = 45)
    
    plt.grid(axis='y') 
    
    #if plot_index==1 or plot_index== 5:
     #   plt.ylabel(ylab)
    #if plot_index > 3:
     #   plt.xlabel('Wavelength [nm]')
  
    return


def _resid_subplot_CP_IP(spec_type,system, plot_index, ylab, percent_limits, df_sys, df_R ,bands):
    ''' suplot routine for residuals'''  
    colors = cm.rainbow(np.linspace(0,1,10)) 
    
    if spec_type == 'nLw':
        resid = [] # res
        resid.append([])
        for i in range(1,8,1):# residual disrtibution in each band
            resid_i = 100*np.array((df_sys[str(bands[i])] -  df_R[str(bands[i])])/df_R[str(bands[i])])
            resid_i = resid_i[~np.isnan(resid_i)]
            resid.append(resid_i)
        resid.append([])
        resid.append([])
    else:
        resid = [] # residual distribution in each band
        for i in range(10):
            resid_i = 100*np.array((df_sys[str(bands[i])] -  df_R[str(bands[i])])/df_R[str(bands[i])])
            resid_i =  resid_i[~np.isnan(resid_i)]
            resid.append(resid_i)

    #
    plt.subplot(1,4,plot_index) 
    plt.title(system)
    plt.plot(np.arange(0,12,1), np.zeros(12),linestyle ='dashed',color='gray')
    bp = plt.boxplot(resid ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90]) 
    plt.xlim(0.5,10.5)
    plt.ylim(-percent_limits, percent_limits)
    for i in range(10):
        bp['boxes'][i].set_facecolor(colors[i])
   
    plt.xticks([1,2,3,4,5,6,7,8,9,10], bands[0:10])
    plt.xticks(rotation = 45)
    
    plt.grid(axis='y') 
    
    #if plot_index==1 or plot_index== 5:
     #   plt.ylabel(ylab)
    #if plot_index > 3:
     #   plt.xlabel('Wavelength [nm]')
  
    return


def _resid_subplot_CP(spec_type,system, plot_index, ylab, percent_limits, df_sys, df_sys_unc ,df_R ,bands):
    ''' suplot routine for residuals'''  
    colors = cm.rainbow(np.linspace(0,1,10)) 
    
    if spec_type == 'nLw':
        resid = [] # res
        resid.append([])
        for i in range(1,8,1):# residual disrtibution in each band
            resid_i = 100*np.array((df_sys[str(bands[i])] -  df_R[str(bands[i])])/df_R[str(bands[i])])
            resid_i = resid_i[~np.isnan(resid_i)]
            resid.append(resid_i)
        resid.append([])
        resid.append([])
    else:
        resid = [] # residual distribution in each band
        for i in range(10):
            resid_i = 100*np.array((df_sys[str(bands[i])] -  df_R[str(bands[i])])/df_R[str(bands[i])])
            resid_i =  resid_i[~np.isnan(resid_i)]
            resid.append(resid_i)

    print(resid)
    plt.subplot(1,4,plot_index) 
    plt.title(system)
    #plt.plot(np.arange(0,12,1), np.zeros(12),linestyle ='dashed',color='gray')
    plt.plot(np.arange(1,11,1), df_sys_unc[0:10],linestyle ='dashed',color='gray',label = 'Median unc')
    plt.plot(np.arange(1,11,1), -df_sys_unc[0:10],linestyle ='dashed',color='gray')
    #if plot_index == 1:  # tunred off for nlw
       # plt.legend()
    bp = plt.boxplot(resid ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90]) 
    plt.xlim(0.5,10.5)
    plt.ylim(-percent_limits, percent_limits)
    for i in range(10):
        bp['boxes'][i].set_facecolor(colors[i])
        
   
    plt.xticks([1,2,3,4,5,6,7,8,9,10], bands[0:10])
    plt.xticks(rotation = 45)
    
    plt.grid(axis='y') 
    
    #if plot_index==1 or plot_index== 5:
     #   plt.ylabel(ylab)
    #if plot_index > 3:
     #   plt.xlabel('Wavelength [nm]')
  
    return



def _unc_subplot_CP(spec_type,system, plot_index, ylab, percent_limits, df_sys ,bands):
    ''' suplot routine for residuals'''  
    colors = cm.rainbow(np.linspace(0,1,10)) 
    
    if spec_type == 'nLw':
        resid = [] # res
        resid.append([])
        for i in range(1,8,1):# residual disrtibution in each band
            resid_i = np.array((df_sys[str(bands[i])]))
            resid_i = resid_i[~np.isnan(resid_i)]
            resid.append(resid_i)
        resid.append([])
        resid.append([])
    else:
        resid = [] # residual distribution in each band
        for i in range(10):
            resid_i =  np.array((df_sys[str(bands[i])]))
            resid_i =  resid_i[~np.isnan(resid_i)]
            resid.append(resid_i)

    plt.subplot(1,4,plot_index) 
    plt.title(system)
    plt.plot(np.arange(0,12,1), np.zeros(12),linestyle ='dashed',color='gray')
    bp = plt.boxplot(resid ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90]) 
    plt.xlim(0.5,10.5)
    plt.ylim(0, percent_limits)
    for i in range(10):
        bp['boxes'][i].set_facecolor(colors[i])
   
    plt.xticks([1,2,3,4,5,6,7,8,9,10], bands[0:10])
    plt.xticks(rotation = 45)
    
    plt.grid(axis='y') 
    
    # if plot_index==1 or plot_index== 5:
    #   plt.ylabel(ylab)
    # if plot_index > 3:
    #   plt.xlabel('Wavelength [nm]')
  
    return


def plot_residuals(spec_type, df_R, df_PML, df_NASA, df_TARTU, df_HEREON, df_RBINS, df_CNR, df_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3'):
    ''' Funtion to plot spectral dependence of % residuals following Tilstone 2020'''  
    
    # QC filtering 
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
    df_RBINS = df_RBINS[Q_mask[Qtype]==1]
    df_CNR = df_CNR[Q_mask[Qtype]==1]
    df_NOAA = df_NOAA[Q_mask[Qtype]==1]
    
    
    df_R = df_R[Q_mask[Qtype]==1]
    
    # spectral reiduals plot
    fig= plt.figure(figsize=(18,12))
    plt.rc('font',size=16)  
    if spec_type == 'Ed':
        ylab = '$E_{d}$ residual [%]'
        plt.suptitle('Percentage residuals for downwelling irradiance: $E_{d}$(0$^{+}$,$\lambda$)')
        percent_limits = 10
        percent_limits_2 = 10
    if spec_type == 'Lsky':
         ylab = '$L_{sky}$ residual [%]'
         plt.suptitle('Percentage residuals for sky radiance: $L_{sky}$(0$^{+}$,$\lambda$)')
         percent_limits = 10
         percent_limits_2 = 10
    if spec_type == 'Lt':
         ylab = '$L_{t}$ residual [%]'
         plt.suptitle('Percentage residuals for upwelling radiance: $L_{t}$(0$^{+}$,$\lambda$)')
         percent_limits = 10
         percent_limits_2 = 10
    if spec_type == 'Rrs':
          ylab = '$R_{rs}$ residual [%]'
          plt.suptitle('Percentage residuals for remote-sensing reflectance: $R_{rs}$($\lambda$)')
          percent_limits = 16
          percent_limits_2 = 16
    if spec_type == 'nLw':
          ylab = '$L_{wn}$ residual [%]'
          percent_limits = 40
          percent_limits_2 = 40
          plt.suptitle('Percentage residuals for normalized water-leaving radiance: $L_{wn}$($\lambda$)')

    xlab = 'Wavelength [nm]'
    fig.supxlabel(xlab)
    fig.supylabel(ylab)
        
    subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML['400'])))
    index = 1
    _resid_subplot(spec_type, subtitle, index, ylab, percent_limits, df_PML, df_R ,bands)
    
    subtitle  = 'HEREON: N = ' + str(np.sum(~np.isnan(df_HEREON['400'])))
    index = 2
    _resid_subplot(spec_type,subtitle, index, ylab,percent_limits, df_HEREON, df_R ,bands)
    
    subtitle  = 'TARTU: N = ' + str(np.sum(~np.isnan(df_TARTU['400'])))
    index = 3
    _resid_subplot(spec_type,subtitle, index, ylab, percent_limits, df_TARTU, df_R ,bands)

   
    if spec_type == 'Lsky' :

        subtitle  = 'RBINS: N = ' + str(np.sum(~np.isnan(df_RBINS['400'])))
        index = 5
        _resid_subplot(spec_type,subtitle, index, ylab,percent_limits_2, df_RBINS, df_R ,bands)
        
        subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
        index = 4
        _resid_subplot(spec_type, subtitle, index, ylab, percent_limits_2, df_NASA, df_R ,bands)
        
        subtitle  = 'CNR: N = ' + str(np.sum(~np.isnan(df_CNR['400'])))
        index = 6
        _resid_subplot(spec_type, subtitle, index, ylab, percent_limits_2, df_CNR, df_R ,bands)
   
    elif spec_type == 'Lt' :
        
        subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
        index = 4
        _resid_subplot(spec_type, subtitle, index, ylab, percent_limits, df_NASA, df_R ,bands)
        
        subtitle  = 'RBINS: N = ' + str(np.sum(~np.isnan(df_RBINS['400'])))
        index = 5
        _resid_subplot(spec_type,subtitle, index, ylab,percent_limits_2, df_RBINS, df_R ,bands)
    
                
        subtitle  = 'CNR: N = ' + str(np.sum(~np.isnan(df_CNR['400'])))
        index = 6
        _resid_subplot(spec_type,subtitle, index, ylab,percent_limits_2, df_CNR, df_R ,bands)
        
    elif spec_type == 'Rrs' :
          
        subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
        index = 4
        _resid_subplot(spec_type, subtitle, index, ylab, percent_limits, df_NASA, df_R ,bands)
        
        subtitle  = 'RBINS: N = ' + str(np.sum(~np.isnan(df_RBINS['400'])))
        index = 5
        _resid_subplot(spec_type,subtitle, index, ylab,percent_limits_2, df_RBINS, df_R ,bands)
    
        subtitle  = 'CNR: N = ' + str(np.sum(~np.isnan(df_CNR['400'])))
        index = 6
        _resid_subplot(spec_type,subtitle, index, ylab,percent_limits_2, df_CNR, df_R ,bands)
   
    else: 

      subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
      index = 4
      _resid_subplot(spec_type, subtitle, index, ylab, percent_limits, df_NASA, df_R ,bands)
      
      subtitle  = 'RBINS: N = ' + str(np.sum(~np.isnan(df_RBINS['400'])))
      index = 5
      _resid_subplot(spec_type,subtitle, index, ylab,percent_limits, df_RBINS, df_R ,bands)
    
      subtitle  = 'CNR: N = ' + str(np.sum(~np.isnan(df_CNR['400'])))
      index = 6
      _resid_subplot(spec_type,subtitle, index, ylab,percent_limits, df_CNR, df_R ,bands)

        
    if spec_type == 'Ed' or  spec_type == 'nLw':
        subtitle  = 'NOAA: N = ' + str(np.sum(~np.isnan(df_NOAA['400'])))
        index = 7
        _resid_subplot(spec_type,subtitle, index, ylab, percent_limits_2, df_NOAA, df_R ,bands)
        
    
    plt.tight_layout()
    
    filename  =  path_output +  '/' + spec_type + '_resiudalsplot.png'
    plt.savefig(filename)
    
    return


def plot_residuals_CP(spec_type, df_R, df_PML, df_NASA, df_TARTU, df_HEREON, df_unc_PML, df_unc_NASA, df_unc_TARTU, df_unc_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3'):
    ''' Funtion to plot spectral dependence of % residuals following Tilstone 2020'''  
    
    # QC filtering 
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
    
    df_unc_PML = df_unc_PML[Q_mask[Qtype]==1]
    df_unc_NASA =  df_unc_NASA[Q_mask[Qtype]==1]
    df_unc_TARTU = df_unc_TARTU[Q_mask[Qtype]==1]
    df_unc_HEREON = df_unc_HEREON[Q_mask[Qtype]==1]
    
    unc_med_PML = np.nanmedian(df_unc_PML.iloc[:,1:-1],0)
    unc_med_NASA = np.nanmean(df_unc_NASA.iloc[:,1:-1],0)
    unc_med_TARTU = np.nanmedian(df_unc_TARTU.iloc[:,1:-1],0)
    unc_med_HEREON = np.nanmedian(df_unc_HEREON.iloc[:,1:-1],0)
    
    df_R = df_R[Q_mask[Qtype]==1]
    
    # spectral reiduals plot
    fig= plt.figure(figsize=(18,6))
    plt.rc('font',size=16)  
    if spec_type == 'Ed':
        ylab = '$E_{d}$ residual [%]'
        plt.suptitle('Community Processor: Percentage residuals for downwelling irradiance: $E_{d}$(0$^{+}$,$\lambda$)')
        percent_limits = 8
    if spec_type == 'Lsky':
         ylab = '$L_{sky}$ residual [%]'
         plt.suptitle('Community Processor: Percentage residuals for sky radiance: $L_{sky}$(0$^{+}$,$\lambda$)')
         percent_limits = 10
    if spec_type == 'Lt':
         ylab = '$L_{t}$ residual [%]'
         plt.suptitle('Community Processor: Percentage residuals for upwelling radiance: $L_{t}$(0$^{+}$,$\lambda$)')
         percent_limits = 10
    if spec_type == 'Rrs':
          ylab = '$R_{rs}$ residual [%]'
          plt.suptitle('Community Processor: Percentage residuals for remote-sensing reflectance: $R_{rs}$($\lambda$)')
          percent_limits = 20
    if spec_type == 'nLw':
          ylab = '$L_{wn}$ residual [%]'
          percent_limits = 16
          plt.suptitle('Community Processor: Percentage residuals for normalized water-leaving radiance: $L_{wn}$($\lambda$)')

    xlab = 'Wavelength [nm]'
    fig.supxlabel(xlab)
    fig.supylabel(ylab)
        
    subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML['400'])))
    index = 1
    _resid_subplot_CP(spec_type, subtitle, index, ylab, percent_limits, df_PML, unc_med_PML, df_R, bands)
    
    subtitle  = 'HEREON: N = ' + str(np.sum(~np.isnan(df_HEREON['400'])))
    index = 2
    _resid_subplot_CP(spec_type,subtitle, index, ylab,percent_limits, df_HEREON, unc_med_HEREON, df_R, bands)
    
    subtitle  = 'TARTU: N = ' + str(np.sum(~np.isnan(df_TARTU['400'])))
    index = 3
    _resid_subplot_CP(spec_type,subtitle, index, ylab, percent_limits, df_TARTU, unc_med_TARTU, df_R, bands)

    subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
    index = 4
    _resid_subplot_CP(spec_type, subtitle, index, ylab, percent_limits, df_NASA, unc_med_NASA, df_R, bands)
        
    plt.tight_layout()
    
    filename  =  path_output +  '/' + spec_type + '_residualsplot_CP.png'
    plt.savefig(filename)
    
    
    return


def plot_residuals_IP(spec_type, df_R, df_PML, df_NASA, df_TARTU, df_HEREON, df_unc_PML, df_unc_NASA, df_unc_TARTU, df_unc_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3'):
    ''' Funtion to plot spectral dependence of % residuals following Tilstone 2020'''  
    
    # QC filtering 
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
    
    df_unc_PML = df_unc_PML[Q_mask[Qtype]==1]
    df_unc_NASA =  df_unc_NASA[Q_mask[Qtype]==1]
    df_unc_TARTU = df_unc_TARTU[Q_mask[Qtype]==1]
    df_unc_HEREON = df_unc_HEREON[Q_mask[Qtype]==1]
    
    unc_med_PML = np.nanmedian(df_unc_PML.iloc[:,1:-1],0)
    unc_med_NASA = np.nanmean(df_unc_NASA.iloc[:,1:-1],0)
    unc_med_TARTU = np.nanmedian(df_unc_TARTU.iloc[:,1:-1],0)
    unc_med_HEREON = np.nanmedian(df_unc_HEREON.iloc[:,1:-1],0)
        
    df_R = df_R[Q_mask[Qtype]==1]
    
    # spectral reiduals plot
    fig= plt.figure(figsize=(18,6))
    plt.rc('font',size=16)  
    if spec_type == 'Ed':
        ylab = '$E_{d}$ residual [%]'
        plt.suptitle('Individual Processor: Percentage residuals for downwelling irradiance: $E_{d}$(0$^{+}$,$\lambda$)')
        percent_limits = 8
    if spec_type == 'Lsky':
         ylab = '$L_{sky}$ residual [%]'
         plt.suptitle('Individual Processor: Percentage residuals for sky radiance: $L_{sky}$(0$^{+}$,$\lambda$)')
         percent_limits = 10
    if spec_type == 'Lt':
         ylab = '$L_{t}$ residual [%]'
         plt.suptitle('Individual Processor: Percentage residuals for upwelling radiance: $L_{t}$(0$^{+}$,$\lambda$)')
         percent_limits = 10
    if spec_type == 'Rrs':
          ylab = '$R_{rs}$ residual [%]'
          plt.suptitle('Individual Processor: Percentage residuals for remote-sensing reflectance: $R_{rs}$(0$^{+}$, $\lambda$)')
          percent_limits = 20
    if spec_type == 'nLw':
          ylab = '$L_{wn}$ residual [%]'
          percent_limits = 16
          plt.suptitle('Individual Processor: Percentage residuals for normalized water-leaving radiance: $L_{wn}$(0$^{+}$, $\lambda$)')

    xlab = 'Wavelength [nm]'
    fig.supxlabel(xlab)
    fig.supylabel(ylab)
        
    subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML['400'])))
    index = 1
    _resid_subplot_CP(spec_type, subtitle, index, ylab, percent_limits, df_PML, unc_med_PML, df_R, bands)
    
    subtitle  = 'HEREON: N = ' + str(np.sum(~np.isnan(df_HEREON['400'])))
    index = 2
    _resid_subplot_CP(spec_type,subtitle, index, ylab,percent_limits, df_HEREON, unc_med_HEREON, df_R, bands)
    
    subtitle  = 'TARTU: N = ' + str(np.sum(~np.isnan(df_TARTU['400'])))
    index = 3
    _resid_subplot_CP(spec_type,subtitle, index, ylab, percent_limits, df_TARTU, unc_med_TARTU, df_R, bands)

    subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
    index = 4
    _resid_subplot_CP(spec_type, subtitle, index, ylab, percent_limits, df_NASA, unc_med_NASA, df_R, bands)
        
    plt.tight_layout()
    
    filename  =  path_output +  '/' + spec_type + '_residualsplot_IP.png'
    plt.savefig(filename)
    
    
    return


def plot_residuals_IP_nounc(spec_type, df_R, df_PML, df_NASA, df_TARTU, df_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3'):
    ''' Funtion to plot spectral dependence of % residuals following Tilstone 2020'''  
    
    # QC filtering 
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
    
    # df_unc_PML = df_unc_PML[Q_mask[Qtype]==1]
    # df_unc_NASA =  df_unc_NASA[Q_mask[Qtype]==1]
    # df_unc_TARTU = df_unc_TARTU[Q_mask[Qtype]==1]
    # df_unc_HEREON = df_unc_HEREON[Q_mask[Qtype]==1]
    
    unc_med_PML = np.nan*np.ones(17)
    unc_med_NASA =  np.nan*np.ones(17)
    unc_med_TARTU =  np.nan*np.ones(17)
    unc_med_HEREON =  np.nan*np.ones(17)     
    df_R = df_R[Q_mask[Qtype]==1]
    
    # spectral reiduals plot
    fig= plt.figure(figsize=(18,6))
    plt.rc('font',size=16)  
    if spec_type == 'Ed':
        ylab = '$E_{d}$ residual [%]'
        plt.suptitle('Individual Processor: Percentage residuals for downwelling irradiance: $E_{d}$(0$^{+}$,$\lambda$)')
        percent_limits = 8
    if spec_type == 'Lsky':
         ylab = '$L_{sky}$ residual [%]'
         plt.suptitle('Individual Processor: Percentage residuals for sky radiance: $L_{sky}$(0$^{+}$,$\lambda$)')
         percent_limits = 10
    if spec_type == 'Lt':
         ylab = '$L_{t}$ residual [%]'
         plt.suptitle('Individual Processor: Percentage residuals for upwelling radiance: $L_{t}$(0$^{+}$,$\lambda$)')
         percent_limits = 10
    if spec_type == 'Rrs':
          ylab = '$R_{rs}$ residual [%]'
          plt.suptitle('Individual Processor: Percentage residuals for remote-sensing reflectance: $R_{rs}$(0$^{+}$, $\lambda$)')
          percent_limits = 20
    if spec_type == 'nLw':
          ylab = '$L_{wn}$ residual [%]'
          percent_limits = 40
          plt.suptitle('Individual Processor: Percentage residuals for normalized water-leaving radiance: $L_{wn}$(0$^{+}$, $\lambda$)')

    xlab = 'Wavelength [nm]'
    fig.supxlabel(xlab)
    fig.supylabel(ylab)
        
    subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML['400'])))
    index = 1
    _resid_subplot_CP(spec_type, subtitle, index, ylab, percent_limits, df_PML, unc_med_PML, df_R, bands)
    
    subtitle  = 'HEREON: N = ' + str(np.sum(~np.isnan(df_HEREON['400'])))
    index = 2
    _resid_subplot_CP(spec_type,subtitle, index, ylab,percent_limits, df_HEREON, unc_med_HEREON, df_R, bands)
    
    subtitle  = 'TARTU: N = ' + str(np.sum(~np.isnan(df_TARTU['400'])))
    index = 3
    _resid_subplot_CP(spec_type,subtitle, index, ylab, percent_limits, df_TARTU, unc_med_TARTU, df_R, bands)

    subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
    index = 4
    _resid_subplot_CP(spec_type, subtitle, index, ylab, percent_limits, df_NASA, unc_med_NASA, df_R, bands)
        
    plt.tight_layout()
    
    filename  =  path_output +  '/' + spec_type + '_residualsplot_IP.png'
    plt.savefig(filename)
    
    
    return


def plot_residuals_CP_vs_IP(spec_type, df_PML, df_NASA, df_TARTU, df_HEREON, df_PML_CP, df_NASA_CP, df_TARTU_CP, df_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3'):
    ''' version of residuals that shows CP verus IP processed data'''  
    
    # QC filtering 
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
    
    df_PML_CP = df_PML_CP[Q_mask[Qtype]==1]
    df_NASA_CP = df_NASA_CP[Q_mask[Qtype]==1]
    df_TARTU_CP = df_TARTU_CP[Q_mask[Qtype]==1]
    df_HEREON_CP = df_HEREON_CP[Q_mask[Qtype]==1]
    
    # spectral reiduals plot
    fig= plt.figure(figsize=(18,6))
    plt.rc('font',size=18)  
    if spec_type == 'Ed':
        ylab = '$E_{d}$ residual, $\Delta$ [%]'
        plt.suptitle('IP vs CP: Percentage residuals for $E_{d}$: $\Delta = 100(X_{IP} - X_{CP})/X_{CP}$')
        percent_limits = 6
    if spec_type == 'Lsky':
         ylab = '$L_{sky}$ residual, $\Delta$  [%]'
         plt.suptitle('IP vs CP: Percentage residuals for $L_{sky}$: $\Delta = 100(X_{IP} - X_{CP})/X_{CP}$')
         percent_limits = 6
    if spec_type == 'Lt':
         ylab = '$L_{t}$ residual, $\Delta$ [%]'
         plt.suptitle('IP vs CP: Percentage residuals for $L_{t}$: $\Delta = 100(X_{IP} - X_{CP})/X_{CP}$')
         percent_limits = 6
    if spec_type == 'Rrs':
          ylab = '$R_{rs}$ residual, $\Delta$  [%]'
          plt.suptitle('IP vs CP: Percentage residuals for $R_{rs}$: $\Delta = 100(X_{IP} - X_{CP})/X_{CP}$')  
          percent_limits = 20
    if spec_type == 'nLw':
          ylab = '$L_{wn}$ residual, $\Delta$  [%]'
          percent_limits = 40
          plt.suptitle('IP vs CP: Percentage residuals for $L_{nw}$: $\Delta = 100(X_{IP} - X_{CP})/X_{CP}$')

    xlab = 'Wavelength [nm]'
    fig.supxlabel(xlab)
    fig.supylabel(ylab)
        
    subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML_CP['400'])))
    index = 1
    _resid_subplot_CP_IP(spec_type, subtitle, index, ylab, percent_limits, df_PML, df_PML_CP ,bands)
    
    subtitle  = 'HEREON: N = ' + str(np.sum(~np.isnan(df_HEREON_CP['400'])))
    index = 2
    _resid_subplot_CP_IP(spec_type,subtitle, index, ylab,percent_limits, df_HEREON, df_HEREON_CP ,bands)
    
    subtitle  = 'TARTU: N = ' + str(np.sum(~np.isnan(df_TARTU_CP['400'])))
    index = 3
    _resid_subplot_CP_IP(spec_type,subtitle, index, ylab, percent_limits, df_TARTU, df_TARTU_CP ,bands)

    subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA_CP['400'])))
    index = 4
    _resid_subplot_CP_IP(spec_type, subtitle, index, ylab, percent_limits, df_NASA, df_NASA_CP ,bands)
    
    plt.tight_layout(pad=1.6)     
    
    filename  =  path_output +  '/' + spec_type + '_resiudalsplot_CP_vs_IP.png'
    plt.savefig(filename)
    
    
    return


def plot_unc_CP(spec_type, df_PML, df_NASA, df_TARTU, df_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3'):
    ''' Funtion to plot spectral dependence of % residuals following Tilstone 2020'''  
    
    # QC filtering 
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
    
    
    # spectral reiduals plot
    fig= plt.figure(figsize=(18,6))
    plt.rc('font',size=16)  
    if spec_type == 'Ed':
        ylab = '$E_{d}$ uncertainty [%]'
        plt.suptitle('Uncertainty distribution (post QC) for downwelling irradiance: $E_{d}$(0$^{+}$,$\lambda$)')
        percent_limits = 25
    if spec_type == 'Lsky':
         ylab = '$L_{sky}$ uncertainty [%]'
         plt.suptitle('Uncertainty distribution (post QC) for sky radiance: $L_{sky}$(0$^{+}$,$\lambda$)')
         percent_limits = 25
    if spec_type == 'Lt':
         ylab = '$L_{t}$ uncertainty [%]'
         plt.suptitle('Uncertainty distribution (post QC) upwelling radiance: $L_{t}$(0$^{+}$,$\lambda$)')
         percent_limits = 25
    if spec_type == 'Rrs':
          ylab = '$R_{rs}$ uncertainty [%]'
          plt.suptitle('Uncertainty distribution (post QC) for remote-sensing reflectance: $R_{rs}$($\lambda$)')
          percent_limits = 25
    if spec_type == 'nLw':
          ylab = '$L_{wn}$ uncertainty [%]'
          percent_limits = 25
          plt.suptitle('Uncertainty distribution (post QC) for normalized water-leaving radiance: $L_{wn}$($\lambda$)')

    xlab = 'Wavelength [nm]'
    fig.supxlabel(xlab)
    fig.supylabel(ylab)
        
    subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML['560'])))
    index = 1
    _unc_subplot_CP(spec_type, subtitle, index, ylab, percent_limits, df_PML,bands)
    
    subtitle  = 'HEREON: N = ' + str(np.sum(~np.isnan(df_HEREON['560'])))
    index = 2
    _unc_subplot_CP(spec_type,subtitle, index, ylab,percent_limits, df_HEREON,bands)
    
    subtitle  = 'TARTU: N = ' + str(np.sum(~np.isnan(df_TARTU['560'])))
    index = 3
    _unc_subplot_CP(spec_type,subtitle, index, ylab, percent_limits, df_TARTU,bands)

    subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['560'])))
    index = 4
    _unc_subplot_CP(spec_type, subtitle, index, ylab, percent_limits, df_NASA, bands)
        
    plt.tight_layout(pad=1.6)
    
    filename  =  path_output +  '/' + spec_type + '_uncplot_CP.png'
    plt.savefig(filename)
    
    return


def tabular_summary(spec_type, df_R, df_PML, df_NASA, df_TARTU, df_HEREON, df_RBINS, bands, path_output, Q_mask, Qtype = 'QC_AOC_3'):
    ''' Funtion to output tabular summary of results based on Tilstone 2020'''    
    
    # QC filtering
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
    df_RBINS = df_RBINS[Q_mask[Qtype]==1]
   # df_CNR = df_CNR[Q_mask[Qtype]==1]
   # df_NOAA = df_NOAA[Q_mask[Qtype]==1]

    df_R = df_R[Q_mask[Qtype]==1]
        
    #  columns for table
    Institution = ['PML', 'NASA', 'TARTU', 'HEREON', 'RBINS'] # 'CNR', 'NOAA']
    
    Sensor = ['Seabird-HyperSAS', 'Seabird-HyperSAS', 'TriOS-RAMSES', 'TriOS-RAMSES', '' ]# '', '']
    
    N_PML = np.sum(~np.isnan(df_PML['400']))
    N_NASA = np.sum(~np.isnan(df_NASA['400']))
    N_TARTU = np.sum(~np.isnan(df_TARTU['400']))
    N_HEREON = np.sum(~np.isnan(df_HEREON['400']))
    N_RBINS = np.sum(~np.isnan(df_RBINS['400']))
    # N_CNR = np.sum(~np.isnan(df_CNR['400']))
    # N_NOAA = np.sum(~np.isnan(df_NOAA['400']))
    N_meas = [N_PML, N_NASA, N_TARTU, N_HEREON, N_RBINS] # N_CNR, N_NOAA
    
    RMSD_442_PML = np.sqrt(np.nanmean((df_PML['442.5'] - df_R['442.5'])**2))
    RMSD_442_NASA = np.sqrt(np.nanmean((df_NASA['442.5'] - df_R['442.5'])**2))
    RMSD_442_TARTU =np.sqrt(np.nanmean((df_TARTU['442.5'] - df_R['442.5'])**2))
    RMSD_442_HEREON = np.sqrt(np.nanmean((df_HEREON['442.5'] - df_R['442.5'])**2))
    RMSD_442_RBINS = np.sqrt(np.nanmean((df_RBINS['442.5'] - df_R['442.5'])**2))
    RMSD_442 = [RMSD_442_PML, RMSD_442_NASA, RMSD_442_TARTU, RMSD_442_HEREON, RMSD_442_RBINS]
    
    RMSD_560_PML =  np.sqrt(np.nanmean((df_PML['560'] - df_R['560'])**2))
    RMSD_560_NASA =  np.sqrt(np.nanmean((df_NASA['560'] - df_R['560'])**2))
    RMSD_560_TARTU =  np.sqrt(np.nanmean((df_TARTU['560'] - df_R['560'])**2))
    RMSD_560_HEREON = np.sqrt(np.nanmean((df_HEREON['560'] - df_R['560'])**2))
    RMSD_560_RBINS = np.sqrt(np.nanmean((df_RBINS['560'] - df_R['560'])**2))
    RMSD_560 = [RMSD_560_PML, RMSD_560_NASA, RMSD_560_TARTU, RMSD_560_HEREON, RMSD_560_RBINS]
    
    RMSD_665_PML =  np.sqrt(np.nanmean((df_PML['665'] - df_R['665'])**2))
    RMSD_665_NASA =  np.sqrt(np.nanmean((df_NASA['665'] - df_R['665'])**2))
    RMSD_665_TARTU =  np.sqrt(np.nanmean((df_TARTU['665'] - df_R['665'])**2))
    RMSD_665_HEREON = np.sqrt(np.nanmean((df_HEREON['665'] - df_R['665'])**2))
    RMSD_665_RBINS = np.sqrt(np.nanmean((df_RBINS['665'] - df_R['665'])**2))
    RMSD_665 = [RMSD_665_PML, RMSD_665_NASA, RMSD_665_TARTU, RMSD_665_HEREON, RMSD_665_RBINS]
    
    RPD_442_PML =  100*np.nanmean((df_PML['442.5'] - df_R['442.5'])/df_R['442.5'])
    RPD_442_NASA = 100*np.nanmean((df_NASA['442.5'] - df_R['442.5'])/df_R['442.5'])
    RPD_442_TARTU =  100*np.nanmean((df_TARTU['442.5'] - df_R['442.5'])/df_R['442.5'])
    RPD_442_HEREON = 100*np.nanmean((df_HEREON['442.5'] - df_R['442.5'])/df_R['442.5'])
    RPD_442_RBINS = 100*np.nanmean((df_RBINS['442.5'] - df_R['442.5'])/df_R['442.5'])
    RPD_442 = [RPD_442_PML, RPD_442_NASA, RPD_442_TARTU, RPD_442_HEREON, RPD_442_RBINS]
    
    RPD_560_PML = 100*np.nanmean((df_PML['560'] - df_R['560'])/df_R['560'])
    RPD_560_NASA =  100*np.nanmean((df_NASA['560'] - df_R['560'])/df_R['560'])
    RPD_560_TARTU =  100*np.nanmean((df_TARTU['560'] - df_R['560'])/df_R['560'])
    RPD_560_HEREON = 100*np.nanmean((df_HEREON['560'] - df_R['560'])/df_R['560'])
    RPD_560_RBINS =  100*np.nanmean((df_RBINS['560'] - df_R['560'])/df_R['560'])
    RPD_560 = [RPD_560_PML, RPD_560_NASA, RPD_560_TARTU, RPD_560_HEREON, RPD_560_RBINS]
    
    
    RPD_665_PML =  100*np.nanmean((df_PML['665'] - df_R['665'])/df_R['665'])
    RPD_665_NASA =  100*np.nanmean((df_NASA['665'] - df_R['665'])/df_R['665'])
    RPD_665_TARTU =  100*np.nanmean((df_TARTU['665'] - df_R['665'])/df_R['665'])
    RPD_665_HEREON =  100*np.nanmean((df_HEREON['665'] - df_R['665'])/df_R['665'])
    RPD_665_RBINS =  100*np.nanmean((df_RBINS['665'] - df_R['665'])/df_R['665'])
    RPD_665 = [RPD_665_PML, RPD_665_NASA, RPD_665_TARTU, RPD_665_HEREON, RPD_665_RBINS]
    
    # convert to df format - 
    summary = pd.DataFrame() 
    summary['Institution'] = Institution 
    summary['Sensor type'] = Sensor
    summary['N'] = N_meas
    
    summary['RMSD 442.5'] = RMSD_442
    summary['RPD 442.5'] = RPD_442
    summary['RMSD 560'] = RMSD_560
    summary['RPD 560'] = RPD_560
    summary['RMSD 665'] = RMSD_665
    summary['RPD 665'] = RPD_665
      
    if spec_type != 'Rrs':
        summary['RMSD 442.5'] = summary['RMSD 442.5'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
        summary['RPD 442.5'] = summary['RPD 442.5'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
        summary['RMSD 560'] = summary['RMSD 560'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
        summary['RPD 560'] = summary['RPD 560'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
        summary['RMSD 665'] = summary['RMSD 665'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
        summary['RPD 665'] = summary['RPD 665'].round(2).apply(lambda x: '{0:g}'.format(float(x)))     
    elif spec_type == 'Rrs': 
       summary['RPD 442.5'] = summary['RPD 442.5'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
       summary['RPD 560'] = summary['RPD 560'].round(2).apply(lambda x: '{0:g}'.format(float(x)))
       summary['RPD 665'] = summary['RPD 665'].round(2).apply(lambda x: '{0:g}'.format(float(x)))

    filename  =  path_output +  '/' + spec_type + '_summary.csv'
    summary.to_csv(filename, na_rep ='NaN', index = False)
    filename2  =  path_output +  '/' + spec_type + '_summary.png' 

    dfi.export(summary.style.hide(axis='index'), filename2)
    
    return summary

