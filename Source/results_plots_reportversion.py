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
import matplotlib.transforms as mtransforms

def _resid_subplot_CP_IP(spec_type,system, panel, plot_index, ylab, percent_limits, df_sys, df_R ,bands):
    ''' suplot routine for residuals'''  
    colors = cm.rainbow(np.linspace(0,1,301)) 
    
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
    plt.subplot(3,4,plot_index) 
    # plt.title(system)
   # plt.plot(np.arange(0,12,1), np.zeros(12),linestyle ='dashed',color='gray')
    bp = plt.boxplot(resid ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90]) 
    plt.xlim(0.5,10.5)
    plt.ylim(-percent_limits, percent_limits)
    for i in range(10):
        c_index = int(round(float(bands[i]))) - 400
        bp['boxes'][i].set_facecolor(colors[c_index])
   
   
    bands_round = []
    for i in range(len(bands)): # round bands to nearest nm
        bands_round.append(str(round(float(bands[i]))))
   
   
    plt.xticks([1,2,3,4,5,6,7,8,9,10], bands_round[0:10])
    plt.xticks(rotation = 45)
    
    
    plt.grid(axis='y') 
    
    if plot_index==9: #or plot_index== 5:
       plt.ylabel(ylab)
    #if plot_index > 3:
     #   plt.xlabel('Wavelength [nm]')
  
    ax=plt.gca()
    plt.text(.05, .95,  panel, ha='left', va='top', transform=ax.transAxes,fontsize=20)  
  
    return


def _resid_subplot_CP(spec_type, system, panel, plot_index, ylab, percent_limits, df_sys, df_sys_unc ,df_R ,bands):
    ''' suplot routine for residuals'''  
    colors = cm.rainbow(np.linspace(0,1,301)) 
    
    
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
    plt.subplot(3,4,plot_index) 
    plt.title(system)
    #plt.plot(np.arange(0,12,1), np.zeros(12),linestyle ='dashed',color='gray')
    plt.plot(np.arange(1,11,1), df_sys_unc[0:10],linestyle ='dashed',color='gray',label = 'Median uncertainty')
    plt.plot(np.arange(1,11,1), -df_sys_unc[0:10],linestyle ='dashed',color='gray')
    if plot_index == 1:  # tunred off for nlw
        plt.legend(loc='lower left', fontsize=15)
    bp = plt.boxplot(resid ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90]) 
    plt.xlim(0.5,10.5)
    plt.ylim(-percent_limits, percent_limits)
    for i in range(10):
        c_index = int(round(float(bands[i]))) - 400
        bp['boxes'][i].set_facecolor(colors[c_index])
   
        
   
    plt.xticks([1,2,3,4,5,6,7,8,9,10],bands[0:10])
    plt.xticks([])
    plt.xticks(rotation = 45)
    
    plt.grid(axis='y') 
        
    if plot_index==1: #or plot_index== 5:
       plt.ylabel(ylab)
       
    #if plot_index==1 or plot_index== 5:
     #   plt.ylabel(ylab)
    #if plot_index > 3:
    #   plt.xlabel('Wavelength [nm]')
 
    # Set both x- and y-axis limits to [0, 10] instead of default [0, 1]
    #  ax=plt.gca()
    ##  ax.text('left', 'top', r'A.', fontsize=18)
  
    ax=plt.gca()
    plt.text(.05, .95,  panel, ha='left', va='top', transform=ax.transAxes,fontsize=20)  
    
    return

def _unc_subplot_CP(spec_type,system, panel, plot_index, ylab, percent_limits, df_sys ,bands):
    ''' suplot routine for residuals'''  
    colors = cm.rainbow(np.linspace(0,1,301)) 
    
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

    plt.subplot(3,4,plot_index) 
    # plt.title(system)
    #  plt.plot(np.arange(0,12,1), np.zeros(12),linestyle ='dashed',color='gray')
    bp = plt.boxplot(resid ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90]) 
    plt.xlim(0.5,10.5)
    plt.ylim(0, percent_limits)
    for i in range(10):
        c_index = int(round(float(bands[i]))) - 400
        bp['boxes'][i].set_facecolor(colors[c_index])
   
    
    bands_round = []
    for i in range(len(bands)): # round bands to nearest nm
        bands_round.append(str(round(float(bands[i]))))
   
    
    plt.xticks([1,2,3,4,5,6,7,8,9,10], bands_round[0:10])
    plt.xticks([])
    plt.xticks(rotation = 45) 
    plt.grid(axis='y') 
    
    if plot_index==5: #or plot_index== 5:
       plt.ylabel(ylab)
       
    # if plot_index==1 or plot_index== 5:
    #   plt.ylabel(ylab)
    # if plot_index > 3:
    #   plt.xlabel('Wavelength [nm]')
    ax=plt.gca()
    plt.text(.05, .95,  panel, ha='left', va='top', transform=ax.transAxes,fontsize=20)  
   
    return



def _resid_subplot_nLw(spec_type, system, panel, plot_index, ylab, percent_limits, df_sys ,df_R,bands):
    ''' suplot routine for residuals'''  
    colors = cm.rainbow(np.linspace(0,1,301)) 
        
    
    if spec_type == 'nLw':
        resid = [] # res
        resid.append([])
        for i in range(1,8,1):# residual disrtibution in each band
            resid_i = 100*np.array((df_sys[str(bands[i])] -  df_R[str(bands[i])])/df_R[str(bands[i])])
            resid_i = resid_i[~np.isnan(resid_i)]
            resid.append(resid_i)
        resid.append([])
        resid.append([])

    print(resid)
    plt.subplot(2,5,plot_index) 
    plt.title(system)


    bp = plt.boxplot(resid ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90]) 
    plt.xlim(1.5,8.5)
    plt.ylim(-percent_limits,percent_limits)
    for i in range(10):
        c_index = int(round(float(bands[i]))) - 400
        bp['boxes'][i].set_facecolor(colors[c_index])
    
    bands_round = []
    for i in range(len(bands)): # round bands to nearest nm
        bands_round.append(str(round(float(bands[i]))))
   
    plt.xticks([2,3,4,5,6,7,8],bands_round[1:8])
    #plt.xticks([])
    plt.xticks(rotation = 45) 
    plt.grid(axis='y') 
    
        
    if plot_index==1 or plot_index== 6:
       plt.ylabel(ylab)
            
    if plot_index < 6:
        plt.xticks([])
    
       
    #if plot_index==1 or plot_index== 5:
     #   plt.ylabel(ylab)
    #if plot_index > 3:
    #   plt.xlabel('Wavelength [nm]')
 
    # Set both x- and y-axis limits to [0, 10] instead of default [0, 1]
    #  ax=plt.gca()
    ##  ax.text('left', 'top', r'A.', fontsize=18)
  
    ax=plt.gca()
    plt.text(.05, .95,  panel, ha='left', va='top', transform=ax.transAxes,fontsize=20)  
    
    return


def _resid_subplot_nLw_SP_HP(spec_type, system, panel, plot_index, ylab, percent_limits, df_sys ,df_R,bands):
    ''' suplot routine for residuals'''  
    
    colors = cm.rainbow(np.linspace(0,1,301)) 
        
    if spec_type == 'nLw':
        resid = [] # res
        resid.append([])
        for i in range(1,8,1):# residual disrtibution in each band
            resid_i = 100*np.array((df_sys[str(bands[i])] -  df_R[str(bands[i])])/df_R[str(bands[i])])
            resid_i = resid_i[~np.isnan(resid_i)]
            resid.append(resid_i)
        resid.append([])
        resid.append([])

    print(resid)
    plt.subplot(1, 2, plot_index) 
    plt.title(system)


    bp = plt.boxplot(resid ,showfliers=True,patch_artist=True, medianprops=dict(color='black'), whis=[10,90]) 
    plt.xlim(1.5,8.5)
    plt.ylim(-percent_limits,percent_limits)
    for i in range(10):
        c_index = int(round(float(bands[i]))) - 400
        bp['boxes'][i].set_facecolor(colors[c_index])
    
    bands_round = []
    for i in range(len(bands)): # round bands to nearest nm
        bands_round.append(str(round(float(bands[i]))))
   
    plt.xticks([2,3,4,5,6,7,8],bands_round[1:8])
    #plt.xticks([])
    plt.xticks(rotation = 45) 
    plt.grid(axis='y') 
    
        
    plt.ylabel(ylab)
    plt.xlabel('Wavelength [nm]')

    
       
    #if plot_index==1 or plot_index== 5:
     #   plt.ylabel(ylab)
    #if plot_index > 3:
    #   plt.xlabel('Wavelength [nm]')
 
    # Set both x- and y-axis limits to [0, 10] instead of default [0, 1]
    #  ax=plt.gca()
    ##  ax.text('left', 'top', r'A.', fontsize=18)
  
    ax=plt.gca()
    plt.text(.05, .95,  panel, ha='left', va='top', transform=ax.transAxes,fontsize=20)  
    
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


def residuals_combined(spec_type, df_R_CP, df_PML_CP, df_NASA_CP, df_TARTU_CP, df_HEREON_CP, df_PML, df_NASA, df_TARTU, df_HEREON, df_unc_PML, df_unc_NASA, df_unc_TARTU, df_unc_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3'):
    ''' Combined plot for 4 reference systems processed with HCP:
        
        Row 1: HCP 4-way baseline experiment (median uncertainties in background)
        Row 2: Uncertainty distrbiutions
        Row 3: HCP versus IP'''  
    
    # QC filtering 
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
    
    df_PML_CP = df_PML_CP[Q_mask[Qtype]==1]
    df_NASA_CP = df_NASA_CP[Q_mask[Qtype]==1]
    df_TARTU_CP = df_TARTU_CP[Q_mask[Qtype]==1]
    df_HEREON_CP = df_HEREON_CP[Q_mask[Qtype]==1]
    
    df_unc_PML = df_unc_PML[Q_mask[Qtype]==1]
    df_unc_NASA =  df_unc_NASA[Q_mask[Qtype]==1]
    df_unc_TARTU = df_unc_TARTU[Q_mask[Qtype]==1]
    df_unc_HEREON = df_unc_HEREON[Q_mask[Qtype]==1]
    
    unc_med_PML = np.nanmedian(df_unc_PML.iloc[:,1:-1],0)
    unc_med_NASA = np.nanmean(df_unc_NASA.iloc[:,1:-1],0)
    unc_med_TARTU = np.nanmedian(df_unc_TARTU.iloc[:,1:-1],0)
    unc_med_HEREON = np.nanmedian(df_unc_HEREON.iloc[:,1:-1],0)
    
    df_R_CP = df_R_CP[Q_mask[Qtype]==1]
    
    
    #row 1 spectral reiduals plot
    fig= plt.figure(figsize=(17,14))
    
    
    
    plt.rc('font',size=18)  
    ylab = '100($X_{HCP}$ - $X_{R}$)/$X_{R}$  [%]'
    if spec_type == 'Ed':
        percent_limits = 8
    if spec_type == 'Lsky':
         percent_limits = 10
    if spec_type == 'Lt':
         percent_limits = 10
    if spec_type == 'Rrs':
        percent_limits = 15
    if spec_type == 'nLw':
        percent_limits = 16
        
    subtitle  = 'HyperSAS: N = ' + str(np.sum(~np.isnan(df_PML['400'])))
    index = 1
    _resid_subplot_CP(spec_type, subtitle, 'A.', index, ylab, percent_limits, df_PML, unc_med_PML, df_R_CP, bands)
  
    subtitle  = 'RAMSES-A: N = ' + str(np.sum(~np.isnan(df_TARTU['400'])))
    index = 2
    _resid_subplot_CP(spec_type,subtitle, 'B.', index, ylab, percent_limits, df_TARTU, unc_med_TARTU, df_R_CP, bands)

    subtitle  = 'RAMSES-B: N = ' + str(np.sum(~np.isnan(df_HEREON['400'])))
    index = 3
    _resid_subplot_CP(spec_type,subtitle, 'C.', index, ylab,percent_limits, df_HEREON, unc_med_HEREON, df_R_CP, bands)

    subtitle  = 'pySAS: N = ' + str(np.sum(~np.isnan(df_NASA['400'])))
    index = 4
    _resid_subplot_CP(spec_type, subtitle,'D.', index, ylab, percent_limits, df_NASA, unc_med_NASA, df_R_CP, bands)
    
 
    if spec_type == 'Ed':
        percent_limits = 15
    if spec_type == 'Lsky':
         percent_limits = 15
    if spec_type == 'Lt':
         percent_limits = 25
    if spec_type == 'Rrs':
         percent_limits = 30

    xlab = 'Wavelength [nm]'
    fig.supxlabel(xlab)
    ylab = 'Uncertainty   [%]'
    
    #fig.supylabel(ylab)
        
  #  subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML['560'])))
    index = 5
    _unc_subplot_CP(spec_type, subtitle,'E.',  index, ylab, percent_limits, df_unc_PML,bands)
    
   # subtitle  = 'TARTU: N = ' + str(np.sum(~np.isnan(df_TARTU['560'])))
    index = 6
    _unc_subplot_CP(spec_type,subtitle, 'F.', index, ylab, percent_limits, df_unc_TARTU,bands)

   # subtitle  = 'HEREON: N = ' + str(np.sum(~np.isnan(df_HEREON['560'])))
    index = 7
    _unc_subplot_CP(spec_type,subtitle,'G.', index, ylab,percent_limits, df_unc_HEREON,bands)
    
  #  subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA['560'])))
    index = 8
    _unc_subplot_CP(spec_type, subtitle,'H.', index, ylab, percent_limits, df_unc_NASA, bands)
    
    if spec_type == 'Ed':
        percent_limits = 2
    if spec_type == 'Lsky':
        percent_limits = 2
    if spec_type == 'Lt':
        percent_limits = 10
    if spec_type == 'Rrs':
        percent_limits = 15
    if spec_type == 'nLw':
        percent_limits = 40

    xlab = 'Wavelength [nm]'
    ylab = '100($X_{IP}$ - $X_{HCP}$)/$X_{HCP}$   [%]'
    fig.supxlabel(xlab)
   # fig.supylabel(ylab)
        
  #  subtitle  = 'PML: N = ' + str(np.sum(~np.isnan(df_PML_CP['400'])))
    index = 9
    _resid_subplot_CP_IP(spec_type, subtitle,'I.', index, ylab, percent_limits, df_PML, df_PML_CP,bands)


    #subtitle  = 'TARTU: N = ' + str(np.sum(~np.isnan(df_TARTU_CP['400'])))
    index = 10
    _resid_subplot_CP_IP(spec_type,subtitle, 'J.',index, ylab, percent_limits, df_TARTU, df_TARTU_CP ,bands)

   # subtitle  = 'HEREON: N = ' + str(np.sum(~np.isnan(df_HEREON_CP['400'])))
    index = 11
    _resid_subplot_CP_IP(spec_type,subtitle, 'K.', index, ylab,percent_limits, df_HEREON, df_HEREON_CP ,bands)
    
   # subtitle  = 'NASA: N = ' + str(np.sum(~np.isnan(df_NASA_CP['400'])))
    index = 12
    _resid_subplot_CP_IP(spec_type, subtitle,'L.', index, ylab, percent_limits, df_NASA, df_NASA_CP ,bands)
    
    plt.tight_layout(pad=1.2)
    
    
    filename  =  path_output +  '/' + spec_type + '_summary.png' 
    plt.savefig(filename, dpi=300)
    

    return


def residuals_combined_nlw(spec_type, df_SEAPRISM, df_NOAA, df_PML, df_NASA, df_NASA2, df_TARTU, df_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3'):
    
    ''' Combined plot for nlw
    #      Row 1: Seaprism baseline
    #      Row 2: Hyperpro'''
    
    # QC filtering 
    df_PML = df_PML[Q_mask[Qtype]==1]
    df_NASA = df_NASA[Q_mask[Qtype]==1]
    df_NASA2 = df_NASA2[Q_mask[Qtype]==1]
    df_TARTU = df_TARTU[Q_mask[Qtype]==1]
    df_HEREON = df_HEREON[Q_mask[Qtype]==1]
  
    df_SEAPRISM = df_SEAPRISM[Q_mask[Qtype]==1]
    df_NOAA =  df_NOAA[Q_mask[Qtype]==1]
    
    df_PML_NOAA_match = df_PML[~np.isnan(df_NOAA['560'])] # filters for NaNs in HP dayta
    df_NASA_NOAA_match =  df_NASA[~np.isnan(df_NOAA['560'])] 
    df_NASA2_NOAA_match = df_NASA2[~np.isnan(df_NOAA['560'])] 
    df_TARTU_NOAA_match = df_TARTU[~np.isnan(df_NOAA['560'])] 
    df_HEREON_NOAA_match = df_HEREON[~np.isnan(df_NOAA['560'])] 
    
    #row 1 spectral reiduals plot
    fig= plt.figure(figsize=(17,12))
    
    plt.rc('font',size=17)  
  
    percent_limits  = 30
    xlab = 'Wavelength [nm]'
    ylab = '100($X_{IP}$ - $X_{R}$)/$X_{R}$: SeaPRISM  [%]'
        
    subtitle  = 'HyperSAS: N = ' + str(np.sum(~np.isnan(df_PML['560'])))
    index = 1
    _resid_subplot_nLw(spec_type, subtitle, 'A.', index, ylab, percent_limits, df_PML, df_SEAPRISM, bands)
  
    subtitle  = 'RAMSES-A: N = ' + str(np.sum(~np.isnan(df_TARTU['560'])))
    index = 2
    _resid_subplot_nLw(spec_type,subtitle, 'B.', index, ylab, percent_limits, df_TARTU, df_SEAPRISM, bands)

    subtitle  = 'RAMSES-B: N = ' + str(np.sum(~np.isnan(df_HEREON['560'])))
    index = 3
    _resid_subplot_nLw(spec_type,subtitle, 'C.', index, ylab,percent_limits, df_HEREON, df_SEAPRISM, bands)

    subtitle  = 'pySAS: N = ' + str(np.sum(~np.isnan(df_NASA['560'])))
    index = 4
    _resid_subplot_nLw(spec_type, subtitle,'D.', index, ylab, percent_limits, df_NASA, df_SEAPRISM, bands)
    
    subtitle  = 'pySAS-v2: N = ' + str(np.sum(~np.isnan(df_NASA2['560'])))
    index = 5
    _resid_subplot_nLw(spec_type, subtitle,'E.', index, ylab, percent_limits, df_NASA2, df_SEAPRISM, bands)
    
    percent_limits  = 15
    xlab = 'Wavelength [nm]'
    ylab = '100($X_{IP}$ - $X_{R}$)/$X_{R}$: HyperPro II [%]'

    fig.supxlabel(xlab)
    
    subtitle  = 'HyperSAS: N = ' + str(np.sum(~np.isnan(df_PML_NOAA_match['560'])))
    index = 6
    _resid_subplot_nLw(spec_type, subtitle, 'F.', index, ylab, percent_limits, df_PML, df_NOAA, bands)
  
    subtitle  =  'RAMSES-A: N = ' + str(np.sum(~np.isnan(df_TARTU_NOAA_match['560'])))
    index = 7
    _resid_subplot_nLw(spec_type,subtitle, 'G.', index, ylab, percent_limits, df_TARTU, df_NOAA, bands)

    subtitle  = 'RAMSES-B: N = ' + str(np.sum(~np.isnan(df_HEREON_NOAA_match['560'])))
    index = 8
    _resid_subplot_nLw(spec_type,subtitle, 'H.', index, ylab,percent_limits, df_HEREON, df_NOAA, bands)

    subtitle  = 'pySAS: N = ' + str(np.sum(~np.isnan(df_NASA_NOAA_match['560'])))
    index = 9
    _resid_subplot_nLw(spec_type, subtitle,'I.', index, ylab, percent_limits, df_NASA, df_NOAA, bands)
    
    subtitle  = 'pySAS: N = ' + str(np.sum(~np.isnan(df_NASA_NOAA_match['560'])))
    index = 10
    _resid_subplot_nLw(spec_type, subtitle,'J.', index, ylab, percent_limits, df_NASA2, df_NOAA, bands)

    

    plt.tight_layout(pad=1.2)
    
    
    filename  =  path_output +  '/' + spec_type + '_summaryzoom.png' 
    plt.savefig(filename, dpi=300)
    

    return



def SB_VS_HP(spec_type, df_SEAPRISM, df_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3'):
    
    ''' Combined plot for nlw
    #      Row 1: Seaprism baseline
    #      Row 2: Hyperpro'''
    
    # QC filtering 

    df_SEAPRISM = df_SEAPRISM[Q_mask[Qtype]==1]
    df_NOAA =  df_NOAA[Q_mask[Qtype]==1]
    
    #row 1 spectral reiduals plot
    fig= plt.figure(figsize=(12,8))
    
    plt.rc('font',size=17)  
  

    plt.subplot(1,2,1)    
    X = np.arange(0,30,1)
    plt.plot(X,X,color='gray')
    colors = cm.rainbow(np.linspace(0,1,301)) 
    bands_round = []
    for i in range(len(bands)): # round bands to nearest nm
        bands_round.append(str(round(float(bands[i]))))
     
    for i in range(1,8,1):
            c_index = int(round(float(bands[i]))) - 400
            plt.scatter(df_SEAPRISM[str(bands[i])], df_NOAA[str(bands[i])], color=colors[c_index], facecolors='none',label=str(bands_round[i]) + ' nm')  


    plt.gca().set_aspect('equal')
    plt.title('N = ' + str(np.sum(~np.isnan(df_NOAA['560']))))    
    plt.legend(loc=4)
    xlab = '$L_{wn}$: SeaPrism [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
    ylab = '$L_{wn}$: HyperPro [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
      
    ax=plt.gca()
    plt.text(.05, .95,  'A.', ha='left', va='top', transform=ax.transAxes,fontsize=20)  
    
    
    percent_limits  = 30
    xlab = 'Wavelength [nm]'
    ylab = '100($X_{HyperPro}$ - $X_{SeaPrism}$)/$X_{SeaPrism}$  [%]'    
    index = 2
    _resid_subplot_nLw_SP_HP(spec_type, '',  'B.', index, ylab, percent_limits, df_NOAA, df_SEAPRISM, bands)
    plt.tight_layout(pad=1.2)
    

    filename  =  path_output +  '/' + spec_type + 'SP_HP_summary.png' 
    plt.savefig(filename, dpi=300)
    
    return


