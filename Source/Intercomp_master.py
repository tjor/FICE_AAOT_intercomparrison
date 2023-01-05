#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 15:36:24 2022

@author: tjor

This script performs the system intercomparrision for the FICE AAOT deploymnet

    
# 1. Baselines

 Ed, Lt, Ls, Rrs: use 4-way mean of NASA, PML, Tartu, HEREON
 nLw: use SEAPRISM-Aeronet as independent reference (Zibordi et al. 2018) or NOAAA

# 2. QC levels
 (i) Zibordi's Aeronet/Seaprism QC 
 Cloudiness filtering (subset passed)

# 3. Methodology follows Tilstone 2020.

    Scatter plots and spectral depedence of residuals as a box plot 
    (show all bands simulatneously)

    Match-up stats - bands 443, 560 and 665 nm
    N - no. of stations used for 
    RPD - relative % difference
    RMSD - root mean square diviation
        
#

"""

import numpy as np
import pandas as pd
import datetime 
import os

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.dates as mdates

import dataframe_image as dfi

import read_IP_data as rd_IP
import read_CP_data as rd_CP
import QC_and_baselines as QC
import results_plots as rp

# this should really be within plot functions - its easiest just to run with a blank arguement in the main
def time_series():
    
    time_start_PML = np.empty(len(Ed_PML), dtype=object)
    time_start_NASA = np.empty(len(Ed_PML), dtype=object)
    time_start_HEREON = np.empty(len(Ed_PML), dtype=object)
    time_start_CNR= np.empty(len(Ed_PML), dtype=object)
    time_start_RBINS = np.empty(len(Ed_PML), dtype=object)
    time_start_TARTU = np.empty(len(Ed_PML), dtype=object)
    
    for i in range(len(Ed_PML)):
        time_start_PML[i] = datetime.datetime.strptime(Ed_PML['time_start'][i],'%Y-%m-%d %H:%M:%S')
        if Ed_NASA['time_start'][i]  != None:
           time_start_NASA[i] = datetime.datetime.strptime(Ed_NASA['time_start'][i],'%Y-%m-%d %H:%M:%S')   
        if str(Ed_RBINS['time_start'][i])[0:4] == '2022':
            time_start_RBINS[i] = datetime.datetime.strptime(Ed_RBINS['time_start'][i],'%Y-%m-%d %H:%M:%S')   
        if Ed_TARTU['time_start'][i] != None:
            time_start_TARTU[i] = datetime.datetime.strptime(Ed_TARTU['time_start'][i],'%Y-%m-%d %H:%M:%S')  
        if Ed_CNR['time_start'][i] != None:
           time_start_CNR[i] = datetime.datetime.strptime(Ed_CNR['time_start'][i],'%Y-%m-%d %H:%M:%S')  
        #if Ed_HEREON['time_start'][i] != None:
         #  time_start_HEREON[i] = datetime.datetime.strptime(Ed_HEREON['time_start'][i],'%Y-%m-%d %H:%M:%S')  
              
    plt.figure(figsize=(18,12))          
    plt.suptitle('Time series from 19/07/2022: no QC')
    plt.subplot(2,3,1)
    plt.title('Ed(560)')
    plt.plot_date(time_start_PML[32:56], Ed_PML['560'][32:56], ms=3 , label='PML',color='orange')
    plt.plot_date(time_start_TARTU[32:56], Ed_TARTU['560'][32:56], ms=3, label='TARTU',color='red')
    plt.plot_date(time_start_PML[32:56], Ed_HEREON['560'][32:56], ms=3, label='HERON',color='pink')
    plt.plot_date(time_start_NASA[32:56], Ed_NASA['560'][32:56], ms=3, label='NASA',color='blue')
    plt.plot_date(time_start_RBINS[32:56], Ed_RBINS['560'][32:56], ms=3, label='RBINS',color='black')
    plt.plot_date(time_start_CNR[32:56], Ed_CNR['560'][32:56], ms=3, label='CNR',color='green')
    plt.legend()
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    plt.xlabel('UTC time [hrs]')
    plt.ylabel('[mW m$^{-2}$ nm$^{-1}$]')
    
    plt.subplot(2,3,2)
    plt.title('Lsky(560)')
    plt.plot_date(time_start_PML[32:56], Lsky_PML['560'][32:56], ms=3 , label='PML',color='orange')
    plt.plot_date(time_start_TARTU[32:56], Lsky_TARTU['560'][32:56], ms=3, label='TARTU',color='red')
    plt.plot_date(time_start_PML[32:56], Lsky_HEREON['560'][32:56], ms=3, label='HERON',color='pink')
    plt.plot_date(time_start_NASA[32:56], Lsky_NASA['560'][32:56], ms=3, label='NASA',color='blue')
    plt.plot_date(time_start_RBINS[32:56], Lsky_RBINS['560'][32:56], ms=3, label='RBINS',color='black')
    plt.plot_date(time_start_CNR[32:56], Lsky_CNR['560'][32:56], ms=3, label='CNR',color='green')
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    plt.ylim(20,60)
    plt.xlabel('UTC time [hrs]')
    plt.ylabel('[mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
    
    plt.subplot(2,3,3)
    plt.title('Lt(560)')
    plt.plot_date(time_start_PML[32:56], Lt_PML['560'][32:56], ms=3 , label='PML',color='orange')
    plt.plot_date(time_start_TARTU[32:56], Lt_TARTU['560'][32:56], ms=3, label='TARTU',color='red')
    plt.plot_date(time_start_PML[32:56], Lt_HEREON['560'][32:56], ms=3, label='HERON',color='pink')
    plt.plot_date(time_start_NASA[32:56], Lt_NASA['560'][32:56], ms=3, label='NASA',color='blue')
    plt.plot_date(time_start_RBINS[32:56], Lt_RBINS['560'][32:56], ms=3, label='RBINS',color='black')
    plt.plot_date(time_start_CNR[32:56], Lt_CNR['560'][32:56], ms=3, label='CNR',color='green')
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    plt.ylim(6,18)
    plt.xlabel('UTC time [hrs]')
    plt.ylabel('[mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
    
    plt.subplot(2,3,4)
    plt.title('Rrs(560)')
    plt.plot_date(time_start_PML[32:56], Rrs_PML['560'][32:56], ms=3 , label='PML',color='orange')
    plt.plot_date(time_start_TARTU[32:56], Rrs_TARTU['560'][32:56], ms=3, label='TARTU',color='red')
    plt.plot_date(time_start_PML[32:56], Rrs_HEREON['560'][32:56], ms=3, label='HERON',color='pink')
    plt.plot_date(time_start_NASA[32:56], Rrs_NASA['560'][32:56], ms=3, label='NASA',color='blue')
    plt.plot_date(time_start_RBINS[32:56], Rrs_RBINS['560'][32:56], ms=3, label='RBINS',color='black')
    plt.plot_date(time_start_CNR[32:56], Rrs_CNR['560'][32:56], ms=3, label='CNR',color='green')
    plt.xlabel('UTC time [hrs]')
    plt.ylim(0.004,0.016)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    plt.xlabel('UTC time [hrs]')
    plt.ylabel('[sr$^{-1}$]')
           
    plt.subplot(2,3,5)
    plt.title('Relative azimuth (absolute)')
    plt.plot_date(time_start_PML[32:56], abs(Ed_PML['azimuth'][32:56]), ms=6 , label='PML, TARTU, HERON',color='red')
    plt.plot_date(time_start_NASA[32:56], abs(Ed_NASA['azimuth'][32:56]), ms=4, label='NASA',color='blue')
    plt.plot_date(time_start_RBINS[32:56], abs(Ed_RBINS['azimuth'][32:56]), ms=3, label='RBINS',color='black')
    plt.plot_date(time_start_CNR[32:56], abs(Ed_CNR['azimuth'][32:56]), ms=3, label='CNR',color='green')
    plt.legend()
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    plt.xlabel('UTC time [hrs]')
    plt.ylabel('deg')
    
    plt.subplot(2,3,6)
    plt.title('Windspeed')
    plt.plot_date(time_start_TARTU[32:56], Ed_TARTU['windspeed'][32:56], ms=3 , label='PML',color='black')
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H'))
    plt.xlabel('UTC time [hrs]')
    plt.ylabel('[ms$^{-1}$]')
    plt.tight_layout(pad=1.8)
    
    filename  =  path_output +  '/19-07-2022_timeseries_noQC.png'
    plt.savefig(filename)
 
    return

def station_time_series_CP():
    
    fig = plt.figure(figsize=(18,12))
    plt.rc('font', size=16)     
    plt.suptitle('Station time series for Community Processor (no QC)')
    
    plt.subplot(2,2,1)
    plt.title('$E_{d}$(400)')
    plt.plot(Ed_PML_CP['400'],label='PML')    
    plt.plot(Ed_HEREON_CP['400'],label='HEREON')
    plt.plot(Ed_TARTU_CP['400'], label='TARTU')
    plt.plot(Ed_NASA_CP['400'], label='NASA')
    #plt.xlabel('Station number')
    plt.ylabel('$E_{d}(400)$ [mW m$^{-2}$ nm$^{-1}$]')
    plt.legend()

    plt.subplot(2,2,2)
    plt.title('$L_{sky}$(400)')
    plt.plot(Lsky_PML_CP['400'],label='PML')    
    plt.plot(Lsky_HEREON_CP['400'],label='HEREON')
    plt.plot(Lsky_TARTU_CP['400'],label='TARTU')
    plt.plot(Lsky_NASA_CP['400'], label='NASA')
   # plt.xlabel('Station number')
    plt.ylabel('$L_{sky}(400)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
    
    plt.subplot(2,2,3)
    plt.title('$L_{t}$(400)')
    plt.plot(Lt_PML_CP['400'],label='PML')    
    plt.plot(Lt_HEREON_CP['400'],label='HEREON')
    plt.plot(Lt_TARTU_CP['400'],label='TARTU')
    plt.plot(Lt_NASA_CP['400'], label='NASA')
    plt.xlabel('Station number')
    plt.ylabel('$L_{t}(400)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
       
    plt.subplot(2,2,4)
    plt.title('$R_{rs}$(400)')
    plt.plot(Rrs_PML_CP['400'],label='PML')    
    plt.plot(Rrs_HEREON_CP['400'],label='HEREON')
    plt.plot(Rrs_TARTU_CP['400'],label='TARTU')
    plt.plot(Rrs_NASA_CP['400'], label='NASA')
    plt.xlabel('Station number')
    plt.ylabel('$L_{t}(400)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
    
    filename  =  path_output +  '/'  '_stationtimeseries_CP.png'
    plt.savefig(filename)

    plt.tight_layout(pad=1.8)

    return


def station_time_series_IP():
    
    fig = plt.figure(figsize=(18,12))
    plt.rc('font', size=16)     
    plt.suptitle('Station time series for Individual Processors (no QC)')
    
    plt.subplot(2,2,1)
    plt.title('$E_{d}$(560)')
    plt.plot(Ed_PML['560'],label='PML')    
    plt.plot(Ed_HEREON['560'],label='HEREON')
    plt.plot(Ed_TARTU['560'], label='TARTU')
    plt.plot(Ed_NASA['560'], label='NASA')
    plt.xlabel('Station number')
    plt.ylabel('$E_{d}(560)$ [mW m$^{-2}$ nm$^{-1}$]')
    plt.legend()

    plt.subplot(2,2,2)
    plt.title('$L_{sky}$(560)')
    plt.plot(Lsky_PML['560'],label='PML')    
    plt.plot(Lsky_HEREON['560'],label='HEREON')
    plt.plot(Lsky_TARTU['560'],label='TARTU')
    plt.plot(Lsky_NASA['560'], label='NASA')
    plt.xlabel('Station number')
    plt.ylabel('$L_{sky}(560)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
    
    plt.subplot(2,2,3)
    plt.title('$L_{t}$(560)')
    plt.plot(Lt_PML['560'],label='PML')    
    plt.plot(Lt_HEREON['560'],label='HEREON')
    plt.plot(Lt_TARTU['560'],label='TARTU')
    plt.plot(Lt_NASA['560'], label='NASA')
    plt.xlabel('Station number')
    plt.ylabel('$L_{t}(560)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
       
    plt.subplot(2,2,4)
    plt.title('$R_{rs}$(560)')
    plt.plot(Rrs_PML['560'],label='PML')    
    plt.plot(Rrs_HEREON['560'],label='HEREON')
    plt.plot(Rrs_TARTU['560'],label='TARTU')
    plt.plot(Rrs_NASA['560'], label='NASA')
    plt.xlabel('Station number')
    plt.ylabel('$L_{t}(560)$ [mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$]')
    
    filename  =  path_output +  '/' + '_stationtimeseries_IP.png'
    plt.savefig(filename)

    plt.tight_layout(pad=1.8)

    return

def spectra_CP_IP():
    'plot to show issues with CP'
    
    bands = [400, 412.5, 442.5,	490, 510, 560, 620,	665, 673.75, 681.25, 708.75, 753.75, 761.25, 764.375, 767.5, 778.75, 865, 885, 900]
    
    plt.figure()
    plt.title('Tartu: station 25: Ed spectrum for OLCI bands')
    plt.plot(bands[:-3], np.array(Ed_TARTU.iloc[25,2:18]),label='IP')
    plt.plot(bands[:-3], np.array(Ed_TARTU_CP.iloc[25,1:17]),label='CP',linestyle='dashed')
    plt.xlabel('wavelength [nm]')
    plt.ylabel('Ed [mW m$^{-2}$ nm$^{-1}]$')
    plt.ylim(700,1600)
    plt.legend()
    
    
    plt.figure()
    plt.title('Tartu: station 30: Ed spectrum for OLCI bands')
    plt.plot(bands[:-3], np.array(Ed_TARTU.iloc[30,2:18]),label='IP')
    plt.plot(bands[:-3], np.array(Ed_TARTU_CP.iloc[30,1:17]),label='CP',linestyle='dashed')
    plt.xlabel('wavelength [nm]')
    plt.ylabel('Ed [mW m$^{-2}$ nm$^{-1}]$')
    plt.ylim(600,1400)
    plt.legend()
    
    plt.figure()
    plt.title('PML: station 25: Ed spectrum for OLCI bands')
    plt.plot(bands[:-3], np.array(Ed_PML.iloc[25,3:19]),label='IP')
    plt.plot(bands[:-3], np.array(Ed_PML_CP.iloc[25,1:17]),label='CP',linestyle='dashed')
    plt.xlabel('wavelength [nm]')
    plt.ylabel('Ed [mW m$^{-2}$ nm$^{-1}]$')
    plt.ylim(700,1600)
    plt.legend()
    
    plt.figure()
    plt.title('PML: station 30: Ed spectrum for OLCI bands')
    plt.plot(bands[:-3], np.array(Ed_PML.iloc[30,3:19]),label='IP')
    plt.plot(bands[:-3], np.array(Ed_PML_CP.iloc[30,1:17]),label='CP',linestyle='dashed')
    plt.xlabel('wavelength [nm]')
    plt.ylabel('Ed [mW m$^{-2}$ nm$^{-1}]$')
    plt.ylim(600,1400)
    plt.legend()

    return

 
if __name__ == '__main__':
    
    # path to data + output
    dir_data = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/DataSubmissions'
    path_output = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/Output'

    # team submissions
    path_PML = dir_data + '/PML/FICE_submission_V4/FRM4SOC_2_FICE_22_AAOT_PML_HSAS_stationsummary_V4_rhoM99_nLw_localChl_noNIR_Earthsuncorrected.csv'
    #path_NASA = dir_data  + '/NASA/FRM4SOC_2_FICE_22_NASA_pySAS_Sentinel3A_rev2.csv'    # NASA - no NIR
    path_NASA = dir_data  + '/NASA/FRM4SOC_2_FICE_22_NASA_pySAS_Sentinel3A_rev3.csv'   # NASA - with NIR
    path_TARTU = dir_data  + '/UniTartu/averages.csv'  
    path_RBINS = dir_data  + '/RBINS/FRM4SOC2_FICE_RBINS_2022-11-25_ViewCorrected_EarthSunCorrected.csv'
    path_HEREON = dir_data + '/HEREON/FRM4SOC2_FICE22_HEREON_DATA_OLCI'
    path_NOAA = dir_data  + '/NOAA/NOAA_Hyperpro_sheet2.csv'
    path_CNR = dir_data + '/CNR/FRM4SOC_2_FICE_22_AAOT_CNR_HYPSTAR_w-rel-az.csv'
    
 
    # addtional data (QC + references)    
    path_QC = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/Aeronet_QC_mask/FRM4SOC-AAOT_V3_ope.txt'
    path_NLW = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/nLw_Zibordireference'    

    # OLCI bands
    bands = [str(400), str(412.5), str(442.5),	str(490), str(510), str(560), str(620),	str(665), str(673.75), str(681.25), str(708.75), str(753.75), str(761.25), str(764.375), str(767.5), str(778.75), str(865), str(885), str(900)]
    
    # Read data Individually Processed data. Each team has own file reader in read_IP.data.py to homogenize data
    Ed_PML, Lsky_PML, Lt_PML, Rrs_PML, Rrs_std_PML, nLw_PML = rd_IP.read_PML_data(path_PML, bands)
    Ed_NASA, Lsky_NASA, Lt_NASA, Rrs_NASA, Rrs_std_NASA, nLw_NASA = rd_IP.read_NASA_data(path_NASA, bands)
    Ed_TARTU, Lsky_TARTU, Lt_TARTU, Rrs_TARTU, Rrs_std_TARTU = rd_IP.read_TARTU_data(path_TARTU, bands) # no nLw yet
    Ed_HEREON, Lsky_HEREON, Lt_HEREON, Rrs_HEREON, Rrs_std_HEREON = rd_IP.read_HEREON_data(path_HEREON, bands) # no nLw yet
    Ed_RBINS, Lsky_RBINS, Lt_RBINS, Rrs_RBINS, Rrs_std_RBINS, nLw_RBINS = rd_IP.read_RBINS_data(path_RBINS, bands)
    Ed_NOAA, Lsky_NOAA, Lt_NOAA, Rrs_NOAA, nLw_NOAA  = rd_IP.read_NOAA_data(path_NOAA, bands, Ed_PML) # PML timestamps used to reference staion no.
    Ed_CNR, Lsky_CNR, Lt_CNR, Rrs_CNR, Rrs_std_CNR, nLw_CNR = rd_IP.read_CNR_data(path_CNR, bands)
   
    nLw_TARTU, scale = rd_IP.nLw_usingPMLBRDF(nLw_PML, Rrs_PML, Rrs_TARTU,'TARTU', bands)
    nLw_HEREON, scale = rd_IP.nLw_usingPMLBRDF(nLw_PML, Rrs_PML, Rrs_HEREON,'HEREON', bands)
    

    # Read Community Procsssed data
    Ed_PML_CP, Ed_unc_PML_CP, Lsky_PML_CP , Lsky_unc_PML_CP, Lt_PML_CP, Lt_unc_PML_CP,  Rrs_PML_CP, Rrs_unc_PML_CP = rd_CP.read_CP_data('PML', Ed_PML, bands) # Ed dataframe used for timestamp
    Ed_TARTU_CP, Ed_unc_TARTU_CP, Lsky_TARTU_CP, Lsky_unc_TARTU_CP, Lt_TARTU_CP, Lt_unc_TARTU_CP,  Rrs_TARTU_CP, Rrs_unc_TARTU_CP = rd_CP.read_CP_data('TARTU', Ed_TARTU, bands)
    Ed_HEREON_CP, Ed_unc_HEREON_CP, Lsky_HEREON_CP, Lsky_unc_HEREON_CP, Lt_HEREON_CP, Lt_unc_HEREON_CP, Rrs_HEREON_CP, Rrs_unc_HEREON_CP = rd_CP.read_CP_data('HEREON', Ed_PML, bands) # PML timestamp used as HEREON did not include this (systems were close to synchronous)
    Ed_NASA_CP, Ed_unc_NASA_CP, Lsky_NASA_CP, Lsky_unc_NASA_CP, Lt_NASA_CP, Lt_unc_NASA_CP, Rrs_NASA_CP, Rrs_unc_NASA_CP = rd_CP.read_CP_data('NASA', Ed_NASA, bands) # PML timestamp used as HEREON did not include this (systems were close to synchronous)
    

    # filter by aziumuth
    Lsky_NASA_CP = QC.filter_by_azimuth_CP(Lsky_PML, Lsky_NASA, Lsky_NASA_CP) # for CP just filter NASA. IP dataframe is used for azimuth
    Lt_NASA_CP = QC.filter_by_azimuth_CP(Lt_PML, Lt_NASA, Lt_NASA_CP)
    Rrs_NASA_CP = QC.filter_by_azimuth_CP(Rrs_PML, Rrs_NASA, Rrs_NASA_CP) 
    
    Lsky_unc_NASA_CP = QC.filter_by_azimuth_CP(Lsky_PML, Lsky_NASA, Lsky_unc_NASA_CP) # for CP just filter NASA. IP dataframe is used for azimuth
    Lt_unc_NASA_CP = QC.filter_by_azimuth_CP(Lt_PML, Lt_NASA, Lt_unc_NASA_CP)
    Rrs_unc_NASA_CP = QC.filter_by_azimuth_CP(Rrs_PML, Rrs_NASA, Rrs_unc_NASA_CP) 
   
    Lsky_NASA, Lsky_RBINS, Lsky_CNR = QC.filter_by_azimuth(Lsky_PML, Lsky_NASA, Lsky_RBINS, Lsky_CNR)
    Lt_NASA, Lt_RBINS, Lt_CNR = QC.filter_by_azimuth(Lt_PML, Lt_NASA, Lt_RBINS, Lt_CNR)
    Rrs_NASA, Rrs_RBINS, Rrs_CNR = QC.filter_by_azimuth(Rrs_PML, Rrs_NASA, Rrs_RBINS, Rrs_CNR) 

    
    # Reference baselines (must be done post azimuth filtering)
    Ed_R = QC.baseline_average_V2('Ed', Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON) # baseline V2 is 4-way mean of PML, NASA, HEREON, TARTU
    Lsky_R = QC.baseline_average_V2('Lsky', Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON) # V1 - is even mixing of seabird/Trios
    Lt_R = QC.baseline_average_V2('Lt', Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON)
    Rrs_R = QC.baseline_average_V2('Rrs', Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON)
    # Rrs_std_R = QC.baseline_average_V2('Rrs', Rrs_std_PML, Rrs_std_NASA, Rrs_std_TARTU, Rrs_std_HEREON)  
    nLw_R = rd_IP.read_Aeronet_nLw(path_NLW, Ed_R, bands)    
    
    # Reference baselines -CP (must be done post azimuth filtering)
    Ed_R_CP = QC.baseline_average_V2_CP('Ed',  Ed_PML, Ed_PML_CP, Ed_NASA_CP, Ed_TARTU_CP, Ed_HEREON_CP) # baseline V2 is 4-way mean of PML, NASA, HEREON, TARTU
    Lsky_R_CP = QC.baseline_average_V2_CP('Lsky', Ed_PML, Lsky_PML_CP, Lsky_NASA_CP, Lsky_TARTU_CP, Lsky_HEREON_CP) # Ed_PML used for windspeed
    Lt_R_CP = QC.baseline_average_V2_CP('Lt', Ed_PML, Lt_PML_CP, Lt_NASA_CP, Lt_TARTU_CP, Lt_HEREON_CP)
    Rrs_R_CP = QC.baseline_average_V2_CP('Rrs',  Ed_PML, Rrs_PML_CP, Rrs_NASA_CP, Rrs_TARTU_CP, Rrs_HEREON_CP)
    
    # Quality control (dervies QC mask accomodating both SeaPrism, and Env conditions QC)
    Q_mask = QC.QC_mask(path_QC, Ed_R, Ed_PML, Lsky_PML, Rrs_PML, Rrs_std_PML, path_output) 
    # QC.plot_dataandQCmasks(Q_mask, Ed_PML, Ed_TARTU, Ed_HEREON, Ed_NASA, Ed_RBINS, Ed_CNR,Ed_NOAA, path_output)
    # QC.plot_dataandQCmasks_4teams(Q_mask, Ed_PML, Ed_TARTU, Ed_HEREON, Ed_NASA, path_output +'/IP_')
    # QC.plot_dataandQCmasks_4teams(Q_mask, Ed_PML_CP, Ed_TARTU_CP, Ed_HEREON_CP, Ed_NASA_CP, path_output +'/CP')
    # QC.azimuth_plot(Ed_PML, Ed_NASA, Ed_RBINS, Ed_CNR, path_output)
    
    ##############################################
    # Scatter & residuals plots: IP data with 7 teams
    ######################################################
    #rp.plot_scatter('Ed', Ed_R, Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_RBINS, Ed_CNR, Ed_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    #rp.plot_scatter('Lsky', Lsky_R, Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_RBINS, Lsky_CNR, Lsky_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_scatter('Lt', Lt_R, Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_RBINS, Lt_CNR, Lt_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_scatter('Rrs', Rrs_R, Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_RBINS,  Rrs_CNR, Rrs_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_scatter('nLw', nLw_R, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, nLw_RBINS, nLw_CNR, nLw_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_scatter('nLw', nLw_NOAA, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, nLw_RBINS, nLw_CNR, nLw_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')  # NOAA as baseline
    
    # rp.plot_residuals('Ed', Ed_R, Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_RBINS, Ed_CNR, Ed_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')    
    # rp.plot_residuals('Lsky', Lsky_R, Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_RBINS, Lsky_CNR, Lsky_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # rp.plot_residuals('Lt', Lt_R, Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_RBINS, Lt_CNR, Lt_NOAA, bands, path_output,  Q_mask, Qtype = 'QC_AOC_3')  
    # rp.plot_residuals('Rrs', Rrs_R, Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_RBINS, Rrs_CNR, Rrs_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_residuals('nLw', nLw_R, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, nLw_RBINS, nLw_CNR, nLw_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_residuals('nLw', nLw_NOAA, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, nLw_RBINS, nLw_CNR, nLw_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') # NOAA as baseline
    
    
    #############################################################
    # scatter & residuals plots - CP data with baseline
    ############################################################
    # rp.plot_scatter_CP('Ed', Ed_R_CP, Ed_PML_CP, Ed_NASA_CP, Ed_TARTU_CP, Ed_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    # rp.plot_scatter_CP('Lsky', Lsky_R_CP, Lsky_PML_CP, Lsky_NASA_CP, Lsky_TARTU_CP, Lsky_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # rp.plot_scatter_CP('Lt', Lt_R, Lt_PML_CP, Lt_NASA_CP, Lt_TARTU_CP, Lt_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
 #   rp.plot_scatter_CP('Rrs', Rrs_R, Rrs_PML_CP, Rrs_NASA_CP, Rrs_TARTU_CP, Rrs_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    
    # rp.plot_residuals_CP('Ed', Ed_R_CP, Ed_PML_CP, Ed_NASA_CP, Ed_TARTU_CP, Ed_HEREON_CP, Ed_unc_PML_CP, Ed_unc_NASA_CP, Ed_unc_TARTU_CP, Ed_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    # rp.plot_residuals_CP('Lsky', Lsky_R_CP, Lsky_PML_CP, Lsky_NASA_CP, Lsky_TARTU_CP, Lsky_HEREON_CP, Lsky_unc_PML_CP, Lsky_unc_NASA_CP, Lsky_unc_TARTU_CP, Lsky_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')                   
    # rp.plot_residuals_CP('Lt', Lt_R_CP, Lt_PML_CP, Lt_NASA_CP, Lt_TARTU_CP, Lt_HEREON_CP, Lt_unc_PML_CP, Lt_unc_NASA_CP, Lt_unc_TARTU_CP, Lt_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')
   # rp.plot_residuals_CP('Rrs', Rrs_R_CP, Rrs_PML_CP, Rrs_NASA_CP, Rrs_TARTU_CP, Rrs_HEREON_CP, Rrs_unc_PML_CP, Rrs_unc_NASA_CP, Rrs_unc_TARTU_CP, Rrs_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')
     
    #rp.plot_unc_CP('Ed', Ed_unc_PML_CP, Ed_unc_NASA_CP, Ed_unc_TARTU_CP, Ed_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')  
    #rp.plot_unc_CP('Lsky', Lsky_unc_PML_CP, Lsky_unc_NASA_CP, Lsky_unc_TARTU_CP, Lsky_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')  
    #rp.plot_unc_CP('Lt', Lt_unc_PML_CP, Lt_unc_NASA_CP, Lt_unc_TARTU_CP, Lt_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')  
   # rp.plot_unc_CP('Rrs', Rrs_unc_PML_CP, Rrs_unc_NASA_CP, Rrs_unc_TARTU_CP, Rrs_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')  
    
    
    ##########################################################################
    # this is `4 system versio'; of IP - Uncertainty from CP is used in plots  
    # #############################################################################
    # rp.plot_scatter_IP('Ed', Ed_R, Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    # rp.plot_scatter_IP('Lsky', Lsky_R, Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # rp.plot_scatter_IP('Lt', Lt_R, Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
   #  rp.plot_scatter_IP('Rrs', Rrs_R, Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
  #  rp.plot_scatter_IP('nLw', nLw_R, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') # SEABIRD
   # rp.plot_scatter_IP('nLw', nLw_NOAA, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') # NOAA
    
    # rp.plot_residuals_IP('Ed', Ed_R, Ed_PML_CP, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_unc_PML_CP, Ed_unc_NASA_CP, Ed_unc_TARTU_CP, Ed_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    # rp.plot_residuals_IP('Lsky', Lsky_R, Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_unc_PML_CP, Lsky_unc_NASA_CP, Lsky_unc_TARTU_CP, Lsky_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')                   
    # rp.plot_residuals_IP('Lt', Lt_R, Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_unc_PML_CP, Lt_unc_NASA_CP, Lt_unc_TARTU_CP, Lt_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')
    rp.plot_residuals_IP('Rrs', Rrs_R, Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_unc_PML_CP, Rrs_unc_NASA_CP, Rrs_unc_TARTU_CP, Rrs_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')
    rp.plot_residuals_IP_nounc('nLw', nLw_R, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')  # Seabird as baseline
    # rp.plot_residuals_IP_nounc('nLw', nLw_NOAA, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') # NOAA as baseline
    
    
    ##########################################################################
    # this is a direct comparison between CP and IP
    # #############################################################################
    # scatter & residuals plots - CP vs IP data with baseline
    # rp.plot_scatter_CP_vs_IP('Ed', Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_PML_CP, Ed_NASA_CP, Ed_TARTU_CP, Ed_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    # rp.plot_scatter_CP_vs_IP('Lsky', Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_PML_CP, Lsky_NASA_CP, Lsky_TARTU_CP, Lsky_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
   # rp.plot_scatter_CP_vs_IP('Lt', Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_PML_CP, Lt_NASA_CP, Lt_TARTU_CP, Lt_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # rp.plot_scatter_CP_vs_IP('Rrs', Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_PML_CP, Rrs_NASA_CP, Rrs_TARTU_CP, Rrs_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    
    # rp.plot_residuals_CP_vs_IP('Ed', Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_PML_CP, Ed_NASA_CP, Ed_TARTU_CP, Ed_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    # rp.plot_residuals_CP_vs_IP('Lsky', Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_PML_CP, Lsky_NASA_CP, Lsky_TARTU_CP, Lsky_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_residuals_CP_vs_IP('Lt', Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_PML_CP, Lt_NASA_CP, Lt_TARTU_CP, Lt_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_residuals_CP_vs_IP('Rrs', Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_PML_CP, Rrs_NASA_CP, Rrs_TARTU_CP, Rrs_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # sub investigation on 19/07
    
    ############################################################################
    # Summary tables - IP data
    #############################################################################
    # Ed_table = rp.tabular_summary('Ed', Ed_R, Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_RBINS, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')    
    # Lsky_table = rp.tabular_summary('Lsky', Lsky_R, Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_RBINS, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # Lt_table = rp.tabular_summary('Lt', Lt_R, Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_RBINS, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # Rrs_table = rp.tabular_summary('Rrs', Rrs_R, Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_RBINS, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # nLw = rp.tabular_summary('nLw', nLw_R, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, nLw_RBINS, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    
    
