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
import results_plots_reportversion as rp

 
if __name__ == '__main__':
    
    # path to data + output
    dir_data = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/DataSubmissions'
    path_output = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/Output_paper/rev_Feb24'

    # team submissions - 
    path_PML = dir_data + '/PML/FICE_submission_V4/FRM4SOC_2_FICE_22_AAOT_PML_HSAS_stationsummary_V4_rhoM99_nLw_localChl_noNIR_Earthsuncorrected.csv'
    path_NASA = dir_data  + '/NASA/FRM4SOC_2_FICE_22_NASA_pySAS_Sentinel3A_rev2.csv'    # NASA - no NIR/M99 - v1 in paper
    path_NASA2 = dir_data  + '/NASA/FRM4SOC_2_FICE_22_NASA_pySAS_Sentinel3A_rev3.csv'   # NASA - with NIR/Zhang17 - v2 in paper
    path_TARTU = dir_data  + '/UniTartu/averages.csv'  
    path_HEREON = dir_data + '/HEREON/FRM4SOC2_FICE22_HEREON_DATA_OLCI'
    path_HEREON_nLw = dir_data + '/HEREON/FRM4SOC2_HEREON_Normalization_revised/FRM4SOC2_FICE22_HEREON_DATA_Lwn_mean_SRF_OLCI.csv'
    path_HEREON_RAFT = dir_data + '/HEREON_RAFT/FRM4SOC2_FICE22_HEREON_RAFT_DATA_OLCI/'
    # path_NOAA = dir_data  + '/NOAA/NOAA_Hyperpro_sheet2.csv'
    path_NOAA  = dir_data + '/NOAA/NOAA_230216_sheet2.csv'
    # path_CNR = dir_data + '/CNR/FRM4SOC_2_FICE_22_AAOT_CNR_HYPSTAR_w-rel-az_ALLDATA.csv'
    path_CNR = dir_data + '/CNR/FRM4SOC_2_FICE_22_AAOT_CNR_HYPSTAR_w-rel-az_FILTERED_for_v_az.csv'
    path_RBINS = dir_data  + '/RBINS/FRM4SOC2_FICE_RBINS_2022-11-25_ViewCorrected_EarthSunCorrected.csv'
    
    # addtional data (QC + references)    
    path_QC = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/Aeronet_QC_mask/FRM4SOC-AAOT_V3_ope.txt'
    path_NLW = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/nLw_Zibordireference/'    

    # HCP output
    dir_CP_class = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/CP/FICE22-Reprocessed_12-23/' 
    dir_CP_FRM = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/CP/FICE22_FRMbranch/'
/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/CP/FICE22_v1.2.4_coserror_revised/Z17_SimSpec
    # OLCI bands
    bands = [str(400), str(412.5), str(442.5),	str(490), str(510), str(560), str(620),	str(665), str(673.75), str(681.25), str(708.75), str(753.75), str(761.25), str(764.375), str(767.5), str(778.75), str(865), str(885), str(900)]
    
    # Read data Individually Processed data. Each team has own file reader in read_IP.data.py to homogenize data
    Ed_PML, Lsky_PML, Lt_PML, Rrs_PML, Rrs_std_PML, nLw_PML = rd_IP.read_PML_data(path_PML, bands)
    Ed_NASA, Lsky_NASA, Lt_NASA, Rrs_NASA, Rrs_std_NASA, nLw_NASA = rd_IP.read_NASA_data(path_NASA, bands)
    Ed_NASA2, Lsky_NASA2, Lt_NASA2, Rrs_NASA2, Rrs_std_NASA2, nLw_NASA2 = rd_IP.read_NASA_data(path_NASA2, bands)
    Ed_TARTU, Lsky_TARTU, Lt_TARTU, Rrs_TARTU, Rrs_std_TARTU = rd_IP.read_TARTU_data(path_TARTU, bands) # no nLw yet
    Ed_HEREON, Lsky_HEREON, Lt_HEREON, Rrs_HEREON, Rrs_std_HEREON = rd_IP.read_HEREON_data(path_HEREON, bands) # no nLw yet
    Ed_RBINS, Lsky_RBINS, Lt_RBINS, Rrs_RBINS, Rrs_std_RBINS, nLw_RBINS = rd_IP.read_RBINS_data(path_RBINS, bands)
    Ed_NOAA, Lsky_NOAA, Lt_NOAA, Rrs_NOAA, nLw_NOAA  = rd_IP.read_NOAA_data_V2(path_NOAA, bands, Ed_PML) # PML timestamps used to reference staion no.
    Ed_CNR, Lsky_CNR, Lt_CNR, Rrs_CNR, Rrs_std_CNR, nLw_CNR = rd_IP.read_CNR_data(path_CNR, bands)
   
    
    nLw_TARTU, scale = rd_IP.nLw_usingPMLBRDF(nLw_PML, Rrs_PML, Rrs_TARTU,'TARTU', bands) # using PML f/q method and re-scaling
    #nLw_HEREON, scale = rd_IP.nLw_usingPMLBRDF(nLw_PML, Rrs_PML, Rrs_HEREON,'HEREON', bands) # using PML f/q method and re-scaling
    nLw_HEREON = pd.read_csv(path_HEREON_nLw,header=None,skiprows =1,usecols=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],delim_whitespace=True) # revised IP submission by HEREON
    nLw_HEREON.columns = bands
      

    # Read Community Procsssed data - convolved versions of files
    # Ed_PML_CP, Ed_unc_PML_CP, Lsky_PML_CP, Lsky_unc_PML_CP, Lt_PML_CP, Lt_unc_PML_CP,  Rrs_PML_CP, Rrs_unc_PML_CP = rd_CP.read_CP_data_convolved('PML', Ed_PML,  dir_CP_class, bands) # Ed dataframe used for timestamp
    # Ed_TARTU_CP, Ed_unc_TARTU_CP, Lsky_TARTU_CP, Lsky_unc_TARTU_CP, Lt_TARTU_CP, Lt_unc_TARTU_CP,  Rrs_TARTU_CP, Rrs_unc_TARTU_CP = rd_CP.read_CP_data('TARTU', Ed_TARTU, bands)
    class_based = False
    Ed_PML_CP_FRM, Ed_unc_PML_CP_FRM, Lsky_PML_CP_FRM , Lsky_unc_PML_CP_FRM, Lt_PML_CP_FRM, Lt_unc_PML_CP_FRM,  Rrs_PML_CP_FRM, Rrs_unc_PML_CP_FRM = rd_CP.read_CP_data('PML', Ed_PML, dir_CP_FRM, bands, class_based) # Ed dataframe used for timestamp
    Ed_TARTU_CP_FRM, Ed_unc_TARTU_CP_FRM, Lsky_TARTU_CP_FRM, Lsky_unc_TARTU_CP_FRM, Lt_TARTU_CP_FRM, Lt_unc_TARTU_CP_FRM,  Rrs_TARTU_CP_FRM, Rrs_unc_TARTU_CP_FRM = rd_CP.read_CP_data('TARTU', Ed_TARTU,  dir_CP_FRM, bands, class_based)
    Ed_HEREON_CP_FRM, Ed_unc_HEREON_CP_FRM, Lsky_HEREON_CP_FRM, Lsky_unc_HEREON_CP_FRM, Lt_HEREON_CP_FRM, Lt_unc_HEREON_CP_FRM, Rrs_HEREON_CP_FRM, Rrs_unc_HEREON_CP_FRM = rd_CP.read_CP_data('HEREON', Ed_PML, dir_CP_FRM, bands, class_based) # PML timestamp used as HEREON did not include this (systems were close to synchronous)
    Ed_NASA_CP_FRM, Ed_unc_NASA_CP_FRM, Lsky_NASA_CP_FRM, Lsky_unc_NASA_CP_FRM, Lt_NASA_CP_FRM, Lt_unc_NASA_CP_FRM, Rrs_NASA_CP_FRM, Rrs_unc_NASA_CP_FRM = rd_CP.read_CP_data('NASA', Ed_NASA,  dir_CP_FRM, bands,  class_based) # 
    
    # these are hyperspectral variants of the read CP function
    class_based = True
    Ed_PML_CP, Ed_unc_PML_CP, Lsky_PML_CP , Lsky_unc_PML_CP, Lt_PML_CP, Lt_unc_PML_CP,  Rrs_PML_CP, Rrs_unc_PML_CP = rd_CP.read_CP_data('PML', Ed_PML, dir_CP_class, bands, class_based) # Ed dataframe used for timestamp
    Ed_TARTU_CP, Ed_unc_TARTU_CP, Lsky_TARTU_CP, Lsky_unc_TARTU_CP, Lt_TARTU_CP, Lt_unc_TARTU_CP,  Rrs_TARTU_CP, Rrs_unc_TARTU_CP = rd_CP.read_CP_data('TARTU', Ed_TARTU,  dir_CP_class, bands, class_based)
    Ed_HEREON_CP, Ed_unc_HEREON_CP, Lsky_HEREON_CP, Lsky_unc_HEREON_CP, Lt_HEREON_CP, Lt_unc_HEREON_CP, Rrs_HEREON_CP, Rrs_unc_HEREON_CP = rd_CP.read_CP_data('HEREON', Ed_PML, dir_CP_class, bands, class_based) # PML timestamp used as HEREON did not include this (systems were close to synchronous)
    Ed_NASA_CP, Ed_unc_NASA_CP, Lsky_NASA_CP, Lsky_unc_NASA_CP, Lt_NASA_CP, Lt_unc_NASA_CP, Rrs_NASA_CP, Rrs_unc_NASA_CP = rd_CP.read_CP_data('NASA', Ed_NASA,  dir_CP_class, bands,  class_based) # 
    
    # filter NASA by aziumuth
    # CP class
    Lsky_NASA_CP = QC.filter_by_azimuth_CP(Lsky_PML, Lsky_NASA, Lsky_NASA_CP) # for CP just filter NASA. IP dataframe is used for azimuth
    Lt_NASA_CP = QC.filter_by_azimuth_CP(Lt_PML, Lt_NASA, Lt_NASA_CP)
    Rrs_NASA_CP = QC.filter_by_azimuth_CP(Rrs_PML, Rrs_NASA, Rrs_NASA_CP) 
    
    Lsky_unc_NASA_CP = QC.filter_by_azimuth_CP(Lsky_PML, Lsky_NASA, Lsky_unc_NASA_CP) # for CP just filter NASA. IP dataframe is used for azimuth
    Lt_unc_NASA_CP = QC.filter_by_azimuth_CP(Lt_PML, Lt_NASA, Lt_unc_NASA_CP)
    Rrs_unc_NASA_CP = QC.filter_by_azimuth_CP(Rrs_PML, Rrs_NASA, Rrs_unc_NASA_CP) 
    
    # CP FRM
    Lsky_NASA_CP_FRM = QC.filter_by_azimuth_CP(Lsky_PML, Lsky_NASA, Lsky_NASA_CP_FRM) # for CP just filter NASA. IP dataframe is used for azimuth
    Lt_NASA_CP_FRM = QC.filter_by_azimuth_CP(Lt_PML, Lt_NASA, Lt_NASA_CP_FRM)
    Rrs_NASA_CP_FRM = QC.filter_by_azimuth_CP(Rrs_PML, Rrs_NASA, Rrs_NASA_CP_FRM) 
    
    Lsky_unc_NASA_CP_FRM = QC.filter_by_azimuth_CP(Lsky_PML, Lsky_NASA, Lsky_unc_NASA_CP_FRM) # for CP just filter NASA. IP dataframe is used for azimuth
    Lt_unc_NASA_CP_FRM = QC.filter_by_azimuth_CP(Lt_PML, Lt_NASA, Lt_unc_NASA_CP_FRM)
    Rrs_unc_NASA_CP_FRM = QC.filter_by_azimuth_CP(Rrs_PML, Rrs_NASA, Rrs_unc_NASA_CP_FRM) 
   
 
    # IP
    Lsky_NASA, Lsky_RBINS, Lsky_CNR = QC.filter_by_azimuth(Lsky_PML, Lsky_NASA, Lsky_RBINS, Lsky_CNR)
    Lt_NASA, Lt_RBINS, Lt_CNR = QC.filter_by_azimuth(Lt_PML, Lt_NASA, Lt_RBINS, Lt_CNR)
    Rrs_NASA, Rrs_RBINS, Rrs_CNR = QC.filter_by_azimuth(Rrs_PML, Rrs_NASA, Rrs_RBINS, Rrs_CNR) 

    Lsky_NASA, Lsky_RBINS = QC.filter_by_azimuth_noCNR(Lsky_PML, Lsky_NASA, Lsky_RBINS)
    Lt_NASA, Lt_RBINS = QC.filter_by_azimuth_noCNR(Lt_PML, Lt_NASA, Lt_RBINS)
    Rrs_NASA, Rrs_RBINS = QC.filter_by_azimuth_noCNR(Rrs_PML, Rrs_NASA, Rrs_RBINS) 
    
    # manual removal of final cast of day on 15, 19, 20th for NOAA HP
    Ed_NOAA.iloc[28,:] = np.nan
    Ed_NOAA.iloc[53,:] = np.nan
    Ed_NOAA.iloc[69,:] = np.nan

    #
    nLw_NOAA.iloc[28,:] = np.nan
    nLw_NOAA.iloc[53,:] = np.nan
    nLw_NOAA.iloc[69,:] = np.nan

    
    # Reference baselines (must be done post azimuth filtering)
    Ed_R = QC.baseline_average_V2('Ed', Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON) # baseline V2 is 4-way mean of PML, NASA, HEREON, TARTU
    Lsky_R = QC.baseline_average_V2('Lsky', Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON) # V1 - is even mixing of seabird/Trios
    Lt_R = QC.baseline_average_V2('Lt', Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON)
    Rrs_R = QC.baseline_average_V2('Rrs', Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON)
    
    
    # Rrs_std_R = QC.baseline_average_V2('Rrs', Rrs_std_PML, Rrs_std_NASA, Rrs_std_TARTU, Rrs_std_HEREON)  
    nLw_SEAPRISM = rd_IP.read_Aeronet_nLw(path_NLW, Ed_R, bands)    
    nLw_RAFT =  rd_IP.read_HEREON_RAFT_data(path_HEREON_RAFT, bands)    
    
    # Reference baselines P (must be done post azimuth filtering)
    Ed_R_CP = QC.baseline_average_V2_CP('Ed',  Ed_PML,  Ed_PML_CP, Ed_NASA_CP, Ed_TARTU_CP, Ed_HEREON_CP) # baseline V2 is 4-way mean of PML, NASA, HEREON, TARTU
    Lsky_R_CP = QC.baseline_average_V2_CP('Lsky', Ed_PML, Lsky_PML_CP, Lsky_NASA_CP, Lsky_TARTU_CP, Lsky_HEREON_CP) # Ed_PML used for windspeed
    Lt_R_CP = QC.baseline_average_V2_CP('Lt', Ed_PML, Lt_PML_CP, Lt_NASA_CP, Lt_TARTU_CP, Lt_HEREON_CP)
    Rrs_R_CP = QC.baseline_average_V2_CP('Rrs', Ed_PML, Rrs_PML_CP, Rrs_NASA_CP, Rrs_TARTU_CP, Rrs_HEREON_CP)
    
    Ed_R_CP_FRM = QC.baseline_average_V2_CP('Ed',  Ed_PML,  Ed_PML_CP_FRM, Ed_NASA_CP_FRM, Ed_TARTU_CP_FRM, Ed_HEREON_CP_FRM)
    Lsky_R_CP_FRM = QC.baseline_average_V2_CP('Lsky', Lsky_PML, Lsky_PML_CP_FRM, Lsky_NASA_CP_FRM, Lsky_TARTU_CP_FRM, Lsky_HEREON_CP_FRM)
    Lt_R_CP_FRM = QC.baseline_average_V2_CP('Lt',  Lt_PML,  Lt_PML_CP_FRM, Lt_NASA_CP_FRM, Lt_TARTU_CP_FRM, Lt_HEREON_CP_FRM)
    Rrs_R_CP_FRM = QC.baseline_average_V2_CP('Rrs', Rrs_PML,  Rrs_PML_CP_FRM, Rrs_NASA_CP_FRM, Rrs_TARTU_CP_FRM, Rrs_HEREON_CP_FRM)
    
    
    # Quality control (dervies QC mask accomodating both SeaPrism, and Env conditions QC)
    Q_mask = QC.QC_mask(path_QC, Ed_R, Ed_PML, Lsky_PML, Rrs_PML, Rrs_std_PML, path_output) 
    
    # QC.plot_dataandQCmasks(Q_mask, Ed_PML, Ed_TARTU, Ed_HEREON, Ed_NASA, Ed_RBINS, Ed_CNR,Ed_NOAA, path_output)
    # QC.plot_dataandQCmasks_4teams(Q_mask, Ed_PML, Ed_TARTU, Ed_HEREON, Ed_NASA, path_output +'/IP_')
    # QC.plot_dataandQCmasks_4teams(Q_mask, Ed_PML_CP, Ed_TARTU_CP, Ed_HEREON_CP, Ed_NASA_CP, path_output +'/CP')
    # QC.azimuth_plot(Ed_PML, Ed_NASA, Ed_RBINS, Ed_CNR, path_output)
    
    ##############################################
    # Scatter & residuals plots: IP data with 7 teams
    ######################################################
    # rp.plot_scatter('Ed', Ed_R, Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_RBINS, Ed_CNR, Ed_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    #rp.plot_scatter('Lsky', Lsky_R, Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_RBINS, Lsky_CNR, Lsky_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_scatter('Lt', Lt_R, Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_RBINS, Lt_CNR, Lt_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_scatter('Rrs', Rrs_R, Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_RBINS,  Rrs_CNR, Rrs_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_scatter('nLw', nLw_R, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, nLw_RBINS, nLw_CNR, nLw_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_scatter('nLw', nLw_NOAA, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, nLw_RBINS, nLw_CNR, nLw_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')  # NOAA as baseline
    
    #rp.plot_residuals('Ed', Ed_R, Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON,Ed_RBINS, Ed_CNR, Ed_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')
    #rp.plot_residuals('Lsky', Lsky_R, Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_RBINS, Lsky_CNR, Lsky_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    #rp.plot_residuals('Lt', Lt_R, Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_RBINS, Lt_CNR, Lt_NOAA, bands, path_output,  Q_mask, Qtype = 'QC_AOC_3')  
    #rp.plot_residuals('Rrs', Rrs_R, Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_RBINS, Rrs_CNR, Rrs_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # rp.plot_residuals('nLw', nLw_R, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, nLw_RBINS, nLw_CNR, nLw_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
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
    # rp.plot_residuals_IP('Rrs', Rrs_R, Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_unc_PML_CP, Rrs_unc_NASA_CP, Rrs_unc_TARTU_CP, Rrs_unc_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')
    #rp.plot_residuals_IP_nounc('nLw', nLw_R, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')  # Seabird as baseline
    # rp.plot_residuals_IP_nounc('nLw', nLw_NOAA, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') # NOAA as baseline
    
    
    ######################################################### 9, r'A.', fontsize=18)
  #################
    # this is a direct comparison between CP and IP
    # #############################################################################
    # scatter & residuals plots - CP vs IP data with baseline
    # rp.plot_scatter_CP_vs_IP('Ed', Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_PML_CP, Ed_NASA_CP, Ed_TARTU_CP, Ed_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    # rp.plot_scatter_CP_vs_IP('Lsky', Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_PML_CP, Lsky_NASA_CP, Lsky_TARTU_CP, Lsky_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # rp.plot_scatter_CP_vs_IP('Lt', Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_PML_CP, Lt_NASA_CP, Lt_TARTU_CP, Lt_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # rp.plot_scatter_CP_vs_IP('Rrs', Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_PML_CP, Rrs_NASA_CP, Rrs_TARTU_CP, Rrs_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    
    # rp.plot_residuals_CP_vs_IP('Ed', Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_PML_CP, Ed_NASA_CP, Ed_TARTU_CP, Ed_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    # rp.plot_residuals_CP_vs_IP('Lsky', Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_PML_CP, Lsky_NASA_CP, Lsky_TARTU_CP, Lsky_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # rp.plot_residuals_CP_vs_IP('Lt', Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_PML_CP, Lt_NASA_CP, Lt_TARTU_CP, Lt_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # rp.plot_residuals_CP_vs_IP('Rrs', Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_PML_CP, Rrs_NASA_CP, Rrs_TARTU_CP, Rrs_HEREON_CP, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # sub investigation on 19/07
    

    ##########################################################################
    # these are summmary ` combined' plots for the OE paper draft: Row 1: HCP, Row 2: Uncertanties, Row 3: HCP versus IP
    # #############################################################################
 
    #rp.residuals_combined('Ed', Ed_R_CP, Ed_PML_CP, Ed_NASA_CP, Ed_TARTU_CP, Ed_HEREON_CP, 
    #                     Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON,
     #                     Ed_unc_PML_CP, Ed_unc_NASA_CP, Ed_unc_TARTU_CP, Ed_unc_HEREON_CP, 
     #                     bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    
   # rp.residuals_combined('Lsky', Lsky_R_CP, Lsky_PML_CP, Lsky_NASA_CP, Lsky_TARTU_CP, Lsky_HEREON_CP, 
    #                         Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, 
     #                        Lsky_unc_PML_CP, Lsky_unc_NASA_CP, Lsky_unc_TARTU_CP, Lsky_unc_HEREON_CP, 
      #                       bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
   
   # rp.residuals_combined('Lt', Lt_R_CP, Lt_PML_CP, Lt_NASA_CP, Lt_TARTU_CP, Lt_HEREON_CP, 
     #                        Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, 
       #                      Lt_unc_PML_CP, Lt_unc_NASA_CP, Lt_unc_TARTU_CP, Lt_unc_HEREON_CP, 
        #                     bands, path_output, Q_mask, Qtype = 'QC_AOC_3')
    
    
   # rp.residuals_combined_('Rrs', Rrs_R_CP, Rrs_PML_CP, Rrs_NASA_CP, Rrs_TARTU_CP, Rrs_HEREON_CP, 
        #                     Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON,
           #                  Rrs_unc_PML_CP, Rrs_unc_NASA_CP, Rrs_unc_TARTU_CP, Rrs_unc_HEREON_CP, 
              #               bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    
   # rp.residuals_combined_nlw('nLw', nLw_SEAPRISM, nLw_NOAA, nLw_PML, nLw_NASA, nLw_NASA2, nLw_TARTU, nLw_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
  #  rp.residuals_combined_nlw_v2('nLw', nLw_SEAPRISM, nLw_NOAA, nLw_RAFT, nLw_PML, nLw_NASA, nLw_NASA2, nLw_TARTU, nLw_HEREON, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')  
    # rp.SB_VS_HP('nLw', nLw_SEAPRISM, nLw_NOAA, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')
    # rp.SB_VS_RAFT('nLw', nLw_SEAPRISM, nLw_RAFT, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')

#    rp.FRM_unc_plot(Ed_unc_PML_CP_FRM, Ed_unc_NASA_CP_FRM, Ed_unc_TARTU_CP_FRM, Ed_unc_HEREON_CP_FRM,
 #                      Lsky_unc_PML_CP_FRM, Lsky_unc_NASA_CP_FRM, Lsky_unc_TARTU_CP_FRM, Lsky_unc_HEREON_CP_FRM,
  #                     Lt_unc_PML_CP_FRM, Lt_unc_NASA_CP_FRM, Lt_unc_TARTU_CP_FRM, Lt_unc_HEREON_CP_FRM,
   #                    Rrs_unc_PML_CP_FRM, Rrs_unc_NASA_CP_FRM, Rrs_unc_TARTU_CP_FRM, Rrs_unc_HEREON_CP_FRM,
    #                   bands, path_output, Q_mask, Qtype = 'QC_AOC_3')


    # uncertainty differences - un-normalized
    # diff_Ed_unc_PML = Ed_unc_PML_CP.iloc[:,1:] - Ed_unc_PML_CP_FRM.iloc[:,1:]
    # diff_Ed_unc_NASA = Ed_unc_NASA_CP.iloc[:,1:] - Ed_unc_NASA_CP_FRM.iloc[:,1:]
    # diff_Ed_unc_TARTU = Ed_unc_TARTU_CP.iloc[:,1:] - Ed_unc_TARTU_CP_FRM.iloc[:,1:]
    # diff_Ed_unc_HEREON = Ed_unc_HEREON_CP.iloc[:,1:] - Ed_unc_HEREON_CP_FRM.iloc[:,1:]
    
    #  diff_Lsky_unc_PML = Lsky_unc_PML_CP.iloc[:,1:] - Lsky_unc_PML_CP_FRM.iloc[:,1:]
    # diff_Lsky_unc_NASA = Lsky_unc_NASA_CP.iloc[:,1:] - Lsky_unc_NASA_CP_FRM.iloc[:,1:]
    # diff_Lsky_unc_TARTU = Lsky_unc_TARTU_CP.iloc[:,1:] - Lsky_unc_TARTU_CP_FRM.iloc[:,1:]
    # diff_Lsky_unc_HEREON = Lsky_unc_HEREON_CP.iloc[:,1:] - Lsky_unc_HEREON_CP_FRM.iloc[:,1:]
     
    # diff_Lt_unc_PML = Lt_unc_PML_CP.iloc[:,1:] - Lt_unc_PML_CP_FRM.iloc[:,1:]
    # diff_Lt_unc_NASA = Lt_unc_NASA_CP.iloc[:,1:] - Lt_unc_NASA_CP_FRM.iloc[:,1:]
    # diff_Lt_unc_TARTU = Lt_unc_TARTU_CP.iloc[:,1:] - Lt_unc_TARTU_CP_FRM.iloc[:,1:]
    # diff_Lt_unc_HEREON = Lt_unc_HEREON_CP.iloc[:,1:] - Lt_unc_HEREON_CP_FRM.iloc[:,1:]
    
   # diff_Rrs_unc_PML = Rrs_unc_PML_CP.iloc[:,1:] - Rrs_unc_PML_CP_FRM.iloc[:,1:]
   # diff_Rrs_unc_NASA = Rrs_unc_NASA_CP.iloc[:,1:] - Rrs_unc_NASA_CP_FRM.iloc[:,1:]
    #diff_Rrs_unc_TARTU = Rrs_unc_TARTU_CP.iloc[:,1:] - Rrs_unc_TARTU_CP_FRM.iloc[:,1:]
    #diff_Rrs_unc_HEREON = Rrs_unc_HEREON_CP.iloc[:,1:] - Rrs_unc_HEREON_CP_FRM.iloc[:,1:]



    #  rp.FRM_uncdiff_plot(diff_Ed_unc_PML, diff_Ed_unc_NASA, diff_Ed_unc_TARTU, diff_Ed_unc_HEREON,  
    #                  diff_Lsky_unc_PML, diff_Lsky_unc_NASA, diff_Lsky_unc_TARTU, diff_Lsky_unc_HEREON, 
    #                 diff_Lt_unc_PML, diff_Lt_unc_NASA, diff_Lt_unc_TARTU, diff_Lt_unc_HEREON, 
    #                diff_Rrs_unc_PML, diff_Rrs_unc_NASA, diff_Rrs_unc_TARTU, diff_Rrs_unc_HEREON, 
           
                    
    #     bands, path_output, Q_mask, Qtype = 'QC_AOC_3')

   # diff_Ed_unc_PML = 100*(Ed_unc_PML_CP.iloc[:,1:] - Ed_unc_PML_CP_FRM.iloc[:,1:])/(Ed_unc_PML_CP_FRM.iloc[:,1:])
   # diff_Ed_unc_NASA = 100*(Ed_unc_NASA_CP.iloc[:,1:] - Ed_unc_NASA_CP_FRM.iloc[:,1:])/(Ed_unc_NASA_CP_FRM.iloc[:,1:])
    #diff_Ed_unc_TARTU = 100*(Ed_unc_TARTU_CP.iloc[:,1:] - Ed_unc_TARTU_CP_FRM.iloc[:,1:])/(Ed_unc_TARTU_CP_FRM.iloc[:,1:])
   # diff_Ed_unc_HEREON = 100*(Ed_unc_HEREON_CP.iloc[:,1:] - Ed_unc_HEREON_CP_FRM.iloc[:,1:])/(Ed_unc_HEREON_CP_FRM.iloc[:,1:])
    
   # diff_Lsky_unc_PML = 100*(Lsky_unc_PML_CP.iloc[:,1:] - Lsky_unc_PML_CP_FRM.iloc[:,1:])/(Lsky_unc_PML_CP_FRM.iloc[:,1:])
   # diff_Lsky_unc_NASA = 100*(Lsky_unc_NASA_CP.iloc[:,1:] - Lsky_unc_NASA_CP_FRM.iloc[:,1:])/(Lsky_unc_NASA_CP_FRM.iloc[:,1:])
    #diff_Lsky_unc_TARTU = 100*(Lsky_unc_TARTU_CP.iloc[:,1:] - Lsky_unc_TARTU_CP_FRM.iloc[:,1:])/(Lsky_unc_TARTU_CP_FRM.iloc[:,1:])
    #diff_Lsky_unc_HEREON = 100*(Lsky_unc_HEREON_CP.iloc[:,1:] - Lsky_unc_HEREON_CP_FRM.iloc[:,1:])/(Lsky_unc_HEREON_CP_FRM.iloc[:,1:])
     
    #diff_Lt_unc_PML = 100*(Lt_unc_PML_CP.iloc[:,1:] - Lt_unc_PML_CP_FRM.iloc[:,1:])/(Lt_unc_PML_CP_FRM.iloc[:,1:])
    #diff_Lt_unc_NASA = 100*(Lt_unc_NASA_CP.iloc[:,1:] - Lt_unc_NASA_CP_FRM.iloc[:,1:])/(Lt_unc_NASA_CP_FRM.iloc[:,1:])
    #diff_Lt_unc_TARTU = 100*(Lt_unc_TARTU_CP.iloc[:,1:] - Lt_unc_TARTU_CP_FRM.iloc[:,1:])/(Lt_unc_TARTU_CP_FRM.iloc[:,1:])
    #diff_Lt_unc_HEREON = 100*(Lt_unc_HEREON_CP.iloc[:,1:] - Lt_unc_HEREON_CP_FRM.iloc[:,1:])/(Lt_unc_HEREON_CP_FRM.iloc[:,1:])
    
   # diff_Rrs_unc_PML = 100*(Rrs_unc_PML_CP.iloc[:,1:] - Rrs_unc_PML_CP_FRM.iloc[:,1:])/(Rrs_unc_PML_CP_FRM.iloc[:,1:])
   # diff_Rrs_unc_NASA = 100*(Rrs_unc_NASA_CP.iloc[:,1:] - Rrs_unc_NASA_CP_FRM.iloc[:,1:])/(Rrs_unc_NASA_CP_FRM.iloc[:,1:])
  #  diff_Rrs_unc_TARTU = 100*(Rrs_unc_TARTU_CP.iloc[:,1:] - Rrs_unc_TARTU_CP_FRM.iloc[:,1:])/(Rrs_unc_TARTU_CP_FRM.iloc[:,1:])
  #  diff_Rrs_unc_HEREON = 100*(Rrs_unc_HEREON_CP.iloc[:,1:] - Rrs_unc_HEREON_CP_FRM.iloc[:,1:])/(Rrs_unc_HEREON_CP_FRM.iloc[:,1:])


    diff_Ed_unc_PML = 200*(Ed_unc_PML_CP.iloc[:,1:] - Ed_unc_PML_CP_FRM.iloc[:,1:])/(Ed_unc_PML_CP.iloc[:,1:] + Ed_unc_PML_CP_FRM.iloc[:,1:])
    diff_Ed_unc_NASA = 200*(Ed_unc_NASA_CP.iloc[:,1:] - Ed_unc_NASA_CP_FRM.iloc[:,1:])/(Ed_unc_NASA_CP.iloc[:,1:] + Ed_unc_NASA_CP_FRM.iloc[:,1:])
    diff_Ed_unc_TARTU = 200*(Ed_unc_TARTU_CP.iloc[:,1:] - Ed_unc_TARTU_CP_FRM.iloc[:,1:])/(Ed_unc_TARTU_CP.iloc[:,1:] + Ed_unc_TARTU_CP_FRM.iloc[:,1:])
    diff_Ed_unc_HEREON = 200*(Ed_unc_HEREON_CP.iloc[:,1:] - Ed_unc_HEREON_CP_FRM.iloc[:,1:])/(Ed_unc_HEREON_CP.iloc[:,1:] + Ed_unc_HEREON_CP_FRM.iloc[:,1:])
    
   
    diff_Lsky_unc_PML = 200*(Lsky_unc_PML_CP.iloc[:,1:] - Lsky_unc_PML_CP_FRM.iloc[:,1:])/(Lsky_unc_PML_CP.iloc[:,1:] + Lsky_unc_PML_CP_FRM.iloc[:,1:])
    diff_Lsky_unc_NASA = 200*(Lsky_unc_NASA_CP.iloc[:,1:] - Lsky_unc_NASA_CP_FRM.iloc[:,1:])/(Lsky_unc_NASA_CP.iloc[:,1:] + Lsky_unc_NASA_CP_FRM.iloc[:,1:])
    diff_Lsky_unc_TARTU = 200*(Lsky_unc_TARTU_CP.iloc[:,1:] - Lsky_unc_TARTU_CP_FRM.iloc[:,1:])/(Lsky_unc_TARTU_CP.iloc[:,1:] + Lsky_unc_TARTU_CP_FRM.iloc[:,1:])
    diff_Lsky_unc_HEREON = 200*(Lsky_unc_HEREON_CP.iloc[:,1:] - Lsky_unc_HEREON_CP_FRM.iloc[:,1:])/(Lsky_unc_HEREON_CP.iloc[:,1:] + Lsky_unc_HEREON_CP_FRM.iloc[:,1:])
     
    diff_Lt_unc_PML = 200*(Lt_unc_PML_CP.iloc[:,1:] - Lt_unc_PML_CP_FRM.iloc[:,1:])/(Lt_unc_PML_CP.iloc[:,1:] + Lt_unc_PML_CP_FRM.iloc[:,1:])
    diff_Lt_unc_NASA = 200*(Lt_unc_NASA_CP.iloc[:,1:] - Lt_unc_NASA_CP_FRM.iloc[:,1:])/(Lt_unc_NASA_CP.iloc[:,1:] + Lt_unc_NASA_CP_FRM.iloc[:,1:])
    diff_Lt_unc_TARTU = 200*(Lt_unc_TARTU_CP.iloc[:,1:] - Lt_unc_TARTU_CP_FRM.iloc[:,1:])/(Lt_unc_TARTU_CP.iloc[:,1:] + Lt_unc_TARTU_CP_FRM.iloc[:,1:])
    diff_Lt_unc_HEREON = 200*(Lt_unc_HEREON_CP.iloc[:,1:] - Lt_unc_HEREON_CP_FRM.iloc[:,1:])/(Lt_unc_HEREON_CP.iloc[:,1:] + Lt_unc_HEREON_CP_FRM.iloc[:,1:])
    
    diff_Rrs_unc_PML = 200*(Rrs_unc_PML_CP.iloc[:,1:] - Rrs_unc_PML_CP_FRM.iloc[:,1:])/(Rrs_unc_PML_CP.iloc[:,1:] + Rrs_unc_PML_CP_FRM.iloc[:,1:])
    diff_Rrs_unc_NASA = 200*(Rrs_unc_NASA_CP.iloc[:,1:] - Rrs_unc_NASA_CP_FRM.iloc[:,1:])/(Rrs_unc_NASA_CP.iloc[:,1:] + Rrs_unc_NASA_CP_FRM.iloc[:,1:])
    diff_Rrs_unc_TARTU = 200*(Rrs_unc_TARTU_CP.iloc[:,1:] - Rrs_unc_TARTU_CP_FRM.iloc[:,1:])/(Rrs_unc_TARTU_CP.iloc[:,1:] + Rrs_unc_TARTU_CP_FRM.iloc[:,1:])
    diff_Rrs_unc_HEREON = 200*(Rrs_unc_HEREON_CP.iloc[:,1:] - Rrs_unc_HEREON_CP_FRM.iloc[:,1:])/(Rrs_unc_HEREON_CP.iloc[:,1:] + Rrs_unc_HEREON_CP_FRM.iloc[:,1:])



    # rp.FRM_uncdiffV2_plot(diff_Ed_unc_PML, diff_Ed_unc_NASA, diff_Ed_unc_TARTU, diff_Ed_unc_HEREON,  
    #                 diff_Lsky_unc_PML, diff_Lsky_unc_NASA, diff_Lsky_unc_TARTU, diff_Lsky_unc_HEREON, 
    #                diff_Lt_unc_PML, diff_Lt_unc_NASA, diff_Lt_unc_TARTU, diff_Lt_unc_HEREON, 
    #                diff_Rrs_unc_PML, diff_Rrs_unc_NASA, diff_Rrs_unc_TARTU, diff_Rrs_unc_HEREON, 
    #                  bands, path_output, Q_mask, Qtype = 'QC_AOC_3')

  #  diff_Ed_unc_PML = 200*(Ed_unc_PML_CP.iloc[:,1:] - Ed_unc_PML_CP_FRM.iloc[:,1:])/(Ed_unc_PML_CP.iloc[:,1:] + Ed_unc_PML_CP_FRM.iloc[:,1:])
   # diff_Ed_unc_NASA = 200*(Ed_unc_NASA_CP.iloc[:,1:] - Ed_unc_NASA_CP_FRM.iloc[:,1:])/(Ed_unc_NASA_CP.iloc[:,1:] + Ed_unc_NASA_CP_FRM.iloc[:,1:])
    #diff_Ed_unc_TARTU = 200*(Ed_unc_TARTU_CP.iloc[:,1:] - Ed_unc_TARTU_CP_FRM.iloc[:,1:])/(Ed_unc_TARTU_CP.iloc[:,1:] + Ed_unc_TARTU_CP_FRM.iloc[:,1:])
    #diff_Ed_unc_HEREON = 200*(Ed_unc_HEREON_CP.iloc[:,1:] - Ed_unc_HEREON_CP_FRM.iloc[:,1:])/(Ed_unc_HEREON_CP.iloc[:,1:] + Ed_unc_HEREON_CP_FRM.iloc[:,1:])
    
   
    # diff_Lsky_unc_PML = 200*(Lsky_unc_PML_CP.iloc[:,1:] - Lsky_unc_PML_CP_FRM.iloc[:,1:])/(Lsky_unc_PML_CP.iloc[:,1:] + Lsky_unc_PML_CP_FRM.iloc[:,1:])
    # diff_Lsky_unc_NASA = 200*(Lsky_unc_NASA_CP.iloc[:,1:] - Lsky_unc_NASA_CP_FRM.iloc[:,1:])/(Lsky_unc_NASA_CP.iloc[:,1:] + Lsky_unc_NASA_CP_FRM.iloc[:,1:])
    # diff_Lsky_unc_TARTU = 200*(Lsky_unc_TARTU_CP.iloc[:,1:] - Lsky_unc_TARTU_CP_FRM.iloc[:,1:])/(Lsky_unc_TARTU_CP.iloc[:,1:] + Lsky_unc_TARTU_CP_FRM.iloc[:,1:])
    # diff_Lsky_unc_HEREON = 200*(Lsky_unc_HEREON_CP.iloc[:,1:] - Lsky_unc_HEREON_CP_FRM.iloc[:,1:])/(Lsky_unc_HEREON_CP.iloc[:,1:] + Lsky_unc_HEREON_CP_FRM.iloc[:,1:])
     
    # diff_Lt_unc_PML = 200*(Lt_unc_PML_CP.iloc[:,1:] - Lt_unc_PML_CP_FRM.iloc[:,1:])/(Lt_unc_PML_CP.iloc[:,1:] + Lt_unc_PML_CP_FRM.iloc[:,1:])
    # diff_Lt_unc_NASA = 200*(Lt_unc_NASA_CP.iloc[:,1:] - Lt_unc_NASA_CP_FRM.iloc[:,1:])/(Lt_unc_NASA_CP.iloc[:,1:] + Lt_unc_NASA_CP_FRM.iloc[:,1:])
    # diff_Lt_unc_TARTU = 200*(Lt_unc_TARTU_CP.iloc[:,1:] - Lt_unc_TARTU_CP_FRM.iloc[:,1:])/(Lt_unc_TARTU_CP.iloc[:,1:] + Lt_unc_TARTU_CP_FRM.iloc[:,1:])
    # diff_Lt_unc_HEREON = 200*(Lt_unc_HEREON_CP.iloc[:,1:] - Lt_unc_HEREON_CP_FRM.iloc[:,1:])/(Lt_unc_HEREON_CP.iloc[:,1:] + Lt_unc_HEREON_CP_FRM.iloc[:,1:])
    
    # diff_Rrs_unc_PML = 200*(Rrs_unc_PML_CP.iloc[:,1:] - Rrs_unc_PML_CP_FRM.iloc[:,1:])/(Rrs_unc_PML_CP.iloc[:,1:] + Rrs_unc_PML_CP_FRM.iloc[:,1:])
    # diff_Rrs_unc_NASA = 200*(Rrs_unc_NASA_CP.iloc[:,1:] - Rrs_unc_NASA_CP_FRM.iloc[:,1:])/(Rrs_unc_NASA_CP.iloc[:,1:] + Rrs_unc_NASA_CP_FRM.iloc[:,1:])
    # diff_Rrs_unc_TARTU = 200*(Rrs_unc_TARTU_CP.iloc[:,1:] - Rrs_unc_TARTU_CP_FRM.iloc[:,1:])/(Rrs_unc_TARTU_CP.iloc[:,1:] + Rrs_unc_TARTU_CP_FRM.iloc[:,1:])
    # diff_Rrs_unc_HEREON = 200*(Rrs_unc_HEREON_CP.iloc[:,1:] - Rrs_unc_HEREON_CP_FRM.iloc[:,1:])/(Rrs_unc_HEREON_CP.iloc[:,1:] + Rrs_unc_HEREON_CP_FRM.iloc[:,1:])


    #
    diff_Ed_PML = 100*(Ed_PML_CP.iloc[:,1:] - Ed_PML_CP_FRM.iloc[:,1:])/(Ed_PML_CP_FRM.iloc[:,1:])
    diff_Ed_NASA = 100*(Ed_NASA_CP.iloc[:,1:] - Ed_NASA_CP_FRM.iloc[:,1:])/(Ed_NASA_CP_FRM.iloc[:,1:])
    diff_Ed_TARTU = 100*(Ed_TARTU_CP.iloc[:,1:] - Ed_TARTU_CP_FRM.iloc[:,1:])/(Ed_TARTU_CP_FRM.iloc[:,1:])
    diff_Ed_HEREON = 100*(Ed_HEREON_CP.iloc[:,1:] - Ed_HEREON_CP_FRM.iloc[:,1:])/(Ed_HEREON_CP_FRM.iloc[:,1:])
   
    diff_Lsky_PML = 100*(Lsky_PML_CP.iloc[:,1:] - Lsky_PML_CP_FRM.iloc[:,1:])/(Lsky_PML_CP_FRM.iloc[:,1:])
    diff_Lsky_NASA = 100*(Lsky_NASA_CP.iloc[:,1:] - Lsky_NASA_CP_FRM.iloc[:,1:])/(Lsky_NASA_CP_FRM.iloc[:,1:])
    diff_Lsky_TARTU = 100*(Lsky_TARTU_CP.iloc[:,1:] - Lsky_TARTU_CP_FRM.iloc[:,1:])/(Lsky_TARTU_CP_FRM.iloc[:,1:])
    diff_Lsky_HEREON = 100*(Lsky_HEREON_CP.iloc[:,1:] - Lsky_HEREON_CP_FRM.iloc[:,1:])/(Lsky_HEREON_CP_FRM.iloc[:,1:])
    
    diff_Lt_PML = 100*(Lt_PML_CP.iloc[:,1:] - Lt_PML_CP_FRM.iloc[:,1:])/(Lt_PML_CP_FRM.iloc[:,1:])
    diff_Lt_NASA = 100*(Lt_NASA_CP.iloc[:,1:] - Lt_NASA_CP_FRM.iloc[:,1:])/(Lt_NASA_CP_FRM.iloc[:,1:])
    diff_Lt_TARTU = 100*(Lt_TARTU_CP.iloc[:,1:] - Lt_TARTU_CP_FRM.iloc[:,1:])/(Lt_TARTU_CP_FRM.iloc[:,1:])
    diff_Lt_HEREON = 100*(Lt_HEREON_CP.iloc[:,1:] - Lt_HEREON_CP_FRM.iloc[:,1:])/(Lt_HEREON_CP_FRM.iloc[:,1:])
   
    diff_Rrs_PML = 100*(Rrs_PML_CP.iloc[:,1:] - Rrs_PML_CP_FRM.iloc[:,1:])/(Rrs_PML_CP_FRM.iloc[:,1:])
    diff_Rrs_NASA = 100*(Rrs_NASA_CP.iloc[:,1:] - Rrs_NASA_CP_FRM.iloc[:,1:])/(Rrs_NASA_CP_FRM.iloc[:,1:])
    diff_Rrs_TARTU = 100*(Rrs_TARTU_CP.iloc[:,1:] - Rrs_TARTU_CP_FRM.iloc[:,1:])/(Rrs_TARTU_CP_FRM.iloc[:,1:])
    diff_Rrs_HEREON = 100*(Rrs_HEREON_CP.iloc[:,1:] - Rrs_HEREON_CP_FRM.iloc[:,1:])/(Rrs_HEREON_CP_FRM.iloc[:,1:])

    #
    rp.FRM_rad_diffplot(diff_Ed_PML, diff_Ed_NASA, diff_Ed_TARTU, diff_Ed_HEREON,  
                 diff_Lsky_PML, diff_Lsky_NASA, diff_Lsky_TARTU, diff_Lsky_HEREON, 
                 diff_Lt_PML, diff_Lt_NASA, diff_Lt_TARTU, diff_Lt_HEREON, 
                 diff_Rrs_PML, diff_Rrs_NASA, diff_Rrs_TARTU, diff_Rrs_HEREON, 
                 bands, path_output, Q_mask, Qtype = 'QC_AOC_3')


    ##########################################################################
    # these are summmary ` combined' plots for the OE paper draft: Row 1: HCP, Row 2:  Class Uncertanties, Row 3: HCP versus IP
    # #############################################################################
 
    
    # currently set to plot sensor specifc
    rp.residuals_combined_V3('Ed',  Ed_R_CP,  Ed_PML_CP,  Ed_NASA_CP,  Ed_TARTU_CP,  Ed_HEREON_CP,
                          Ed_R_CP_FRM, Ed_PML_CP_FRM, Ed_NASA_CP_FRM, Ed_TARTU_CP_FRM, Ed_HEREON_CP_FRM, 
                          Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON,
                          diff_Ed_PML,  diff_Ed_NASA, diff_Ed_TARTU, diff_Ed_HEREON,
                          Ed_unc_PML_CP, Ed_unc_NASA_CP, Ed_unc_TARTU_CP, Ed_unc_HEREON_CP, 
                          Ed_unc_PML_CP_FRM, Ed_unc_NASA_CP_FRM, Ed_unc_TARTU_CP_FRM, Ed_unc_HEREON_CP_FRM, 
                          diff_Ed_unc_PML,  diff_Ed_unc_NASA, diff_Ed_unc_TARTU, diff_Ed_unc_HEREON,
                          bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    
    rp.residuals_combined_V3('Lsky', Lsky_R_CP, Lsky_PML_CP,  Lsky_NASA_CP,  Lsky_TARTU_CP,  Lsky_HEREON_CP,
                             Lsky_R_CP_FRM, Lsky_PML_CP_FRM, Lsky_NASA_CP_FRM, Lsky_TARTU_CP_FRM, Lsky_HEREON_CP_FRM, 
                             Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, 
                             diff_Lsky_PML,  diff_Lsky_NASA, diff_Lsky_TARTU, diff_Lsky_HEREON,
                             Lsky_unc_PML_CP, Lsky_unc_NASA_CP, Lsky_unc_TARTU_CP, Lsky_unc_HEREON_CP, 
                             Lsky_unc_PML_CP_FRM, Lsky_unc_NASA_CP_FRM, Lsky_unc_TARTU_CP_FRM, Lsky_unc_HEREON_CP_FRM, 
                             diff_Lsky_unc_PML,   diff_Lsky_unc_NASA,  diff_Lsky_unc_TARTU ,  diff_Lsky_unc_HEREON,
                             bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
   
    rp.residuals_combined_V3('Lt', Lt_R_CP, Lt_PML_CP,  Lt_NASA_CP,  Lt_TARTU_CP,  Lt_HEREON_CP, 
                             Lt_R_CP_FRM, Lt_PML_CP_FRM, Lt_NASA_CP_FRM, Lt_TARTU_CP_FRM, Lt_HEREON_CP_FRM, 
                             Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, 
                             diff_Lt_PML,  diff_Lt_NASA, diff_Lt_TARTU, diff_Lt_HEREON,
                             Lt_unc_PML_CP, Lt_unc_NASA_CP, Lt_unc_TARTU_CP, Lt_unc_HEREON_CP, 
                             Lt_unc_PML_CP_FRM, Lt_unc_NASA_CP_FRM, Lt_unc_TARTU_CP_FRM, Lt_unc_HEREON_CP_FRM, 
                             diff_Lt_unc_PML,   diff_Lt_unc_NASA,  diff_Lt_unc_TARTU ,  diff_Lt_unc_HEREON,
                             bands, path_output, Q_mask, Qtype = 'QC_AOC_3')
    
    rp.residuals_combined_V3('Rrs',  Rrs_R_CP, Rrs_PML_CP, Rrs_NASA_CP,  Rrs_TARTU_CP,  Rrs_HEREON_CP, 
                             Rrs_R_CP_FRM, Rrs_PML_CP_FRM, Rrs_NASA_CP_FRM, Rrs_TARTU_CP_FRM, Rrs_HEREON_CP_FRM, 
                             Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON,
                             diff_Rrs_PML,  diff_Rrs_NASA, diff_Rrs_TARTU, diff_Rrs_HEREON,
                             Rrs_unc_PML_CP, Rrs_unc_NASA_CP, Rrs_unc_TARTU_CP, Rrs_unc_HEREON_CP, 
                             Rrs_unc_PML_CP_FRM, Rrs_unc_NASA_CP_FRM, Rrs_unc_TARTU_CP_FRM, Rrs_unc_HEREON_CP_FRM, 
                             diff_Rrs_unc_PML,   diff_Rrs_unc_NASA,  diff_Rrs_unc_TARTU,  diff_Rrs_unc_HEREON,
                             bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   


    #rp.residuals_combined_V2('Ed', Ed_R_CP, Ed_PML_CP, Ed_NASA_CP, Ed_TARTU_CP, Ed_HEREON_CP, 
     #                     Ed_PML_CP_FRM, Ed_NASA_CP_FRM, Ed_TARTU_CP_FRM, Ed_HEREON_CP_FRM,
      #                    Ed_unc_PML_CP, Ed_unc_NASA_CP, Ed_unc_TARTU_CP, Ed_unc_HEREON_CP, 
       #                   Ed_unc_PML_CP_FRM, Ed_unc_NASA_CP_FRM, Ed_unc_TARTU_CP_FRM, Ed_unc_HEREON_CP_FRM, 
        #                  diff_Ed_unc_PML,    diff_Ed_unc_NASA,    diff_Ed_unc_TARTU ,   diff_Ed_unc_HEREON,
         #                 bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
    
   # rp.residuals_combined_V2('Lsky', Lsky_R_CP, Lsky_PML_CP, Lsky_NASA_CP, Lsky_TARTU_CP, Lsky_HEREON_CP, 
    #                         Lsky_PML_CP_FRM, Lsky_NASA_CP_FRM, Lsky_TARTU_CP_FRM, Lsky_HEREON_CP_FRM, 
     #                        Lsky_unc_PML_CP, Lsky_unc_NASA_CP, Lsky_unc_TARTU_CP, Lsky_unc_HEREON_CP, 
      #                       Lsky_unc_PML_CP_FRM, Lsky_unc_NASA_CP_FRM, Lsky_unc_TARTU_CP_FRM, Lsky_unc_HEREON_CP_FRM, 
       #                      diff_Lsky_unc_PML,   diff_Lsky_unc_NASA,  diff_Lsky_unc_TARTU ,  diff_Lsky_unc_HEREON,
        #                     bands, path_output, Q_mask, Qtype = 'QC_AOC_3')   
   
  #  rp.residuals_combined_V2('Lt', Lt_R_CP, Lt_PML_CP, Lt_NASA_CP, Lt_TARTU_CP, Lt_HEREON_CP, 
   #                          Lt_PML_CP_FRM, Lt_NASA_CP_FRM, Lt_TARTU_CP_FRM, Lt_HEREON_CP_FRM, 
    #                         Lt_unc_PML_CP, Lt_unc_NASA_CP, Lt_unc_TARTU_CP, Lt_unc_HEREON_CP, 
     #                        Lt_unc_PML_CP_FRM, Lt_unc_NASA_CP_FRM, Lt_unc_TARTU_CP_FRM, Lt_unc_HEREON_CP_FRM, 
      #                       diff_Lt_unc_PML,   diff_Lt_unc_NASA,  diff_Lt_unc_TARTU ,  diff_Lt_unc_HEREON,
       #                      bands, path_output, Q_mask, Qtype = 'QC_AOC_3')
    
    
    #rp.residuals_combined_V2('Rrs', Rrs_R_CP, Rrs_PML_CP, Rrs_NASA_CP, Rrs_TARTU_CP, Rrs_HEREON_CP, 
     #                        Rrs_PML_CP_FRM, Rrs_NASA_CP_FRM, Rrs_TARTU_CP_FRM, Rrs_HEREON_CP_FRM,
      #                       Rrs_unc_PML_CP, Rrs_unc_NASA_CP, Rrs_unc_TARTU_CP, Rrs_unc_HEREON_CP, 
       #                      Rrs_unc_PML_CP_FRM, Rrs_unc_NASA_CP_FRM, Rrs_unc_TARTU_CP_FRM, Rrs_unc_HEREON_CP_FRM, 
        #                     diff_Rrs_unc_PML,   diff_Rrs_unc_NASA,  diff_Rrs_unc_TARTU,  diff_Rrs_unc_HEREON,
         #                    bands, path_output, Q_mask, Qtype = 'QC_AOC_3')  


    ############################################################################
    # Summary tables - IP data
    #############################################################################
    # Ed_table = rp.tabular_summary('Ed', Ed_R, Ed_PML, Ed_NASA, Ed_TARTU, Ed_HEREON, Ed_RBINS, bands, path_output, Q_mask, Qtype = 'QC_AOC_3')    
    # Lsky_table = rp.tabular_summary('Lsky', Lsky_R, Lsky_PML, Lsky_NASA, Lsky_TARTU, Lsky_HEREON, Lsky_RBINS, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # Lt_table = rp.tabular_summary('Lt', Lt_R, Lt_PML, Lt_NASA, Lt_TARTU, Lt_HEREON, Lt_RBINS, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # Rrs_table = rp.tabular_summary('Rrs', Rrs_R, Rrs_PML, Rrs_NASA, Rrs_TARTU, Rrs_HEREON, Rrs_RBINS, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
    # nLw = rp.tabular_summary('nLw', nLw_R, nLw_PML, nLw_NASA, nLw_TARTU, nLw_HEREON, nLw_RBINS, bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
  
    #filename  =  path_output +  '/'  + 'Ed_stations_NOAA.csv'
    #Ed_NOAA.to_csv(filename, na_rep ='NaN', index = False)
    
    #filename  =  path_output +  '/'  + 'Ed_stations_PML.csv'
    #Ed_PML.to_csv(filename, na_rep ='NaN', index = False)
    
    #filename  =  path_output +  '/'  + 'Ed_stations_Reference.csv'
    #Ed_R.to_csv(filename, na_rep ='NaN', index = False)
    
    # filename  =  path_output +  '/'  + 'nLw_stations_NOAA.csv'
    # nLw_NOAA.to_csv(filename, na_rep ='NaN', index = False)
    
    #filename  =  path_output +  '/'  + 'nLw_stations_SEAPRISM.csv'
    #nLw_SEAPRISM.to_csv(filename, na_rep ='NaN', index = False)

    #filename  =  path_output +  '/'  + 'QC_mask.csv'
    #Q_mask.to_csv(filename, na_rep ='NaN', index = False)



