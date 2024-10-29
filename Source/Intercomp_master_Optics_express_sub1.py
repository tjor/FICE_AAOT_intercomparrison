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
import results_plots_Optics_express as rp

 
if __name__ == '__main__':
    
    # path to data + output
    dir_data = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/DataSubmissions'
    #path_output = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/Output_paper/nLw_final_Z17_SimSpec/'
    path_output = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/Output_paper/EdLtLsktRrs_final_M99_noNIR/' # use for Zhang 17 condfig

    # team submissions  
    path_PML = dir_data + '/PML/FICE_submission_V4/FRM4SOC_2_FICE_22_AAOT_PML_HSAS_stationsummary_V4_rhoM99_nLw_localChl_noNIR_Earthsuncorrected.csv'
    path_NASA = dir_data  + '/NASA/FRM4SOC_2_FICE_22_NASA_pySAS_Sentinel3A_rev2.csv'    # NASA - no NIR/M99 - v1 in paper
    path_NASA2 = dir_data  + '/NASA/FRM4SOC_2_FICE_22_NASA_pySAS_Sentinel3A_rev3.csv'   # NASA - with NIR/Zhang17 - v2 in paper
    path_TARTU = dir_data  + '/UniTartu/averages.csv'  
    path_HEREON = dir_data + '/HEREON/FRM4SOC2_FICE22_HEREON_DATA_OLCI'
    path_HEREON_nLw = dir_data + '/HEREON/FRM4SOC2_HEREON_Normalization_revised/FRM4SOC2_FICE22_HEREON_DATA_Lwn_mean_SRF_OLCI.csv'
    path_HEREON_RAFT = dir_data + '/HEREON_RAFT/FRM4SOC2_FICE22_HEREON_RAFT_DATA_OLCI/'
    # path_NOAA = dir_data  + '/NOAA/NOAA_Hyperpro_sheet2.csv' # redundant version
    path_NOAA  = dir_data + '/NOAA/NOAA_230216_sheet2.csv'
    # path_CNR = dir_data + '/CNR/FRM4SOC_2_FICE_22_AAOT_CNR_HYPSTAR_w-rel-az_ALLDATA.csv' # redundate version
    path_CNR = dir_data + '/CNR/FRM4SOC_2_FICE_22_AAOT_CNR_HYPSTAR_w-rel-az_FILTERED_for_v_az.csv'
    path_RBINS = dir_data  + '/RBINS/FRM4SOC2_FICE_RBINS_2022-11-25_ViewCorrected_EarthSunCorrected.csv'
    
    # addtional data (QC + references) 
    path_QC = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/Aeronet_QC_mask/FRM4SOC-AAOT_V3_ope.txt'
    #path_NLW = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/nLw_Zibordireference/'    # L1.5
    path_NLW = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/nLw_Zibordireference_L2/' # L2
    # HCP output 
    dir_CP_class = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/CP/FICE22-Reprocessed_12-23/' # M99, no NIR
    dir_CP_FRM = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/CP/FICE22_v1.2.4_coserror_revised/M99_NoNIR/' # M99, no NIR - FRM - used in Ed, Lt, Lsky, Rrs plots
    #dir_CP_FRM = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/CP/FICE22_v1.2.4_coserror_revised/Z17_SimSpec/' # Z17, Sim spec.
    # dir_CP_FRM = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/CP/FICE22-Reprocessed_12-23/' # M99, no NIR
    # dir_CP_FRM = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/CP/FICE22_FRMbranch/' # M99. no NIR
    
    # OLCI bands #
    bands = [str(400), str(412.5), str(442.5),	str(490), str(510), str(560), str(620),	str(665), str(673.75), str(681.25), str(708.75), str(753.75), str(761.25), str(764.375), str(767.5), str(778.75), str(865), str(885), str(900)]
    
    # Read data Individually Processed (IP) data. Each team has own file reader in read_IP.data.py
    Ed_PML, Lsky_PML, Lt_PML, Rrs_PML, Rrs_std_PML, nLw_PML = rd_IP.read_PML_data(path_PML, bands)
    Ed_NASA, Lsky_NASA, Lt_NASA, Rrs_NASA, Rrs_std_NASA, nLw_NASA = rd_IP.read_NASA_data(path_NASA, bands)
    Ed_NASA2, Lsky_NASA2, Lt_NASA2, Rrs_NASA2, Rrs_std_NASA2, nLw_NASA2 = rd_IP.read_NASA_data(path_NASA2, bands)
    Ed_TARTU, Lsky_TARTU, Lt_TARTU, Rrs_TARTU, Rrs_std_TARTU = rd_IP.read_TARTU_data(path_TARTU, bands) # no nLw yet
    Ed_HEREON, Lsky_HEREON, Lt_HEREON, Rrs_HEREON, Rrs_std_HEREON = rd_IP.read_HEREON_data(path_HEREON, bands) # no nLw yet
    Ed_RBINS, Lsky_RBINS, Lt_RBINS, Rrs_RBINS, Rrs_std_RBINS, nLw_RBINS = rd_IP.read_RBINS_data(path_RBINS, bands)
    Ed_NOAA, Lsky_NOAA, Lt_NOAA, Rrs_NOAA, nLw_NOAA = rd_IP.read_NOAA_data_V2(path_NOAA, bands, Ed_PML) # PML timestamps used to reference staion no.
    Ed_CNR, Lsky_CNR, Lt_CNR, Rrs_CNR, Rrs_std_CNR, nLw_CNR = rd_IP.read_CNR_data(path_CNR, bands)
   
    
    nLw_TARTU, scale = rd_IP.nLw_usingPMLBRDF(nLw_PML, Rrs_PML, Rrs_TARTU,'TARTU', bands) # using PML f/q method and re-scaling
    # nLw_HEREON, scale = rd_IP.nLw_usingPMLBRDF(nLw_PML, Rrs_PML, Rrs_HEREON,'HEREON', bands) # using PML f/q method and re-scaling
    nLw_HEREON = pd.read_csv(path_HEREON_nLw,header=None,skiprows =1,usecols=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],delim_whitespace=True) # revised IP submission by HEREON
    nLw_HEREON.columns = bands
      

    # Read Community Procsssed data  - FRM/SS branch  Versions without nLw 
    class_based = False
    # Ed_PML_CP_FRM, Ed_unc_PML_CP_FRM, Lsky_PML_CP_FRM, Lsky_unc_PML_CP_FRM, Lt_PML_CP_FRM, Lt_unc_PML_CP_FRM,  Rrs_PML_CP_FRM, Rrs_unc_PML_CP_FRM = rd_CP.read_CP_data('PML', Ed_PML, dir_CP_FRM, bands, class_based) # Ed dataframe used for timestamp
    #Ed_TARTU_CP_FRM, Ed_unc_TARTU_CP_FRM, Lsky_TARTU_CP_FRM, Lsky_unc_TARTU_CP_FRM, Lt_TARTU_CP_FRM, Lt_unc_TARTU_CP_FRM,  Rrs_TARTU_CP_FRM, Rrs_unc_TARTU_CP_FRM = rd_CP.read_CP_data('TARTU', Ed_TARTU,  dir_CP_FRM, bands, class_based)
    # Ed_HEREON_CP_FRM, Ed_unc_HEREON_CP_FRM, Lsky_HEREON_CP_FRM, Lsky_unc_HEREON_CP_FRM, Lt_HEREON_CP_FRM, Lt_unc_HEREON_CP_FRM, Rrs_HEREON_CP_FRM, Rrs_unc_HEREON_CP_FRM = rd_CP.read_CP_data('HEREON', Ed_PML, dir_CP_FRM, bands, class_based) # PML timestamp used as HEREON did not include this (systems were close to synchronous)  
    #Ed_NASA_CP_FRM, Ed_unc_NASA_CP_FRM, Lsky_NASA_CP_FRM, Lsky_unc_NASA_CP_FRM, Lt_NASA_CP_FRM, Lt_unc_NASA_CP_FRM, Rrs_NASA_CP_FRM, Rrs_unc_NASA_CP_FRM = rd_CP.read_CP_data('NASA', Ed_NASA,  dir_CP_FRM, bands,  class_based) # 

    dir_CP_FRM = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/CP/FICE22_v1.2.4_coserror_revised/Z17_SimSpec/' # Load Sim-spec first for HyperPro plot - label nLw fields with explict Z17
    Ed_PML_CP_FRM, Ed_unc_PML_CP_FRM, Lsky_PML_CP_FRM, Lsky_unc_PML_CP_FRM, Lt_PML_CP_FRM, Lt_unc_PML_CP_FRM,  Rrs_PML_CP_FRM, Rrs_unc_PML_CP_FRM, nLw_PML_CP_FRM_Z17, nLw_unc_PML_CP_FRM_Z17 = rd_CP.read_CP_data_withnLw('PML', Ed_PML, dir_CP_FRM, bands, class_based) # Ed dataframe used for timestamp
    Ed_TARTU_CP_FRM, Ed_unc_TARTU_CP_FRM, Lsky_TARTU_CP_FRM, Lsky_unc_TARTU_CP_FRM, Lt_TARTU_CP_FRM, Lt_unc_TARTU_CP_FRM,  Rrs_TARTU_CP_FRM, Rrs_unc_TARTU_CP_FRM, nLw_TARTU_CP_FRM_Z17, nLw_unc_TARTU_CP_FRM_Z17 = rd_CP.read_CP_data_withnLw('TARTU', Ed_TARTU,  dir_CP_FRM, bands, class_based)
    Ed_HEREON_CP_FRM, Ed_unc_HEREON_CP_FRM, Lsky_HEREON_CP_FRM, Lsky_unc_HEREON_CP_FRM, Lt_HEREON_CP_FRM, Lt_unc_HEREON_CP_FRM, Rrs_HEREON_CP_FRM, Rrs_unc_HEREON_CP_FRM, nLw_HEREON_CP_FRM_Z17, nLw_unc_HEREON_CP_FRM_Z17 = rd_CP.read_CP_data_withnLw('HEREON', Ed_PML, dir_CP_FRM, bands, class_based) # PML timestamp used as HEREON did not include this (systems were close to synchronous)
    Ed_NASA_CP_FRM, Ed_unc_NASA_CP_FRM, Lsky_NASA_CP_FRM, Lsky_unc_NASA_CP_FRM, Lt_NASA_CP_FRM, Lt_unc_NASA_CP_FRM, Rrs_NASA_CP_FRM, Rrs_unc_NASA_CP_FRM, nLw_NASA_CP_FRM_Z17, nLw_unc_NASA_CP_FRM_Z17  = rd_CP.read_CP_data_withnLw('NASA', Ed_NASA,  dir_CP_FRM, bands,  class_based) # 

    dir_CP_FRM = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/CP/FICE22_v1.2.4_coserror_revised/M99_NoNIR/' # all fields apart from nLW_Z17 will now be overwritten
    Ed_PML_CP_FRM, Ed_unc_PML_CP_FRM, Lsky_PML_CP_FRM, Lsky_unc_PML_CP_FRM, Lt_PML_CP_FRM, Lt_unc_PML_CP_FRM,  Rrs_PML_CP_FRM, Rrs_unc_PML_CP_FRM, nLw_PML_CP_FRM, nLw_unc_PML_CP_FRM = rd_CP.read_CP_data_withnLw('PML', Ed_PML, dir_CP_FRM, bands, class_based) # Ed dataframe used for timestamp
    Ed_TARTU_CP_FRM, Ed_unc_TARTU_CP_FRM, Lsky_TARTU_CP_FRM, Lsky_unc_TARTU_CP_FRM, Lt_TARTU_CP_FRM, Lt_unc_TARTU_CP_FRM,  Rrs_TARTU_CP_FRM, Rrs_unc_TARTU_CP_FRM, nLw_TARTU_CP_FRM, nLw_unc_TARTU_CP_FRM = rd_CP.read_CP_data_withnLw('TARTU', Ed_TARTU,  dir_CP_FRM, bands, class_based)
    Ed_HEREON_CP_FRM, Ed_unc_HEREON_CP_FRM, Lsky_HEREON_CP_FRM, Lsky_unc_HEREON_CP_FRM, Lt_HEREON_CP_FRM, Lt_unc_HEREON_CP_FRM, Rrs_HEREON_CP_FRM, Rrs_unc_HEREON_CP_FRM, nLw_HEREON_CP_FRM, nLw_unc_HEREON_CP_FRM = rd_CP.read_CP_data_withnLw('HEREON', Ed_PML, dir_CP_FRM, bands, class_based) # PML timestamp used as HEREON did not include this (systems were close to synchronous)
    Ed_NASA_CP_FRM, Ed_unc_NASA_CP_FRM, Lsky_NASA_CP_FRM, Lsky_unc_NASA_CP_FRM, Lt_NASA_CP_FRM, Lt_unc_NASA_CP_FRM, Rrs_NASA_CP_FRM, Rrs_unc_NASA_CP_FRM, nLw_NASA_CP_FRM, nLw_unc_NASA_CP_FRM  = rd_CP.read_CP_data_withnLw('NASA', Ed_NASA,  dir_CP_FRM, bands,  class_based) # 

    
    # Read Community Procsssed data  - FRM/SS branch  Versions with nLw 
    Ed_PML_CP_FRM, Ed_unc_PML_CP_FRM, Lsky_PML_CP_FRM, Lsky_unc_PML_CP_FRM, Lt_PML_CP_FRM, Lt_unc_PML_CP_FRM,  Rrs_PML_CP_FRM, Rrs_unc_PML_CP_FRM, nLw_PML_CP_FRM, nLw_unc_PML_CP_FRM = rd_CP.read_CP_data_withnLw('PML', Ed_PML, dir_CP_FRM, bands, class_based) # Ed dataframe used for timestamp
    Ed_TARTU_CP_FRM, Ed_unc_TARTU_CP_FRM, Lsky_TARTU_CP_FRM, Lsky_unc_TARTU_CP_FRM, Lt_TARTU_CP_FRM, Lt_unc_TARTU_CP_FRM,  Rrs_TARTU_CP_FRM, Rrs_unc_TARTU_CP_FRM, nLw_TARTU_CP_FRM, nLw_unc_TARTU_CP_FRM = rd_CP.read_CP_data_withnLw('TARTU', Ed_TARTU,  dir_CP_FRM, bands, class_based)
    Ed_HEREON_CP_FRM, Ed_unc_HEREON_CP_FRM, Lsky_HEREON_CP_FRM, Lsky_unc_HEREON_CP_FRM, Lt_HEREON_CP_FRM, Lt_unc_HEREON_CP_FRM, Rrs_HEREON_CP_FRM, Rrs_unc_HEREON_CP_FRM, nLw_HEREON_CP_FRM, nLw_unc_HEREON_CP_FRM = rd_CP.read_CP_data_withnLw('HEREON', Ed_PML, dir_CP_FRM, bands, class_based) # PML timestamp used as HEREON did not include this (systems were close to synchronous)
    Ed_NASA_CP_FRM, Ed_unc_NASA_CP_FRM, Lsky_NASA_CP_FRM, Lsky_unc_NASA_CP_FRM, Lt_NASA_CP_FRM, Lt_unc_NASA_CP_FRM, Rrs_NASA_CP_FRM, Rrs_unc_NASA_CP_FRM, nLw_NASA_CP_FRM, nLw_unc_NASA_CP_FRM  = rd_CP.read_CP_data_withnLw('NASA', Ed_NASA,  dir_CP_FRM, bands,  class_based) # 


    # Class branch
    class_based = True
    Ed_PML_CP, Ed_unc_PML_CP, Lsky_PML_CP , Lsky_unc_PML_CP, Lt_PML_CP, Lt_unc_PML_CP,  Rrs_PML_CP, Rrs_unc_PML_CP = rd_CP.read_CP_data('PML', Ed_PML, dir_CP_class, bands, class_based) # Ed dataframe used for timestamp
    Ed_TARTU_CP, Ed_unc_TARTU_CP, Lsky_TARTU_CP, Lsky_unc_TARTU_CP, Lt_TARTU_CP, Lt_unc_TARTU_CP,  Rrs_TARTU_CP, Rrs_unc_TARTU_CP = rd_CP.read_CP_data('TARTU', Ed_TARTU,  dir_CP_class, bands, class_based)
    Ed_HEREON_CP, Ed_unc_HEREON_CP, Lsky_HEREON_CP, Lsky_unc_HEREON_CP, Lt_HEREON_CP, Lt_unc_HEREON_CP, Rrs_HEREON_CP, Rrs_unc_HEREON_CP = rd_CP.read_CP_data('HEREON', Ed_PML, dir_CP_class, bands, class_based) # PML timestamp used as HEREON did not include this (systems were close to synchronous)
    Ed_NASA_CP, Ed_unc_NASA_CP, Lsky_NASA_CP, Lsky_unc_NASA_CP, Lt_NASA_CP, Lt_unc_NASA_CP, Rrs_NASA_CP, Rrs_unc_NASA_CP = rd_CP.read_CP_data('NASA', Ed_NASA,  dir_CP_class, bands,  class_based) # 
    
    
    # filter NASA by aziumuth - PML system used as reference
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
    nLw_NASA_CP_FRM = QC.filter_by_azimuth_CP(Rrs_PML, Rrs_NASA, nLw_NASA_CP_FRM) 
    
    Lsky_unc_NASA_CP_FRM = QC.filter_by_azimuth_CP(Lsky_PML, Lsky_NASA, Lsky_unc_NASA_CP_FRM) # for CP just filter NASA. IP dataframe is used for azimuth
    Lt_unc_NASA_CP_FRM = QC.filter_by_azimuth_CP(Lt_PML, Lt_NASA, Lt_unc_NASA_CP_FRM)
    Rrs_unc_NASA_CP_FRM = QC.filter_by_azimuth_CP(Rrs_PML, Rrs_NASA, Rrs_unc_NASA_CP_FRM) 
    nLw_unc_NASA_CP_FRM_unc = QC.filter_by_azimuth_CP(Rrs_PML, Rrs_NASA, nLw_unc_NASA_CP_FRM) 
    
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
  
    # Reference baselines P (must be done post azimuth filtering)
    Ed_R_CP = QC.baseline_average_V2_CP('Ed',  Ed_PML,  Ed_PML_CP, Ed_NASA_CP, Ed_TARTU_CP, Ed_HEREON_CP) # baseline V2 is 4-way mean of PML, NASA, HEREON, TARTU
    Lsky_R_CP = QC.baseline_average_V2_CP('Lsky', Ed_PML, Lsky_PML_CP, Lsky_NASA_CP, Lsky_TARTU_CP, Lsky_HEREON_CP) # Ed_PML used for windspeed
    Lt_R_CP = QC.baseline_average_V2_CP('Lt', Ed_PML, Lt_PML_CP, Lt_NASA_CP, Lt_TARTU_CP, Lt_HEREON_CP)
    Rrs_R_CP = QC.baseline_average_V2_CP('Rrs', Ed_PML, Rrs_PML_CP, Rrs_NASA_CP, Rrs_TARTU_CP, Rrs_HEREON_CP)
    
    Ed_R_CP_FRM = QC.baseline_average_V2_CP('Ed',  Ed_PML,  Ed_PML_CP_FRM, Ed_NASA_CP_FRM, Ed_TARTU_CP_FRM, Ed_HEREON_CP_FRM)
    Lsky_R_CP_FRM = QC.baseline_average_V2_CP('Lsky', Lsky_PML, Lsky_PML_CP_FRM, Lsky_NASA_CP_FRM, Lsky_TARTU_CP_FRM, Lsky_HEREON_CP_FRM)
    Lt_R_CP_FRM = QC.baseline_average_V2_CP('Lt',  Lt_PML,  Lt_PML_CP_FRM, Lt_NASA_CP_FRM, Lt_TARTU_CP_FRM, Lt_HEREON_CP_FRM)
    Rrs_R_CP_FRM = QC.baseline_average_V2_CP('Rrs', Rrs_PML,  Rrs_PML_CP_FRM, Rrs_NASA_CP_FRM, Rrs_TARTU_CP_FRM, Rrs_HEREON_CP_FRM)
 

    
    path_NLW = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/nLw_Zibordireference_L2/' #
    nLw_SEAPRISM = rd_IP.read_Aeronet_nLw(path_NLW, Ed_R, bands)    
    nLw_RAFT =  rd_IP.read_HEREON_RAFT_data(path_HEREON_RAFT, bands)
    

    Q_mask = QC.QC_mask(path_QC, Ed_R, Ed_PML, Lsky_PML, Rrs_PML, Rrs_std_PML, path_output)   # L1.5 QC used in FICE report
     #Q_mask = QC.QC_mask_L2(path_QC, Ed_R, nLw_SEAPRISM,  path_output) # L2 QC - directly using SeapPRISM nlW passes 0- gives similar to above!


    # baseline plots
    class_based=True
    rp.residualbaseline_plot(Ed_R_CP, Ed_PML_CP, Ed_NASA_CP, Ed_TARTU_CP, Ed_HEREON_CP,
                              Lsky_R_CP, Lsky_PML_CP, Lsky_NASA_CP, Lsky_TARTU_CP, Lsky_HEREON_CP,
                              Lt_R_CP, Lt_PML_CP, Lt_NASA_CP, Lt_TARTU_CP, Lt_HEREON_CP,
                              Rrs_R_CP, Rrs_PML_CP, Rrs_NASA_CP, Rrs_TARTU_CP, Rrs_HEREON_CP,
                              Ed_unc_PML_CP, Ed_unc_NASA_CP, Ed_unc_TARTU_CP, Ed_unc_HEREON_CP,
                              Lsky_unc_PML_CP, Lsky_unc_NASA_CP, Lsky_unc_TARTU_CP, Lsky_unc_HEREON_CP,
                              Lt_unc_PML_CP, Lt_unc_NASA_CP, Lt_unc_TARTU_CP, Lt_unc_HEREON_CP,
                              Rrs_unc_PML_CP, Rrs_unc_NASA_CP, Rrs_unc_TARTU_CP, Rrs_unc_HEREON_CP,
                              bands, path_output, Q_mask, class_based=True, Qtype = 'QC_AOC_3')

   # class_based=True
    rp.residualbaseline_plot(Ed_R_CP_FRM, Ed_PML_CP_FRM, Ed_NASA_CP_FRM, Ed_TARTU_CP_FRM, Ed_HEREON_CP_FRM,
                             Lsky_R_CP_FRM, Lsky_PML_CP_FRM, Lsky_NASA_CP_FRM, Lsky_TARTU_CP_FRM, Lsky_HEREON_CP_FRM,
                              Lt_R_CP_FRM, Lt_PML_CP_FRM, Lt_NASA_CP_FRM, Lt_TARTU_CP_FRM, Lt_HEREON_CP_FRM,
                             Rrs_R_CP_FRM, Rrs_PML_CP_FRM, Rrs_NASA_CP_FRM, Rrs_TARTU_CP_FRM, Rrs_HEREON_CP_FRM,
                              Ed_unc_PML_CP_FRM, Ed_unc_NASA_CP_FRM, Ed_unc_TARTU_CP_FRM, Ed_unc_HEREON_CP_FRM,
                              Lsky_unc_PML_CP, Lsky_unc_NASA_CP_FRM, Lsky_unc_TARTU_CP_FRM, Lsky_unc_HEREON_CP_FRM,
                              Lt_unc_PML_CP_FRM, Lt_unc_NASA_CP_FRM, Lt_unc_TARTU_CP_FRM, Lt_unc_HEREON_CP_FRM,
                              Rrs_unc_PML_CP_FRM, Rrs_unc_NASA_CP_FRM, Rrs_unc_TARTU_CP_FRM, Rrs_unc_HEREON_CP_FRM,
                              bands, path_output, Q_mask, class_based=False, Qtype = 'QC_AOC_3')


    # Uncertainty plots
    class_based= True
    rp.unc_plot(Ed_unc_PML_CP, Ed_unc_NASA_CP, Ed_unc_TARTU_CP, Ed_unc_HEREON_CP,
                       Lsky_unc_PML_CP, Lsky_unc_NASA_CP, Lsky_unc_TARTU_CP, Lsky_unc_HEREON_CP,
                       Lt_unc_PML_CP, Lt_unc_NASA_CP, Lt_unc_TARTU_CP, Lt_unc_HEREON_CP,
                       Rrs_unc_PML_CP, Rrs_unc_NASA_CP, Rrs_unc_TARTU_CP, Rrs_unc_HEREON_CP,
                       bands, path_output, Q_mask, class_based, Qtype = 'QC_AOC_3')
     
    class_based= False
    rp.unc_plot(Ed_unc_PML_CP_FRM, Ed_unc_NASA_CP_FRM, Ed_unc_TARTU_CP_FRM, Ed_unc_HEREON_CP_FRM,
                        Lsky_unc_PML_CP_FRM, Lsky_unc_NASA_CP_FRM, Lsky_unc_TARTU_CP_FRM, Lsky_unc_HEREON_CP_FRM,
                       Lt_unc_PML_CP_FRM, Lt_unc_NASA_CP_FRM, Lt_unc_TARTU_CP_FRM, Lt_unc_HEREON_CP_FRM,
                        Rrs_unc_PML_CP_FRM, Rrs_unc_NASA_CP_FRM, Rrs_unc_TARTU_CP_FRM, Rrs_unc_HEREON_CP_FRM,
                        bands, path_output, Q_mask, class_based, Qtype = 'QC_AOC_3')
       
    # difference of Class and Sensor Specific
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

    rp.FRM_rad_diffplot(diff_Ed_PML, diff_Ed_NASA, diff_Ed_TARTU, diff_Ed_HEREON,  
                 diff_Lsky_PML, diff_Lsky_NASA, diff_Lsky_TARTU, diff_Lsky_HEREON, 
                 diff_Lt_PML, diff_Lt_NASA, diff_Lt_TARTU, diff_Lt_HEREON, 
                 diff_Rrs_PML, diff_Rrs_NASA, diff_Rrs_TARTU, diff_Rrs_HEREON, 
                 bands, path_output, Q_mask, Qtype = 'QC_AOC_3')


    # difference of Class and Sensor Specific uncertainty - in % terms (disregard?)
    # diff_Ed_unc_PML = 200*(Ed_unc_PML_CP.iloc[:,1:] - Ed_unc_PML_CP_FRM.iloc[:,1:])/(Ed_unc_PML_CP.iloc[:,1:] + Ed_unc_PML_CP_FRM.iloc[:,1:])
    # diff_Ed_unc_NASA = 200*(Ed_unc_NASA_CP.iloc[:,1:] - Ed_unc_NASA_CP_FRM.iloc[:,1:])/(Ed_unc_NASA_CP.iloc[:,1:] + Ed_unc_NASA_CP_FRM.iloc[:,1:])
    # diff_Ed_unc_TARTU = 200*(Ed_unc_TARTU_CP.iloc[:,1:] - Ed_unc_TARTU_CP_FRM.iloc[:,1:])/(Ed_unc_TARTU_CP.iloc[:,1:] + Ed_unc_TARTU_CP_FRM.iloc[:,1:])
    # diff_Ed_unc_HEREON = 200*(Ed_unc_HEREON_CP.iloc[:,1:] - Ed_unc_HEREON_CP_FRM.iloc[:,1:])/(Ed_unc_HEREON_CP.iloc[:,1:] + Ed_unc_HEREON_CP_FRM.iloc[:,1:])
    
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
    #diff_Rrs_unc_HEREON = 200*(Rrs_unc_HEREON_CP.iloc[:,1:] - Rrs_unc_HEREON_CP_FRM.iloc[:,1:])/(Rrs_unc_HEREON_CP.iloc[:,1:] + Rrs_unc_HEREON_CP_FRM.iloc[:,1:])


    #rp.FRM_uncdiffV2_plot(diff_Ed_unc_PML, diff_Ed_unc_NASA, diff_Ed_unc_TARTU, diff_Ed_unc_HEREON,  
     #                diff_Lsky_unc_PML, diff_Lsky_unc_NASA, diff_Lsky_unc_TARTU, diff_Lsky_unc_HEREON, 
      #               diff_Lt_unc_PML, diff_Lt_unc_NASA, diff_Lt_unc_TARTU, diff_Lt_unc_HEREON, 
       #              diff_Rrs_unc_PML, diff_Rrs_unc_NASA, diff_Rrs_unc_TARTU, diff_Rrs_unc_HEREON, 
        #             bands, path_output, Q_mask, Qtype = 'QC_AOC_3')

    # uncertainty differences - un-normalized
    diff_Ed_unc_PML = Ed_unc_PML_CP.iloc[:,1:] - Ed_unc_PML_CP_FRM.iloc[:,1:]
    diff_Ed_unc_NASA = Ed_unc_NASA_CP.iloc[:,1:] - Ed_unc_NASA_CP_FRM.iloc[:,1:]
    diff_Ed_unc_TARTU = Ed_unc_TARTU_CP.iloc[:,1:] - Ed_unc_TARTU_CP_FRM.iloc[:,1:]
    diff_Ed_unc_HEREON = Ed_unc_HEREON_CP.iloc[:,1:] - Ed_unc_HEREON_CP_FRM.iloc[:,1:]
    
    diff_Lsky_unc_PML = Lsky_unc_PML_CP.iloc[:,1:] - Lsky_unc_PML_CP_FRM.iloc[:,1:]
    diff_Lsky_unc_NASA = Lsky_unc_NASA_CP.iloc[:,1:] - Lsky_unc_NASA_CP_FRM.iloc[:,1:]
    diff_Lsky_unc_TARTU = Lsky_unc_TARTU_CP.iloc[:,1:] - Lsky_unc_TARTU_CP_FRM.iloc[:,1:]
    diff_Lsky_unc_HEREON = Lsky_unc_HEREON_CP.iloc[:,1:] - Lsky_unc_HEREON_CP_FRM.iloc[:,1:]
     
    diff_Lt_unc_PML = Lt_unc_PML_CP.iloc[:,1:] - Lt_unc_PML_CP_FRM.iloc[:,1:]
    diff_Lt_unc_NASA = Lt_unc_NASA_CP.iloc[:,1:] - Lt_unc_NASA_CP_FRM.iloc[:,1:]
    diff_Lt_unc_TARTU = Lt_unc_TARTU_CP.iloc[:,1:] - Lt_unc_TARTU_CP_FRM.iloc[:,1:]
    diff_Lt_unc_HEREON = Lt_unc_HEREON_CP.iloc[:,1:] - Lt_unc_HEREON_CP_FRM.iloc[:,1:]
    
    diff_Rrs_unc_PML = Rrs_unc_PML_CP.iloc[:,1:] - Rrs_unc_PML_CP_FRM.iloc[:,1:]
    diff_Rrs_unc_NASA = Rrs_unc_NASA_CP.iloc[:,1:] - Rrs_unc_NASA_CP_FRM.iloc[:,1:]
    diff_Rrs_unc_TARTU = Rrs_unc_TARTU_CP.iloc[:,1:] - Rrs_unc_TARTU_CP_FRM.iloc[:,1:]
    diff_Rrs_unc_HEREON = Rrs_unc_HEREON_CP.iloc[:,1:] - Rrs_unc_HEREON_CP_FRM.iloc[:,1:]


    rp.FRM_uncdiff_plot(diff_Ed_unc_PML, diff_Ed_unc_NASA, diff_Ed_unc_TARTU, diff_Ed_unc_HEREON,  
                      diff_Lsky_unc_PML, diff_Lsky_unc_NASA, diff_Lsky_unc_TARTU, diff_Lsky_unc_HEREON, 
                     diff_Lt_unc_PML, diff_Lt_unc_NASA, diff_Lt_unc_TARTU, diff_Lt_unc_HEREON, 
                    diff_Rrs_unc_PML, diff_Rrs_unc_NASA, diff_Rrs_unc_TARTU, diff_Rrs_unc_HEREON, 
         bands, path_output, Q_mask, Qtype = 'QC_AOC_3')



    # difference of Class and IP
    diff_Ed_PML = 100*(Ed_PML.iloc[:,3:] - Ed_PML_CP.iloc[:,1:])/(Ed_PML_CP.iloc[:,1:])
    diff_Ed_NASA = 100*(Ed_NASA.iloc[:,3:] - Ed_NASA_CP.iloc[:,1:])/(Ed_NASA_CP.iloc[:,1:])
    diff_Ed_TARTU = 100*(Ed_TARTU.iloc[:,2:] - Ed_TARTU_CP.iloc[:,1:])/(Ed_TARTU_CP.iloc[:,1:])
    diff_Ed_HEREON = 100*(Ed_HEREON.iloc[:,0:] - Ed_HEREON_CP.iloc[:,1:])/(Ed_HEREON_CP.iloc[:,1:])
   
    diff_Lsky_PML = 100*(Lsky_PML.iloc[:,3:] - Lsky_PML_CP.iloc[:,1:])/(Lsky_PML_CP.iloc[:,1:])
    diff_Lsky_NASA = 100*(Lsky_NASA.iloc[:,3:] - Lsky_NASA_CP.iloc[:,1:])/(Lsky_NASA_CP.iloc[:,1:])
    diff_Lsky_TARTU = 100*(Lsky_TARTU.iloc[:,2:] - Lsky_TARTU_CP.iloc[:,1:])/(Lsky_TARTU_CP.iloc[:,1:])
    diff_Lsky_HEREON = 100*(Lsky_HEREON.iloc[:,0:] - Lsky_HEREON_CP.iloc[:,1:])/(Lsky_HEREON_CP.iloc[:,1:])
    
    diff_Lt_PML = 100*(Lt_PML.iloc[:,3:] - Lt_PML_CP.iloc[:,1:])/(Lt_PML_CP.iloc[:,1:])
    diff_Lt_NASA = 100*(Lt_NASA.iloc[:,3:] - Lt_NASA_CP.iloc[:,1:])/(Lt_NASA_CP.iloc[:,1:])
    diff_Lt_TARTU = 100*(Lt_TARTU.iloc[:,2:] - Lt_TARTU_CP.iloc[:,1:])/(Lt_TARTU_CP.iloc[:,1:])
    diff_Lt_HEREON = 100*(Lt_HEREON.iloc[:,0:] - Lt_HEREON_CP.iloc[:,1:])/(Lt_HEREON_CP.iloc[:,1:])
   
    diff_Rrs_PML = 100*(Rrs_PML.iloc[:,3:] - Rrs_PML_CP.iloc[:,1:])/(Rrs_PML_CP.iloc[:,1:])
    diff_Rrs_NASA = 100*(Rrs_NASA.iloc[:,3:] - Rrs_NASA_CP.iloc[:,1:])/(Rrs_NASA_CP.iloc[:,1:])
    diff_Rrs_TARTU = 100*(Rrs_TARTU.iloc[:,2:] - Rrs_TARTU_CP.iloc[:,1:])/(Rrs_TARTU_CP.iloc[:,1:])
    diff_Rrs_HEREON = 100*(Rrs_HEREON.iloc[:,0:] - Rrs_HEREON_CP.iloc[:,1:])/(Rrs_HEREON_CP.iloc[:,1:])

    rp.FRM_rad_diffplot_IP(diff_Ed_PML, diff_Ed_NASA, diff_Ed_TARTU, diff_Ed_HEREON,  
                 diff_Lsky_PML, diff_Lsky_NASA, diff_Lsky_TARTU, diff_Lsky_HEREON, 
                 diff_Lt_PML, diff_Lt_NASA, diff_Lt_TARTU, diff_Lt_HEREON, 
                 diff_Rrs_PML, diff_Rrs_NASA, diff_Rrs_TARTU, diff_Rrs_HEREON, 
                 bands, path_output, Q_mask, Qtype = 'QC_AOC_3')


   # obsolete versions of nLw plot functions
   #  rp.residuals_combined_nlw_v2('nLw', nLw_SEAPRISM, nLw_NOAA, nLw_RAFT, 
    #                              nLw_PML_CP_FRM, nLw_NASA_CP_FRM, nLw_TARTU_CP_FRM, nLw_HEREON_CP_FRM, bands, 
   #                               path_output, Q_mask, Qtype = 'QC_AOC_3') 
  
  #   rp.unc_plot_nLw(nLw_SEAPRISM, nLw_NOAA, nLw_RAFT, 
  #                               nLw_unc_PML_CP_FRM, nLw_unc_NASA_CP_FRM, nLw_unc_TARTU_CP_FRM, nLw_unc_HEREON_CP_FRM,
   #                              bands, path_output, Q_mask, Qtype = 'QC_AOC_3')
                               
    # rp.residuals_unc_nlw('nLw',  nLw_SEAPRISM, nLw_NOAA,
     #                            nLw_PML_CP_FRM, nLw_NASA_CP_FRM, nLw_TARTU_CP_FRM, nLw_HEREON_CP_FRM, 
      #                           nLw_unc_PML_CP_FRM, nLw_unc_NASA_CP_FRM, nLw_unc_TARTU_CP_FRM, nLw_unc_HEREON_CP_FRM,
       #                          bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
       
    rp.residuals_unc_nlw_SeaPrism('nLw',  nLw_SEAPRISM, 
                                     nLw_PML_CP_FRM, nLw_NASA_CP_FRM, nLw_TARTU_CP_FRM, nLw_HEREON_CP_FRM, 
                                     nLw_unc_PML_CP_FRM, nLw_unc_NASA_CP_FRM, nLw_unc_TARTU_CP_FRM, nLw_unc_HEREON_CP_FRM,
                                     bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
      
    rp.residuals_unc_nlw_hyperpro('nLw', nLw_NOAA,
                                 nLw_PML_CP_FRM, nLw_NASA_CP_FRM, nLw_TARTU_CP_FRM, nLw_HEREON_CP_FRM, 
                                 nLw_unc_PML_CP_FRM, nLw_unc_NASA_CP_FRM, nLw_unc_TARTU_CP_FRM, nLw_unc_HEREON_CP_FRM,
                                 nLw_PML_CP_FRM_Z17, nLw_NASA_CP_FRM_Z17, nLw_TARTU_CP_FRM_Z17, nLw_HEREON_CP_FRM_Z17, 
                                 nLw_unc_PML_CP_FRM_Z17, nLw_unc_NASA_CP_FRM_Z17, nLw_unc_TARTU_CP_FRM_Z17, nLw_unc_HEREON_CP_FRM_Z17,
                                 bands, path_output, Q_mask, Qtype = 'QC_AOC_3') 
       