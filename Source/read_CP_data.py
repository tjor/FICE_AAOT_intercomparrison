#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 09:58:07 2022

Script that contains file reader functions for AAOT intercomparrision excercise
CP = `Community Processed' 

Includes spectral convolution to OLCI bands

@authors: aderu - ACRI (orginal file readers) 
          tjor - PML (modifications to compare with IP processed data)
         
"""

import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py
import datetime

from scipy.interpolate import interp1d

import netCDF4 as nc
import read_IP_data as rd_IP


def read_hdf(file, grp, var):
    ' function to read hdf file and save in pandas dataframe format'
    f = h5py.File(file, 'r')
    data = f[grp][var][()]
    
    return pd.DataFrame(data)

def read_L1B_hdf(file, grp, var):
    ' function to read L1B hdf file and save in pandas dataframe format'
    f = h5py.File(file, 'r')
    data = pd.DataFrame(f[grp][var][()])
    time = convert_datetime(data['Datetag'], data['Timetag2'])
    
    return time, data

def convert_datetime(date, time):
    ' converts cp time to date-time format'
    # 'DATETAG'  - Vector of length y, YYYYDOY.0, where DOY is the UTC sequential day of the year
    # 'TIMETAG2' - Vector of length y, HHMMSSUUU UTC (hour, minute, second, millisecond)
    time_str =  [str(x)+str(y) for x, y in zip(date,time)]
    time_str = time_str[0]
   
    return datetime.datetime.strptime(time_str, "%Y%j.0%H%M%S%f.0") 


def read_hyperspectral(institute, var):
    '''Reads hyperspectral data from each institute: specify institute (PML, NASA, TARTU, HEREON),
      and var: LI,LT, ES,nLw Rrs'''
      
    if "LI" in var:
        ref_name = "Li_mean"
        ref_std = "Li_standard_deviation"
        grp = "RADIANCE"
        grp_unc = "UNCERTAINTY_BUDGET"
        grp_std = "RADIANCE"    
        var_unc = "LI_unc"
        unit = "mW/sr.m2.nm"
        
    elif "LT" in var:
        ref_name = "Lt_mean"
        ref_std = "Lt_standard_deviation"
        grp = "RADIANCE"
        grp_unc = "UNCERTAINTY_BUDGET"
        grp_std = "RADIANCE"  
        var_unc = "LT_unc"
        unit = "mW/sr.m2.nm"
        
    elif "ES" in var:
        ref_name = "Es_mean"
        ref_std = "Es_standard_deviation"
        grp = "IRRADIANCE"
        grp_unc = "UNCERTAINTY_BUDGET"
        grp_std = "IRRADIANCE"  
        var_unc = "ES_unc"
        unit = "mW/m2.nm"
        
    elif "nLw" in var:
        ref_name = "nLw_mean"
        ref_std = "nLw_standard_deviation"
        grp = "REFLECTANCE"
        grp_unc = "REFLECTANCE"
        grp_std = None
        var_unc = "nLw_HYPER_unc"
        unit = "mW/m2.nm"
        
    elif "Rrs" in var:
        ref_name = "reflectance_mean"
        ref_std = "reflectance_standard_deviation"
        grp = "REFLECTANCE"
        grp_unc = "REFLECTANCE"
        var_unc = "Rrs_HYPER_unc"
        grp_std = None
        unit = "1/sr"
    
    # locate hdf files for each institute
    dir_data = "/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/CP/"  + institute 
    files = glob.glob(dir_data + "/*.hdf")
    
    # extract no. of timestamps (==stations) and no. of wavelengths
    root = h5py.File(files[0], 'r') # read 1st file to read wavelength
    data = pd.DataFrame(root[grp][var][()])
    wv = np.array([float(x) for x in data.columns[2:] ])
    N_wv = len(wv)
    N_stat = len(files)   
    
    if institute != 'NASA': # includes uncertainty
        # define empty data matrices
        cp_data = np.nan*np.ones([N_stat, N_wv])
        cp_unc = np.nan*np.ones([N_stat, N_wv])
        time = np.empty(N_stat, dtype=object)
        
        for i in range(len(files)):
            
            root = h5py.File(files[i], 'r')
            data = pd.DataFrame(root[grp][var][()]) # spectral data
            data_unc = pd.DataFrame(root[grp_unc][var_unc][()]) # uncertainty of spectral data 
            time[i] = convert_datetime(data['Datetag'], data['Timetag2'])
            root.close()
        
            # convert unit to align Seabird and trios in the same unit
            if "Rrs" not in var:
                scale_unit = 1e-3/1e-4 
                cp_data[i,:] = data.iloc[0][2:]
                cp_unc[i, :] = data_unc.iloc[0][:] 
            else:
                cp_data[i,:] = data.iloc[0][2:]
                cp_unc[i,:] = data_unc.iloc[0][2:] 
                
            if institute == 'PML' and "Rrs" not in var:
                cp_data[i,:] = 10*cp_data[i,:]
                cp_unc[i,:] = 10*cp_unc[i,:]
    
    else: # does not include uncertainty (just NASA)
        
        cp_data = np.nan*np.ones([N_stat, N_wv])
        time = np.empty(N_stat, dtype=object)
        cp_unc = None
        
        for i, f in enumerate(sorted(files)):
           
            root = h5py.File(f, 'r')
            data = pd.DataFrame(root[grp][var][()]) # spectral data
            time[i] = convert_datetime(data['Datetag'], data['Timetag2'])
            root.close()
        
            # convert unit to align Seabird and trios in the same unit
            if "Rrs" not in var:
                scale_unit = 1e-3/1e-4 
                cp_data[i,:] = scale_unit*data.iloc[0][2:]
            else:
                cp_data[i,:] = data.iloc[0][2:]
              
    return cp_data, cp_unc, time, wv

def unpack_srf(dir_srf):
    '''Function to unpack OLCI SRF in np array format'''
   
    srf_data = nc.Dataset(dir_srf + '/S3A_OL_SRF_20160713_mean_rsr.nc4') # load OLCI SRF
    srf_bands =  srf_data['nominal_centre_wavelength'][0:19]  # band centres up to 900 nm
    # band_width = srf_data['bandwidth_fwhm'][0:16] - FWHM of SRF - not needed
    srf_wv = np.array(srf_data['mean_spectral_response_function_wavelength'])
    srf = np.array(srf_data['mean_spectral_response_function'])

    return srf, srf_wv, srf_bands

def hyperspec_to_OLCIbands(S, time, wv, srf, srf_bands, srf_wv):
    ''' Function to downsample hyperspectral data to OLCI bands: input spectra 
   in dimenions of time * wavelength - this version of funciton returns
   data frame format'''
   
    D = np.nan*np.zeros([len(time),len(srf_bands)]) # matrix format for down-sampled spectra - time * spectral band

    # Loop for spectral down-sampling (outputs np array D)
    for j in range(len(S)): # loop over timestamps
        if np.sum(np.isnan(S[j,:])) < len(S[j,:]): # checks spectra exists ()
            wv_j = wv[np.isnan(S[j,:]) == 0] # remove wl bins with nan padding (rq. for interpolation functions)
            S_j = S[j, np.isnan(S[j,:]) == 0]    
            for k in range(len(srf_bands)): # loop over spectral bands
                first_wv = wv[int(np.where(np.isnan(S[j,:]) == 0)[0][0])]  # first and last wavelengths
                last_wv = wv[int(np.where(np.isnan(S[j,:]) == 0)[0][-1])]
                if np.min(srf_wv[k]) > first_wv and np.max(srf_wv[k]) < last_wv: # tests data exists in band
                    interp_funct = interp1d(wv_j, S_j, kind = 'cubic')  # interpolate to same wv interval as OLCI
                    S_j_interp = interp_funct(srf_wv[k]) # oversample on same wv range as OLCI SRF 
                    S_j_k = np.sum(S_j_interp*srf[k])/np.sum(srf[k]) # convolution of kth band for jth timestamp with OLCI SRF
                    D[j,k] =  S_j_k # fill data matrix
    
    # Loop for pandas dataframe format
    D_df = pd.DataFrame(index = time) # pandas data frame format for down-sampled spectra
    for k in range(len(srf_bands)):
           D_df[str(srf_bands[k])] = D[:,k] 
                
    return D_df


# timestamp matching sub routine
def nearest(items, pivot): 
    'function to locate nearest values and indicies in list x relative to a pivot (used to locate nearest timestamps)'
    
    nearest_value = min(items, key=lambda x: abs(x - pivot)) # finds nearest value    
    nearest_index= items.index(min(items, key=lambda x: abs(x - pivot))) # finds nearest index
   
    return nearest_value, nearest_index

def time_match(df_CP, df_IP, bands):
    'Sub routine to find station index for CP data based on IP dataframe & re-format the same as IP data'
    
    bands_CP = df_CP.keys() # slight difference in CP keys  - decimal point for integers
    
    time_start_CP = [datetime.datetime.strptime(str(df_CP.index[i])[0:19],'%Y-%m-%d %H:%M:%S')  for i in range(len(df_CP))]
    time_start_IP = np.empty(len(df_IP), dtype=object)
    for i in range(len(df_IP)):
        if df_IP['time_start'][i]  != None:
           time_start_IP[i] = datetime.datetime.strptime(df_IP['time_start'][i],'%Y-%m-%d %H:%M:%S')  

    # match CP to IP timestamps
    time_start_CP_matching = np.nan*np.ones(len(df_IP),dtype = object) # mathching timestamps
    spec_data_matching = np.nan*np.ones([len(df_IP),len(df_IP.columns)]) # matching spectra
    tol = 10*60   # 10 mins buffer     
    for i in range(len(time_start_IP)):
         if time_start_IP[i] != None:
             nearest_time, nearest_index = nearest(time_start_CP, time_start_IP[i])
             delta_t = abs(time_start_IP[i] - nearest_time) 
             # print(delta_t)
             # print(nearest_index)
             if delta_t.total_seconds() < tol:
                 time_start_CP_matching[i] = str(time_start_CP[nearest_index]) 
                 for j in range(len(df_CP.columns)-1):
                     spec_data_matching[i,j] = df_CP[str(bands_CP[j])][nearest_index]
      
    # re-convert to CP 
    df = pd.DataFrame() 
    df['time_start'] =  time_start_CP_matching
    for j in range(len(df_CP.columns)-1): # considers bands < 800 nm              
        df[str(bands[j])] = spec_data_matching[:,j] 
           
    return df


def time_match_allvars(Ed_IP, Ed_CP, Lsky_CP, Lt_CP, Rrs_CP, bands):
    ' Applies time matching to all variables. Ed_IP is used for timestamps'
    ' Function is called once for each institute'
        
    Ed_CP_matching = time_match(Ed_CP, Ed_IP, bands)
    Lsky_CP_matching = time_match(Lsky_CP, Ed_IP, bands)
    Lt_CP_matching = time_match(Lt_CP, Ed_IP, bands)
    Rrs_CP_matching = time_match(Rrs_CP, Ed_IP, bands)

    return  Ed_CP_matching,  Lsky_CP_matching,   Lt_CP_matching,   Rrs_CP_matching


def read_CP_data(institute, Ed_IP, bands):
    ' Applies time matching to all variables. Ed_IP is used for timestamps'
    ' Function is called once for each institute. NASA currently does not have uncertainty - this is stored as None'
    
     # Unpack SRF
    dir_srf = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/Source/HSAS/HSAS_Processing/FICE_Formatting' # contains OLCI SRF and template
    srf, srf_wv, srf_bands = unpack_srf(dir_srf)
    
    if institute == 'TARTU':
        institute = 'UT' # UT = TARTU for hdf file ID
    
    # read hyperspectral CP data
    Ed, Ed_unc, time, wv = read_hyperspectral(institute, 'ES_HYPER')
    Lsky, Lsky_unc, time, wv = read_hyperspectral(institute, 'LI_HYPER')    
    Lt, Lt_unc, time, wv = read_hyperspectral(institute, 'LT_HYPER')                                                                      
    Rrs, Rrs_unc, time, wv = read_hyperspectral(institute, 'Rrs_HYPER')
    
    # spectral convolution - also applies to uncertainty
    Ed = hyperspec_to_OLCIbands(Ed, time, wv, srf, srf_bands, srf_wv)  
    Lsky = hyperspec_to_OLCIbands(Lsky, time, wv, srf, srf_bands, srf_wv)  
    Lt = hyperspec_to_OLCIbands(Lt, time, wv, srf, srf_bands, srf_wv) 
    Rrs = hyperspec_to_OLCIbands(Rrs, time, wv, srf, srf_bands, srf_wv) 
    
    if institute != 'NASA':
        Ed_unc = hyperspec_to_OLCIbands(Ed_unc, time, wv, srf, srf_bands, srf_wv)  
        Lsky_unc = hyperspec_to_OLCIbands(Lsky_unc, time, wv, srf, srf_bands, srf_wv)  
        Lt_unc = hyperspec_to_OLCIbands(Lt_unc, time, wv, srf, srf_bands, srf_wv) 
        Rrs_unc = hyperspec_to_OLCIbands(Rrs_unc, time, wv, srf, srf_bands, srf_wv) 
    
    # Timestamp mathcing to IP-processed data
    Ed, Lsky, Lt, Rrs = time_match_allvars(Ed_IP, Ed, Lsky, Lt, Rrs, bands)
    if institute != 'NASA':
        Ed_unc, Lsky_unc, Lt_unc,  Rrs_unc = time_match_allvars(Ed_IP, Ed_unc, Lsky_unc, Lt_unc, Rrs_unc, bands)
          
    return  Ed, Ed_unc, Lsky, Lsky_unc, Lt, Lt_unc, Rrs, Rrs_unc



if __name__ == '__main__':
    
    # code in main was used to test CP file reading. When running via intercomp_master
    # `read_CP_data' is used to call data
    
    # unpack OLCI SRF
    dir_srf = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/Source/HSAS/HSAS_Processing/FICE_Formatting' # contains OLCI SRF and template
    srf, srf_wv, srf_bands = unpack_srf(dir_srf)


    # Read Hyperspectral data: PML - note: all L2 wavelengths + timestamps are the same
    Ed_PML, Ed_unc_PML, time_PML, wv_PML = read_hyperspectral('PML', 'ES_HYPER')
    Lsky_PML, Lsky_unc_PML, time_PML, wv_PML = read_hyperspectral('PML', 'LI_HYPER')
    Lt_PML, Lt_unc_PML, time_PML, wv_PML = read_hyperspectral('PML', 'LT_HYPER')                                                                      
    Rrs_PML, Rrs_unc_PML, time_PML, wv_PML = read_hyperspectral('PML', 'Rrs_HYPER')
    #nLw_PML, nLw_unc_PML, time_PML, wv_PML = read_hyperspectral('PML', 'nLw_HYPER')
    
    # Read Hyperspectral data: TARTU - note: all L2 wavelengths + timestamps are the same
    Ed_TARTU, Ed_unc_TARTU, time_TARTU, wv_TARTU = read_hyperspectral('UT', 'ES_HYPER')
    Lsky_TARTU, Lsky_unc_TARTU, time_TARTU, wv_TARTU = read_hyperspectral('UT', 'LI_HYPER')
    Lt_TARTU, Lt_unc_TARTU, time_TARTU, wv_TARTU = read_hyperspectral('UT', 'LT_HYPER')                                                                      
    Rrs_TARTU, Rrs_unc_TARTU, time_TARTU, wv_TARTU = read_hyperspectral('UT', 'Rrs_HYPER')
    # nLw_TARTU, nLw_unc_TARTU, time_TARTU, wv_TARTU = read_hyperspectral('UT', 'nLw_HYPER')
    
    # Read Hyperspectral data: HEREON - note: all L2 wavelengths + timestamps are the same
    Ed_HEREON, Ed_unc_HEREON, time_HEREON, wv_HEREON = read_hyperspectral('HEREON', 'ES_HYPER')
    Lsky_HEREON, Lsky_unc_HEREON, time_HEREON, wv_HEREON = read_hyperspectral('HEREON', 'LI_HYPER')
    Lt_HEREON, Lt_unc_HEREON, time_HEREON, wv_HEREON  = read_hyperspectral('HEREON', 'LT_HYPER')                                                                      
    Rrs_HEREON, Rrs_unc_HEREON, time_HEREON, wv_HEREONU = read_hyperspectral('HEREON', 'Rrs_HYPER')
    # nLw_HEREON, nLw_unc_HEREON, time_HEREON, wv_TARTU = read_hyperspectral('UT', 'nLw_HYPER')
    
    # Read Hyperspectral data: NASA - note: all L2 wavelengths + timestamps are the same
    Ed_NASA, Ed_unc_NASA, time_NASA, wv_NASA = read_hyperspectral('NASA', 'ES_HYPER')
    Lsky_NASA, Lsky_unc_NASA, time_NASA, wv_NASA = read_hyperspectral('NASA', 'LI_HYPER')
    Lt_NASA, Lt_unc_NASA, time_NASA, wv_NASA  = read_hyperspectral('NASA', 'LT_HYPER')                                                                      
    Rrs_NASA, Rrs_unc_NASA, time_NASA, wv_NASA = read_hyperspectral('NASA', 'Rrs_HYPER')
    # nLw_HEREON, nLw_unc_HEREON, time_HEREON, wv_TARTU = read_hyperspectral('UT', 'nLw_HYPER')


    # spectral convolution: PML (also applies to uncertainty )
    Ed_PML = hyperspec_to_OLCIbands(Ed_PML, time_PML, wv_PML, srf, srf_bands, srf_wv)  
    Lsky_PML = hyperspec_to_OLCIbands(Lsky_PML, time_PML, wv_PML, srf, srf_bands, srf_wv)  
    Lt_PML = hyperspec_to_OLCIbands(Lt_PML, time_PML, wv_PML, srf, srf_bands, srf_wv) 
    Rrs_PML = hyperspec_to_OLCIbands(Rrs_PML, time_PML, wv_PML, srf, srf_bands, srf_wv) 
    
    Ed_unc_PML = hyperspec_to_OLCIbands(Ed_unc_PML, time_PML, wv_PML, srf, srf_bands, srf_wv)  
    Lsky_unc_PML = hyperspec_to_OLCIbands(Lsky_unc_PML, time_PML, wv_PML, srf, srf_bands, srf_wv)  
    Lt_unc_PML = hyperspec_to_OLCIbands(Lt_unc_PML, time_PML, wv_PML, srf, srf_bands, srf_wv) 
    Rrs_unc_PML = hyperspec_to_OLCIbands(Rrs_unc_PML, time_PML, wv_PML, srf, srf_bands, srf_wv) 
    
    # spectral convolution TARTU (also applies to uncertainty )
    Ed_TARTU = hyperspec_to_OLCIbands(Ed_TARTU, time_TARTU, wv_TARTU, srf, srf_bands, srf_wv)  
    Lsky_TARTU = hyperspec_to_OLCIbands(Lsky_TARTU, time_TARTU, wv_TARTU, srf, srf_bands, srf_wv)  
    Lt_TARTU = hyperspec_to_OLCIbands(Lt_TARTU, time_TARTU, wv_TARTU, srf, srf_bands, srf_wv)  
    Rrs_TARTU = hyperspec_to_OLCIbands(Rrs_TARTU, time_TARTU, wv_TARTU, srf, srf_bands, srf_wv) 
    
    Ed_unc_TARTU = hyperspec_to_OLCIbands(Ed_unc_TARTU, time_TARTU, wv_TARTU, srf, srf_bands, srf_wv)  
    Lsky_unc_TARTU = hyperspec_to_OLCIbands(Lsky_unc_TARTU, time_TARTU, wv_TARTU, srf, srf_bands, srf_wv)  
    Lt_unc_TARTU = hyperspec_to_OLCIbands(Lt_unc_TARTU, time_TARTU, wv_TARTU, srf, srf_bands, srf_wv)  
    Rrs_unc_TARTU = hyperspec_to_OLCIbands(Rrs_unc_TARTU, time_TARTU, wv_TARTU, srf, srf_bands, srf_wv) 
      
    # spectral convolution: HEREON (also applies to uncertainty )
    Ed_HEREON = hyperspec_to_OLCIbands(Ed_HEREON, time_HEREON, wv_HEREON, srf, srf_bands, srf_wv)  
    Lsky_HEREON = hyperspec_to_OLCIbands(Lsky_HEREON, time_HEREON, wv_HEREON, srf, srf_bands, srf_wv)  
    Lt_HEREON = hyperspec_to_OLCIbands(Lt_HEREON, time_HEREON, wv_HEREON, srf, srf_bands, srf_wv)  
    Rrs_HEREON = hyperspec_to_OLCIbands(Rrs_HEREON, time_HEREON, wv_HEREON, srf, srf_bands, srf_wv)  
   
    Ed_unc_HEREON = hyperspec_to_OLCIbands(Ed_unc_HEREON, time_HEREON, wv_HEREON, srf, srf_bands, srf_wv)  
    Lsky_unc_HEREON = hyperspec_to_OLCIbands(Lsky_unc_HEREON, time_HEREON, wv_HEREON, srf, srf_bands, srf_wv)  
    Lt_unc_HEREON = hyperspec_to_OLCIbands(Lt_unc_HEREON, time_HEREON, wv_HEREON, srf, srf_bands, srf_wv)  
    Rrs_unc_HEREON = hyperspec_to_OLCIbands(Rrs_unc_HEREON, time_HEREON, wv_HEREON, srf, srf_bands, srf_wv)  
    
    # spectral convolution NASA (no uncertainty, yet )
    Ed_NASA = hyperspec_to_OLCIbands(Ed_NASA, time_NASA, wv_NASA, srf, srf_bands, srf_wv)  
    Lsky_NASA = hyperspec_to_OLCIbands(Lsky_NASA, time_NASA, wv_NASA, srf, srf_bands, srf_wv)  
    Lt_NASA = hyperspec_to_OLCIbands(Lt_NASA, time_NASA, wv_NASA, srf, srf_bands, srf_wv)  
    Rrs_NASA = hyperspec_to_OLCIbands(Rrs_NASA, time_NASA, wv_NASA, srf, srf_bands, srf_wv) 
    
    # Read IP data (used in timestamp matching)
   
    # path to data + output
    dir_data = '/data/datasets/cruise_data/active/FRM4SOC_2/FICE22/System_intercomparrison/DataSubmissions'

    # team submissions
    path_PML = dir_data + '/PML/FICE_submission_V2/FRM4SOC_2_FICE_22_AAOT_PML_HSAS_stationsummary_L1timesWindspeeds_phi_output.csv'
    path_NASA = dir_data + '/NASA/FRM4SOC_2_FICE_22_AAOT_pySAS_Sentinel3A_rev1.csv'
    path_TARTU = dir_data + '/UniTartu/averages.csv'  
    path_HEREON = dir_data + '/HEREON/FRM4SOC2_FICE22_HEREON_DATA_OLCI'
   
    # OLCI bands
    bands = [str(400), str(412.5), str(442.5),	str(490), str(510), str(560), str(620),	str(665), str(673.75), str(681.25), str(708.75), str(753.75), str(761.25), str(764.375), str(767.5), str(778.75), str(865), str(885), str(900)]
 
    Ed_PML_IP, Lsky_PML_IP, Lt_PML_IP, Rrs_PML_IP, Rrs_std_PML_IP, nLw_PML_IP = rd_IP.read_PML_data(path_PML, bands)
    Ed_NASA_IP, Lsky_NASA_IP, Lt_NASA_IP, Rrs_NASA_IP, Rrs_std_NASA_IP, nLw_NASA_IP = rd_IP.read_NASA_data(path_NASA, bands)
    Ed_TARTU_IP, Lsky_TARTU_IP, Lt_TARTU_IP, Rrs_TARTU_IP, Rrs_std_TARTU_IP, nLw_TARTU_IP = rd_IP.read_TARTU_data(path_TARTU, bands)
    Ed_HEREON_IP, Lsky_HEREON_IP, Lt_HEREON_IP, Rrs_HEREON_IP, Rrs_std_HEREON_IP, nLw_HEREON_IP = rd_IP.read_HEREON_data(path_HEREON, bands)


    # time-stamp matching
    Ed_PML, Lsky_PML,  Lt_PML, Rrs_PML = time_match_allvars(Ed_PML_IP, Ed_PML, Lsky_PML, Lt_PML, Rrs_PML, bands)
    Ed_unc_PML, Lsky_unc_PML, Lt_unc_PML,  Rrs_PML_unc = time_match_allvars(Ed_PML_IP, Ed_unc_PML, Lsky_unc_PML, Lt_unc_PML, Rrs_unc_PML, bands)
            
    Ed_TARTU, Lsky_TARTU,  Lt_TARTU,  Rrs_TARTU = time_match_allvars(Ed_TARTU_IP, Ed_TARTU, Lsky_TARTU, Lt_TARTU, Rrs_TARTU, bands)
    Ed_unc_TARTU, Lsky_unc_TARTU,  Lt_unc_TARTU,  Rrs_unc_TARTU = time_match_allvars(Ed_TARTU_IP, Ed_unc_TARTU, Lsky_unc_TARTU, Lt_unc_TARTU, Rrs_unc_TARTU, bands)
            
    Ed_HEREON, Lsky_HEREON,  Lt_HEREON,  Rrs_HEREON = time_match_allvars(Ed_PML_IP, Ed_HEREON, Lsky_HEREON, Lt_HEREON, Rrs_HEREON, bands) # PML timestamp used for HEREON (not available!) - these systems were approximately synchronous
    Ed_unc_HEREON, Lsky_unc_HEREON,  Lt_unc_HEREON,  Rrs_unc_HEREON = time_match_allvars(Ed_PML_IP, Ed_unc_HEREON, Lsky_unc_HEREON, Lt_unc_HEREON, Rrs_unc_HEREON, bands)  #PML timestamp used for HEREON (not available!) - these systems were approximately synchronous
    
    Ed_NASA, Lsky_NASA,  Lt_NASA,  Rrs_NASA = time_match_allvars(Ed_NASA_IP, Ed_NASA, Lsky_NASA, Lt_NASA, Rrs_NASA, bands)
    