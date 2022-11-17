#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 09:58:07 2022

Script that contains file reader functions for AAOT intercomparrision excercise
IP = `Individually Processed' (i.e. submitted by each team in Sep 2022)'

Th was previously in intercomparrison master

@author: tjor
"""


from csv import reader
import numpy as np
import pandas as pd
import datetime 
import os

# timestamp matching sub routine
def nearest(items, pivot): 
    'function to locate nearest values and indicies in list x relative to a pivot (used to locate nearest timestamps)'
    
    nearest_value = min(items, key=lambda x: abs(x - pivot)) # finds nearest value    
    nearest_index= items.index(min(items, key=lambda x: abs(x - pivot))) # finds nearest index
   
    return nearest_value, nearest_index

# PML read functions
def read_PML_data(path, bands):
    'Wrappers to read spectra in dataframe format'
    
    Ed = read_PML_spec(path, 'Es_mean', bands) 
    Lsky = read_PML_spec(path, 'Li_mean', bands) 
    Lt = read_PML_spec(path, 'Lt_mean', bands)
    Rrs = read_PML_spec(path, 'reflectance_mean', bands) 
    Rrs_std = read_PML_spec(path, 'reflectance_standard_deviation', bands) 
    nLw = read_PML_spec(path, 'nLw_mean', bands) 
    
    return Ed, Lsky, Lt, Rrs, Rrs_std, nLw

def read_PML_spec(path, S ,bands): 
    '''Returns dataframe for spectral type S, OLCI bands as columns, indexed by station
    ''''''PML has slightly different read function than other systems as we require all timestamps 
    (used as reference time)'''

    
    with open(path, 'r') as read_obj:
     
        csv_reader = reader(read_obj)
        
        station = np.arange(0,78)
        time_start = np.empty(78, dtype=object)
        windspeed = np.nan*np.ones(78)
        spec_data = np.nan*np.ones([78,19]) 
        azimuth = np.nan*np.ones(78)
    
        # fill up data matrix
        i = 0
        for row in csv_reader:
            if i == 70:
                station[i] = 70
                time_start[i] = '2022-07-20 14:00:00' # missing staion- TS added (from Tartu)
                windspeed[i] = 2.7
                azimuth[i] =  90
                i=i+1
            if(row[0]) == S: 
                    station[i] = row[1]
                    if station[i] != 2: # cleaning for duplicate
                        timestamp_i =  datetime.datetime.strptime(row[3], '%Y-%m-%d %H:%M:%S')
                        time_start[i] = str(timestamp_i)
                        windspeed[i] = row[4]
                        azimuth[i] = row[5]
                        spec_data[i,:] = row[6:25]
                    else:
                        time_start[i]= '2022-07-14 09:00:00' # duplicate station (from Tartu)
                        windspeed[i] = 1.2
                        azimuth[i] =  135
                    i = i + 1

        # convert to df format - columns are station, time, windspeed, spectra in each band
        df = pd.DataFrame(index = station) 
        df['time_start'] = time_start
        df['windspeed'] = windspeed
        df['azimuth'] = azimuth
        for i in range(len(bands)):
                df[str(bands[i])] = spec_data[:,i] 
       
    return df
    

# NASA read functions 
def read_NASA_data(path, bands):
    'Wrappers to read spectra in dataframe format'
    
    Ed = read_NASA_spec(path, 'ed_average', bands) 
    Lsky = read_NASA_spec(path, 'lsky_average', bands) 
    Lt = read_NASA_spec(path, 'lu_average', bands)
    Rrs = read_NASA_spec(path, 'reflectance_average', bands) 
    Rrs_std = read_NASA_spec(path, 'reflectance_uncertainty', bands)#  - unc not std?
    nLw = read_NASA_spec(path, 'nlw_average', bands) 

    return Ed, Lsky, Lt, Rrs, Rrs_std, nLw

def read_NASA_spec(path, S, bands): 
    'Returns dataframe for spectral type S, OLCI bands as columns, indexed by station'
    with open(path, 'r') as read_obj:
     
        csv_reader = reader(read_obj)  

        data = np.nan*np.ones([78,19]) 
        time_start = np.empty(78, dtype=object)
        windspeed = np.nan*np.ones(78)
        azimuth = np.nan*np.ones(78)
        station = np.arange(0,78)
       
        # fill up data matrix
        for row in csv_reader:
            if(row[0]) == S: 
                i = int(row[1]) 
                timestamp_i = datetime.datetime.strptime(row[3][0:19], '%Y-%m-%dT%H:%M:%S')
                time_start[i] = str(timestamp_i)
                windspeed[i] = row[4]
                azimuth[i] = row[5]
                data[i, :] = row[6:25]
             
         
        # convert to df format - columns are station, time, windspeed, spectra in each band
        df = pd.DataFrame(index = station) 
        df['time_start'] = time_start
        df['windspeed'] = windspeed
        df['azimuth'] = azimuth
        
        for i in range(len(bands)):
                if S =='ed_average' or S == 'lsky_average' or S == 'lu_average' or S == 'nlw_average':
                    df[str(bands[i])] = 10*data[:,i] 
                else:
                    df[str(bands[i])] = data[:,i] 
   
    return df


# Tartu read functions 
def read_TARTU_data(path, bands):
    'Wrappers to read spectra in dataframe format'
    
    Ed = read_TARTU_spec(path, 'ramses3_above_ed_average', bands) 
    Lsky = read_TARTU_spec(path, 'ramses3_above_lsky_average', bands)
    Lt = read_TARTU_spec(path, 'ramses3_above_lu_average', bands)
    Rrs = read_TARTU_spec(path, 'ramses3_above_reflectance_average', bands) 
    Rrs_std = read_TARTU_spec(path, 'ramses3_above_reflectance_standardDeviation', bands)
#    nLw = read_TARTU_spec(path, 'ramses3_above_nlw_average', bands) # yet to arrive
    
    return Ed, Lsky, Lt, Rrs, Rrs_std

def read_TARTU_spec(path, S, bands): 
    'Returns dataframe for spectral type S, OLCI bands as columns, indexed by station'
    'Timestamp is read in np64 format'
    
    with open(path, 'r') as read_obj:
     
        csv_reader = reader(read_obj)
        
        station = np.arange(0,78)
        time_start = np.empty(78,  dtype=object)
        windspeed = np.nan*np.ones(78)
        spec_data = np.nan*np.ones([78,19]) 
 
        # fill up data matrix
        for row in csv_reader:
            if(row[0]) == S: 
                i = int(row[3]) 
                timestamp_i = datetime.datetime.strptime(row[5][0:19], '%Y-%m-%dT%H:%M:%S')
                time_start[i] = str(timestamp_i)
                windspeed[i] = row[6]
                spec_data[i,:] = row[10:29]
   

        # convert to df format - columns are station, time, windspeed, spectra in each band
        df = pd.DataFrame(index = station) 
        df['time_start'] = time_start
        df['windspeed'] = windspeed
        for i in range(len(bands)):
            if S ==  'ramses3_above_reflectance_average':
                df[str(bands[i])] = spec_data[:,i]/np.pi
            elif  S ==  'ramses3_above_reflectance_standardDeviation':
                df[str(bands[i])] = spec_data[:,i]/np.sqrt(np.pi)
            else:
                df[str(bands[i])] = spec_data[:,i]
        
    return df


# RBINS read functions 
def read_RBINS_data(path, bands):
    'Wrappers to read spectra in dataframe format'
    
    Ed = read_RBINS_spec(path, 'ed_average', bands) 
    Lsky = read_RBINS_spec(path, 'lsky_average', bands)
    Lt = read_RBINS_spec(path, 'lu_average', bands)
    Rrs = read_RBINS_spec(path, 'reflectance_average', bands) 
    Rrs_std = read_RBINS_spec(path, 'reflectance_standardDeviation', bands) 
    nLw = read_RBINS_spec(path, 'nlw_average', bands) 
    
    return Ed, Lsky, Lt, Rrs, Rrs_std, nLw


def read_RBINS_spec(path, S, bands): 
    'Returns dataframe for spectral type S, OLCI bands as columns, indexed by station'
    with open(path, 'r') as read_obj:
     
        csv_reader = reader(read_obj)  

        spec_data = np.nan*np.ones([78,19]) 
        time_start = np.nan*np.ones(78,dtype=object)
        windspeed = np.nan*np.ones(78)
        azimuth = np.nan*np.ones(78)
        station = np.arange(0,78)
       
        # fill up data matrix
        for row in csv_reader:
            if(row[0]) == S: 
                if str(row[1][0]) != 'P': #
                    i = int(row[1][0:2])
                    timestamp_i = datetime.datetime.strptime(row[3][0:19], '%Y-%m-%dT%H:%M:%S')
                    if timestamp_i.hour < 4: # accounts for 12 hr clock
                       timestamp_i =  timestamp_i + datetime.timedelta(hours=12)
                    time_start[i] = str(timestamp_i)
                    windspeed[i] = row[4]
                    azimuth[i] = 135
                    spec_data[i,:] = row[6:25]

            
        # convert to df format - columns are station, time, windspeed, spectra in each band
        df = pd.DataFrame(index = station) 
        df['time_start'] = time_start
        df['windspeed'] = windspeed
        df['azimuth'] = azimuth
        for i in range(len(bands)):
            if S == 'reflectance_average':
                df[str(bands[i])] = spec_data[:,i]/np.pi
            elif  S ==  'ramses3_above_reflectance_standardDeviation':
                df[str(bands[i])] = spec_data[:,i]/np.sqrt(np.pi)
            else:
                df[str(bands[i])] = spec_data[:,i]
  
    return df


# HEREON read functions 
def read_HEREON_data(path, bands):
    'Wrappers to read spectra in dataframe format'
    
    Ed = read_HEREON_spec(path, 'ed_average', bands) 
    Lsky = read_HEREON_spec(path, 'lsky_average', bands)
    Lt = read_HEREON_spec(path, 'lt_average', bands)
    Rrs = read_HEREON_spec(path, 'reflectance_average', bands) 
    Rrs_std = read_HEREON_spec(path, 'reflectance_standardDeviation', bands) 
  #  nLw = read_HEREON_spec(path, 'lw_average', bands) 
    
    return Ed, Lsky, Lt, Rrs, Rrs_std


def read_HEREON_spec(path, S, bands): 
    ' read function for HEREON FILE FORMAT: First Line: Wavelength; 1) Rrs_median, 2) Rrs_std, 3) Rrs_average, 4) Ed_median, 5) Ed_std, 6) Ed_average, 7) Lu_median, 8) Lu_std, 9) Lu_average, 10) Lsky_median, 11) Lsky_std, 12) Lsky_average, 13) Lw_median, 14) Lw_std, 15) Lw_average'
   
    files = os.listdir(path)
    data = np.nan*np.ones([78,19]) 
    station = np.arange(0,78)
       
    if   S == 'ed_average':
        row_index = 6
    elif S == 'lt_average':
        row_index = 9
    elif S == 'lsky_average':
        row_index = 12
    elif S == 'reflectance_average':
        row_index = 3
    elif S == 'reflectance_standardDeviation':
            row_index = 2
    elif S == 'lw_average':
        row_index = 15

    for i in range(len(files)):
        station[i] = i
        text_i = np.loadtxt(path + '/' + files[i], skiprows = 1)
        data[i,:] = text_i[row_index,1:-1] 
        
        df = pd.DataFrame(index = station) 
        for i in range(len(bands)):
                df[str(bands[i])] = data[:,i] 
    
    return df

# NOAA read functions 
def read_NOAA_data(path, bands, Ed_PML):
    'Wrappers to read spectra in dataframe format - Ed_PML used for reference time'
    
    Ed = read_NOAA_spec(path, 'Es_average', bands) 
    Lsky = read_NOAA_spec(path, 'Lsky_average', bands) 
    Lt =  read_NOAA_spec(path, 'Lt_average', bands) 
    Rrs =  read_NOAA_spec(path, 'Rrs_average', bands) 
    nLw = read_NOAA_spec(path, 'nlw_average', bands)

    Ed = time_match_NOAA(Ed, Ed_PML, bands)
    nLw = time_match_NOAA(nLw, Ed_PML, bands)

    return Ed, Lsky, Lt, Rrs, nLw


def read_NOAA_spec(path, S ,bands): 
    'Returns dataframe for spectral type S, OLCI bands as columns, indexed by station'
    
    with open(path, 'r') as read_obj:
       
        # fill up data matrix  - case ed
        if  S == 'Es_average':
                
            csv_reader = reader(read_obj)  
            data = np.nan*np.ones([78,16]) 
            time_start = np.empty(78, dtype=object)
            # station = np.arange(0,78)
            i = 0
            for row in csv_reader:
                if(row[0]) == S: 
                    # print(row[6:22])
                    timestamp_i = datetime.datetime.strptime(row[3][0:19], '%Y-%m-%dT%H:%M:%S')
                    time_start[i] = str(timestamp_i)
                    data[i,:] = row[6:22]
                    i = i +1
                
            # convert to df format - 
            df = pd.DataFrame() 
            df['time_start'] = time_start
            
            for i in range(len(bands) - 3): # considers bands < 800 nm              
                df[str(bands[i])] = 10*data[:,i] 
           
       # fill up data matrix - case nLw or Lw
        elif S == 'nlw_average' or S == 'lw_average':
                        
            csv_reader = reader(read_obj)  
            data = np.nan*np.ones([78,8]) 
            time_start = np.empty(78, dtype=object)
            # station = np.arange(0,78)
            i = 0
            for row in csv_reader:
                if(row[0]) == S: 
                    print(row[6:14])
                    timestamp_i = datetime.datetime.strptime(row[3][0:19], '%Y-%m-%dT%H:%M:%S')
                    time_start[i] = str(timestamp_i)
                    data[i, :] = row[6:14]
                    i = i +1
                
            # convert to df format - 
            df = pd.DataFrame() 
            df['time_start'] = time_start
            
            for i in range(8): # considers bands < 800 nm              
                df[str(bands[i])] = 10*data[:,i] 
                
        else:
            df = pd.DataFrame() 
            
    return df

def time_match_NOAA(df_NOAA, Ed_PML, bands):
    'Sub routine to find station index for NOAA, and then reformat index of dataframe - assumes 10 min buffer'
    
    time_start_PML = [datetime.datetime.strptime(Ed_PML['time_start'][i],'%Y-%m-%d %H:%M:%S')  for i in range(len(Ed_PML))] 
    time_start_NOAA = [datetime.datetime.strptime(df_NOAA['time_start'][i],'%Y-%m-%d %H:%M:%S')  for i in range(72)]        # NOAA just has 72 timestamps     
    
    time_start_NOAA_matching = np.nan*np.ones(78,dtype = object)
    spec_data_matching = np.nan*np.ones([len(df_NOAA),len(df_NOAA.columns)]) 
    tol = 10*60       
    for i in range(len(time_start_PML)):
         nearest_time, nearest_index = nearest(time_start_NOAA, time_start_PML[i])
         delta_t = abs(time_start_PML[i] - nearest_time) 
         #   print(delta_t)
         #  print(nearest_index)
         if delta_t.total_seconds() < tol:
             time_start_NOAA_matching[i] = str(time_start_NOAA[nearest_index]) 
             for j in range(len(df_NOAA.columns)-1):
                 spec_data_matching[i,j] = df_NOAA[str(bands[j])][nearest_index]
      
    # convert to df format - columns are station, time, windspeed, spectra in each band
    df = pd.DataFrame() 
    df['time_start'] =  time_start_NOAA_matching
    for j in range(len(df_NOAA.columns)-1): # considers bands < 800 nm              
        df[str(bands[j])] =  spec_data_matching[:,j] 
           
    return df

# CNR read functions 
def read_CNR_data(path, bands):
    'Wrappers to read spectra in dataframe format'
    
    Ed = read_CNR_spec(path, 'ed_average', bands) 
    Lsky = read_CNR_spec(path, 'lsky_average', bands) 
    Lt = read_CNR_spec(path, 'lu_average', bands)
    Rrs = read_CNR_spec(path, 'reflectance_average', bands) 
    Rrs_std = read_CNR_spec(path, 'reflectance_stdev', bands)
    nlw = read_CNR_spec(path, 'nlw', bands) # - unc not std?

    return Ed, Lsky, Lt, Rrs, Rrs_std, nlw

def read_CNR_spec(path, S, bands): 
    'Returns dataframe for spectral type S, OLCI bands as columns, indexed by station'
    with open(path, 'r') as read_obj:
     
        csv_reader = reader(read_obj)  

        data = np.nan*np.ones([78,19]) 
        time_start = np.empty(78, dtype=object)
        azimuth = np.nan*np.ones(78)
        station = np.arange(0,78)
       
        # fill up data matrix
        for row in csv_reader:
            if(row[0]) == S: 
                station = row[1]
                if station != 'None':
                    i = int(row[1]) 
                    timestamp_i = datetime.datetime.strptime(row[3][0:19], '%Y-%m-%dT%H:%MZ')
                    time_start[i] = str(timestamp_i)
                    azimuth[i]= row[6]
                    data[i, :] = row[7:26]
             
         
        # convert to df format - columns are station, time, windspeed, spectra in each band
        df = pd.DataFrame() 
        df['time_start'] = time_start
        df['azimuth'] = azimuth
        if S == 'reflectance_average' or S == 'reflectance_stdev':
            for i in range(len(bands)):
                  df[str(bands[i])] = data[:,i]/np.pi 
        else:
            for i in range(len(bands)):
                  df[str(bands[i])] = data[:,i] 
                  
             
    return df

def read_Aeronet_nLw(path_NLW, Ed_R, bands):
    ''' function to read nLw 1.5 data - assumes similar format to baseline _R dataframes
    Based on nearest neighbour interpolation within a +/- 10 minute tolerance to be counted.
    Note : 667 nm band is currently interpreted as 665 nm OLCI band (may need interpolating later on)
    Units are mw/(cm^2 sr micro-m) - multiply by 10 to get  mw/(m^2 sr nm)'''
    
    # extract nLw data in data frame formate
    files = os.listdir(path_NLW)
    data = pd.DataFrame() 
    for i in range(len(files)): 
        data_i = pd.read_csv(path_NLW + '/' + files[i], skiprows = 0)
        data.append(data_i)
        data = pd.concat([data, data_i]) #
    data_AOC = data.reset_index(drop=True)
    
    # convert timestamps to datetime formats - Ed_R is used for station timestamps
    date_nLw = [datetime.datetime.strptime(str(data_AOC['Date(dd-mm-yyyy)'][i]),'%d:%m:%Y') for i in range(len(data_AOC))] # combine dates and times ino single timestamp
    time_nLw = [datetime.datetime.strptime(str(data_AOC['Time(hh:mm:ss)'][i]),'%H:%M:%S') for i in range(len(data_AOC))]
    timestamp_nLw = [datetime.datetime.combine(date_nLw[i].date(), time_nLw[i].time()) for i in range(len(data_AOC))] 
    timestamp_station = [datetime.datetime.strptime(Ed_R['time_start'][i],'%Y-%m-%d %H:%M:%S') + datetime.timedelta(0,2,30) for i in range(len(Ed_R))]       
         
    # fill up nLw data matrix based on time stamp matchung (similar approach to QC) - 
    keys_nLw = ['Lwn_f/Q[412nm]','Lwn_f/Q[443nm]','Lwn_f/Q[490nm]','Lwn_f/Q[510nm]','Lwn_f/Q[560nm]','Lwn_f/Q[620nm]','Lwn_f/Q[667nm]']
    nLw_data = np.nan*np.zeros([78, len(keys_nLw)])
    station = np.arange(0,78,1)
    tol = 10*60 # as with QC flags +/- 10 minutes is default tolerance (should approximately match QC)
    for i in range(len(timestamp_station)):# Loop over FICE timestamps to find good AOC flag with +/- tolerance
        nearest_time, nearest_index = nearest(timestamp_nLw, timestamp_station[i])
        delta_t = abs(timestamp_station[i] - nearest_time) # time difference
        if delta_t.total_seconds() < tol:
            for j in range(len(keys_nLw)):
                nLw_data[i,j] = 10*data_AOC[keys_nLw[j]][i] # factor 10 for unit conversion
                
    # fill up data frame with similar format as baseline references
    df_R = pd.DataFrame(index = station) 
    for j in range(len(keys_nLw)):
            df_R[str(bands[j+1])] = nLw_data[:,j]
    df_R['time_start']  = Ed_R['time_start'] 
    
    return df_R


def nLw_usingPMLBRDF(nLw_PML, Rrs_PML, Rrs_ins, institution, bands):
    ''' function to derrive TARTU & HEREON nLW using PML BRDF
    Based on re-arranging eq. 4 in Tiltstone 2020'''
    
  
    station = np.arange(0,78)
    nLw_ins = pd.DataFrame(index = station) 
    
    scale_factor = nLw_PML.iloc[:,3:23]/Rrs_PML.iloc[:,3:23] # == BRDF*F_O
    
    if institution == 'TARTU':
        nLw_ins =  scale_factor*Rrs_ins.iloc[:,2:21]
    elif institution == 'HEREON':
        nLw_ins =  scale_factor*Rrs_ins
    
    return nLw_ins, scale_factor
    

