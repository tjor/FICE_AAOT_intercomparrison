#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 10:19:31 2022

Data analysis scripts for FICe

- Baseline computations
- QC masks and azimuth filtering

@author: tjor
"""

import numpy as np
import pandas as pd
import datetime 

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# timestamp matching sub routine
def nearest(items, pivot): 
    'function to locate nearest values and indicies in list x relative to a pivot (used to locate nearest timestamps)'
    
    nearest_value = min(items, key=lambda x: abs(x - pivot)) # finds nearest value    
    nearest_index= items.index(min(items, key=lambda x: abs(x - pivot))) # finds nearest index
   
    return nearest_value, nearest_index

# baseline averages
def baseline_average(spec_type, df_PML, df_NASA, df_TARTU, df_HEREON):
    'Computes reference baseline (_R) based on weighthed mean of PML, NASA, HEREON, TARTU.'                             
    
    if spec_type == 'Ed':
        # first, average seabird reference systems
        df_seabird = ( 
            # combine dataframes into a single dataframe
            pd.concat([df_PML, df_NASA])
            # replace 0 values with nan to exclude them from mean calculation
            .replace(0, np.nan)
            .reset_index()
            # group by the row within the original dataframe
            .groupby("index")
            # calculate the mean
            .mean()
        )
        
    elif spec_type == 'Lsky' or  spec_type == 'Lt' or spec_type == 'Rrs':
        # first, average seabird reference systems
        df_seabird = df_PML
       
    # second, average trios reference systems
    df_trios = (
        # combine dataframes into a single dataframe
        pd.concat([df_TARTU, df_HEREON])
        # replace 0 values with nan to exclude them from mean calculation
        .replace(0, np.nan)
        .reset_index()
        # group by the row within the original dataframe
        .groupby("index")
        # calculate the mean
        .mean()
    )
    
    # average trios and seabird
    df_R = (
        # combine dataframes into a single dataframe
        pd.concat([df_seabird, df_trios])
        # replace 0 values with nan to exclude them from mean calculation
        .replace(0, np.nan)
        .reset_index()
        # group by the row within the original dataframe
        .groupby("index")
        # calculate the mean
        .mean()
    )
    
    # used to calculate number used in mean - require need 3 or more systems for mean to be defined
    df_N = pd.concat([df_PML, df_NASA, df_TARTU, df_HEREON])
    N_mean = np.zeros(78)
    for i in range(len(df_R)):
        no_of_records = 0
        for j in range(4):
            if np.isnan(df_N['400'][i].iloc[j]) == False:
                no_of_records =  no_of_records + 1
            N_mean[i] = no_of_records
            
    # PML used for timestamps and windspeeds (PML is complete dataset  so copied to reference)    
    df_R['time_start']  = df_PML['time_start'] 
    df_R['windspeed']  = df_PML['windspeed']     
    df_R['N_mean']  = N_mean    
    
    return df_R


def baseline_average_V2(spec_type, df_PML, df_NASA, df_TARTU, df_HEREON):
      'Computes reference baselines - assummes even weighting of each submitted team'''
    
      df_R = ( 
            # combine dataframes into a single dataframe
            pd.concat([df_PML, df_NASA, df_TARTU, df_HEREON])
            # replace 0 values with nan to exclude them from mean calculation
            .replace(0, np.nan)
            .reset_index()
            # group by the row within the original dataframe
            .groupby("index")
            # calculate the mean
            .mean()
        )
       
       # used to calculate number used in mean - require need 3 or more systems for mean to be defined
      df_N =  pd.concat([df_PML, df_NASA, df_TARTU, df_HEREON])
      N_mean = np.zeros(78)
      for i in range(len(df_R)):
           no_of_records = 0
           for j in range(4):
               if np.isnan(df_N['400'][i].iloc[j]) == False:
                   no_of_records =  no_of_records + 1
               N_mean[i] = no_of_records
               
       # PML used for timestamps and windspeeds (PML is complete dataset  so copied to reference)    
      df_R['time_start']  = df_PML['time_start'] 
      df_R['windspeed']  = df_PML['windspeed']     
      df_R['N_mean']  = N_mean    
         
      return df_R
  
def baseline_average_V2_CP(spec_type, Ed_PML,df_PML, df_NASA, df_TARTU, df_HEREON):
          'Computes reference baselines - assummes even weighting of each submitted team'''
        
          df_R = ( 
                # combine dataframes into a single dataframe
                pd.concat([df_PML, df_NASA, df_TARTU, df_HEREON])
                # replace 0 values with nan to exclude them from mean calculation
                .replace(0, np.nan)
                .reset_index()
                # group by the row within the original dataframe
                .groupby("index")
                # calculate the mean
                .mean()
            )
           
           # used to calculate number used in mean - require need 3 or more systems for mean to be defined
          df_N =  pd.concat([df_PML, df_NASA, df_TARTU, df_HEREON])
          N_mean = np.zeros(78)
          for i in range(len(df_R)):
               no_of_records = 0
               for j in range(4):
                   if np.isnan(df_N['400'][i].iloc[j]) == False:
                       no_of_records =  no_of_records + 1
                   N_mean[i] = no_of_records
                   
           # Ed PML used for timestamps and windspeeds (PML is complete dataset  so copied to reference)    
          df_R['time_start']  = Ed_PML['time_start'] 
          df_R['windspeed']  = Ed_PML['windspeed']     
          df_R['N_mean']  = N_mean    
             
          return df_R
  
    
def QC_mask(path_QC, Ed_R, Ed_PML, Lsky_PML, Rrs_PML, Rrs_std_PML, path_output):
    
    ''' Function to derive QC mask considering: (i) Aeronet QC, (ii) Filtering based
    on cloudiness and CV_rrs (averaged on OLCI bands 1-6). 1==good data, 0==bad data
   
    - Method (i) icludes centering data timestamp to middle of 5 min recording interval.
    The threshold rank_AOC = AOC['Rank'] # rank 0.6 or 1 is a QC pass
   
    - Ed_R, Lsky_R, Rrs_R, Rrs_std_R  (4-way reference values) are used in filtering step 
    for method (ii).
    
    - tol is the maximum time (+/-) that a good Aeronet OC timestamp can be from a FICE 
    data record for the fICE
    

    '''   
    
    # (i) QC using Aeronet QC flags from Zibordiplt.figure(figsize=(12,12))
    
    # Read AOC data and convert FICE timestamp to dt format
    AOC_data = pd.read_csv(path_QC, delim_whitespace=True, skiprows=[0])   # QC mask from aeronet
    timestamp_AOC = [datetime.datetime.strptime(AOC_data['Lev_1.5_record_ID'][i][5:24],'%d:%m:%Y-%H:%M:%S') for i in range(len(AOC_data))]  
    rank_AOC = AOC_data['Rank'] 
    timestamp_data = [datetime.datetime.strptime(Ed_R['time_start'][i],'%Y-%m-%d %H:%M:%S') + datetime.timedelta(0,2,30) for i in range(len(Ed_R))] 
    
    # Loop over FICE timestamps to find good AOC flag with +/- tolerance
    tol = 10*60 # +/- tolerance in seconds - +/1 10 minutes used as default
    QC_AOC = np.zeros(78) # AOC_QC mask: zero == failed QC, one = passed QC
    QC_AOC_3 = np.zeros(78)  # AOC_QC with at least 3 systems present
    for i in range(len(timestamp_data)):
      nearest_time, nearest_index = nearest(timestamp_AOC, timestamp_data[i]) # finds nearest AOC timestamp and index to FICE record
      delta_t = abs(timestamp_data[i] - nearest_time) # time difference
      if delta_t.total_seconds() < tol and rank_AOC[nearest_index] >= 0.6: # test if time difference and rank criteria are satisfied
          QC_AOC[i] = 1 # 
      if delta_t.total_seconds() < tol and rank_AOC[nearest_index] >=0.6 and Ed_R['N_mean'][i] > 2:
          QC_AOC_3[i] = 1 # 
    # (ii) Cloudiness index QC  
    
    CI = np.pi*Lsky_PML['400']/Ed_PML['400'] # cloudiness index from Simis 2013. Less than  0.4 is clear or light haze
    CV_Rrs_specav = (1/6)*100*(Rrs_std_PML['400']/Rrs_PML['400'] + Rrs_std_PML['412.5']/Rrs_PML['412.5']+ Rrs_std_PML['442.5']/Rrs_PML['442.5']+ Rrs_std_PML['490']/Rrs_PML['490']+ Rrs_std_PML['510']/Rrs_PML['510'] + Rrs_std_PML['560']/Rrs_PML['560'] + Rrs_std_PML['560']/Rrs_PML['560']) # baverage on bands 1-6
    
    #
    plt.figure(figsize=(10, 8))
    plt.rc('font', size=16)   
    plt.title('Dependence of Rrs variability on cloudiness (station average values)')
     # colors = cm.jet(np.linspace(0, 1, len(df_overall)))
    #for i in range(len(df_overall)):
    plt.scatter(CV_Rrs_specav,CI,color='blue')
    #  plt.legend(labels = df_overall['station'], ncol=2,fontsize=18)
    plt.xlabel('CV[$R_{rs}$]: mean for OLCI bands 1-6 [%]')
    plt.ylabel('$\pi$$L_{sky}$(400)/$E_{d}$(400): (cloudiness index)')
    plt.ylim(0,1)
    plt.xlim(0,10)

    filename  =  path_output +  '/' + 'Method2QC.png'
    plt.savefig(filename)


    # Loop over FICE timestamps to find good condtions flag
    QC_CI = np.zeros(78) # mask that passes CI QC
    QC_CI_3 = np.zeros(78) # mask that passes CI QC and at least 3 systems present
    for i in range(len(timestamp_data)):
        if CI[i] < 0.4 and CV_Rrs_specav[i] < 2.5:
            QC_CI[i] = 1
        if CI[i] < 0.4 and CV_Rrs_specav[i] < 2.5 and  Ed_R['N_mean'][i] > 2:
            QC_CI_3[i] = 1
            
    
    #  fill up dataframe with QC fields present 
    station = np.array(Ed_R.index)
    df = pd.DataFrame(index = station)  
   # df['time_start'] =  Ed_R['time_start']
    #df['rank_AOC'] = rank_AOC
    #df['QC_AOC'] = QC_AOC
    df['QC_AOC_3'] = QC_AOC_3 # default AOC QC
    #df['Cloudiness_Index'] = CI
    #df['windspeed'] =  Ed_R['windspeed'] 
    #df['CV_Rrs_specAv'] = CV_Rrs_specav
    #df['QC_CI'] = QC_CI # 
    df['QC_CI_3'] = QC_CI_3 # default C.I QC
    
    
    filename  =  path_output +  '/' + 'AAOT_QCmask.csv'
    df.to_csv(filename, na_rep ='NaN')
    
    return df


def plot_dataandQCmasks(Q_mask, Ed_PML, Ed_TARTU, Ed_HEREON, Ed_NASA, Ed_RBINS, Ed_CNR,Ed_NOAA, path_output):
   
    '''plot to show data and QC masks for AAOT: (i) Data no QC, (ii) QC, 
    (iii) Data with AOC QC, (vi) Data with cloudiness index w'''
    # 
    
    
    # 1. plot of all data records
    # creates plotting mask ,Z,  where records exist
    Z = np.zeros([7, 78])
    for i in range(78):
        for j in range(78):
            if Ed_NOAA.index[j] == i and ~np.isnan(np.array(Ed_NOAA['400'])[j]) == 1:
                Z[0,i] = 1
            if Ed_CNR.index[j] == i and ~np.isnan(np.array(Ed_CNR['400'])[j]) == 1: # replace with CNR
                Z[1,i] = 1
            if Ed_RBINS.index[j] == i and ~np.isnan(np.array(Ed_RBINS['400'])[j]) == 1:
                Z[2,i] = 1
            if Ed_NASA.index[j] == i and ~np.isnan(np.array(Ed_NASA['400'])[j]) == 1:
                Z[3,i] = 1
            if Ed_HEREON.index[j] == i and ~np.isnan(np.array(Ed_HEREON['400'])[j]) == 1:
                Z[4,i] = 1
            if Ed_TARTU.index[j] == i and ~np.isnan(np.array(Ed_TARTU['400'])[j]) == 1:
                Z[5,i] = 1
            if Ed_PML.index[j] == i and ~np.isnan(np.array(Ed_PML['400'])[j]) == 1:
                Z[6,i] = 1
    
    # figure for stations
    plt.figure(figsize=(12,8))
    plt.rc('font', size=18)       
    cmap1 = matplotlib.cm.get_cmap('coolwarm')
    plt.pcolor(Z,edgecolors='k',cmap=cmap1)
    plt.title('AAOT station mask')
    
    plt.xlabel('Station number')
    ylabels = ['NOAA', 'CNR','RBINS','NASA','HEREON', 'TARTU','PML']
    plt.gca().set_yticklabels(ylabels)
    plt.gca().set_yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5])
    
    red_patch = mpatches.Patch(color=cmap1(1.001), label='Data')
    blue_patch = mpatches.Patch(color=cmap1(0), label='No data')
    plt.legend(handles=[red_patch,blue_patch],loc=3)
    filename  =  path_output +  '/' + 'stationmask.png'
    plt.savefig(filename)
    
    # 2. plot of QC
    # creates plotting mask ,Z,  where records exist
    Z = np.zeros([2, 78])
    Z[0,:] = Q_mask['QC_CI_3']
    Z[1,:] = Q_mask['QC_AOC_3']
  
    # figure for stations
    plt.figure(figsize=(14,6))
    plt.rc('font', size=18)      
    plt.pcolor(Z,edgecolors='k',cmap='coolwarm')
    plt.title('Quality control mask')
    ylabels = ['2. Env.' + '\n' + '& CV[$R_{rs}$]', '1. AOC:'  +'\n' +  'SEAPRISM']
    plt.xlabel('Station number')
    plt.gca().set_yticklabels(ylabels)
    plt.gca().set_yticks([0.5,1.5])
    
    red_patch = mpatches.Patch(color=cmap1(1.001), label='Pass QC')
    blue_patch = mpatches.Patch(color=cmap1(0), label='Fail QC')
    plt.legend(handles=[red_patch,blue_patch],loc=4)  
    
    filename  =  path_output +  '/' + 'QCmask.png'
    plt.savefig(filename,dpi=300)
    
    # 3. plot of data wwith AOC QC
    
    # creates plotting mask ,Z, where records exist
    Z = np.zeros([7, 78])
    for i in range(78):
        for j in range(78):
            if Ed_NOAA.index[j] == i and ~np.isnan(np.array(Ed_NOAA['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[0,i] = 1
            if Ed_CNR.index[j] == i and ~np.isnan(np.array(Ed_CNR['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[1,i] = 1    
            if Ed_RBINS.index[j] == i and ~np.isnan(np.array(Ed_RBINS['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[2,i] = 1
            if Ed_NASA.index[j] == i and ~np.isnan(np.array(Ed_NASA['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[3,i] = 1
            if Ed_HEREON.index[j] == i and ~np.isnan(np.array(Ed_HEREON['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[4,i] = 1
            if Ed_TARTU.index[j] == i and ~np.isnan(np.array(Ed_TARTU['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[5,i] = 1
            if Ed_PML.index[j] == i and ~np.isnan(np.array(Ed_PML['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[6,i] = 1
    
    # figure for stations
    plt.figure(figsize=(12,8))
    plt.rc('font', size=18)         
    plt.pcolor(Z,edgecolors='k',cmap='coolwarm')
    plt.title('AAOT station mask:' + '\n' +  'AOC-SEAPRISM QC and $\geq$ 3 reference systems with data')
    plt.xlabel('Station number')
    ylabels = ['NOAA','CNR', 'RBINS','NASA','HEREON', 'TARTU','PML']
    plt.gca().set_yticklabels(ylabels)
    plt.gca().set_yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5])
    
    red_patch = mpatches.Patch(color=cmap1(1.001), label='Data')
    blue_patch = mpatches.Patch(color=cmap1(0), label='No data')
    plt.legend(handles=[red_patch,blue_patch],loc=3)
    
    filename  =  path_output +  '/' + 'stationmask_QC1.png'
    plt.savefig(filename)
    
    # 3. plot of data with Cloudiness Index QC

    # creates plotting mask ,Z,  where records exist
    Z = np.zeros([7, 78])
    for i in range(78):
        for j in range(78):
            if Ed_NOAA.index[j] == i and ~np.isnan(np.array(Ed_NOAA['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[0,i] = 1
            if Ed_CNR.index[j] == i and ~np.isnan(np.array(Ed_CNR['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[1,i] = 1
            if Ed_RBINS.index[j] == i and ~np.isnan(np.array(Ed_RBINS['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[2,i] = 1
            if Ed_NASA.index[j] == i and ~np.isnan(np.array(Ed_NASA['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[3,i] = 1
            if Ed_HEREON.index[j] == i and ~np.isnan(np.array(Ed_HEREON['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[4,i] = 1
            if Ed_TARTU.index[j] == i and ~np.isnan(np.array(Ed_TARTU['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[5,i] = 1
            if Ed_PML.index[j] == i and ~np.isnan(np.array(Ed_PML['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[6,i] = 1
    
    # figure for stations
    plt.figure(figsize=(12,8))
    plt.rc('font', size=18)         
    plt.pcolor(Z,edgecolors='k',cmap='coolwarm')
    plt.title('AAOT station mask:' + '\n' + 'C.I and CV[$R_{rs}$] QC and $\geq$ 3 reference systems with data')
    
    plt.xlabel('Station number')
    ylabels = ['NOAA','CNR', 'RBINS','NASA','HEREON', 'TARTU','PML']
    plt.gca().set_yticklabels(ylabels)
    plt.gca().set_yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5])
    
    red_patch = mpatches.Patch(color=cmap1(1.001), label='Data')
    blue_patch = mpatches.Patch(color=cmap1(0), label='No data')
    plt.legend(handles=[red_patch,blue_patch],loc=3)
    
    filename  =  path_output +  '/' + 'stationmask_QC2.png'
    plt.savefig(filename)
    
    return

def plot_dataandQCmasks_4teams(Q_mask, Ed_PML, Ed_TARTU, Ed_HEREON, Ed_NASA, path_output):
   
    '''plot to show data and QC masks for AAOT: (i) Data no QC, (ii) QC, 
    (iii) Data with AOC QC, (vi) Data with cloudiness index w'''
    # 
    
    
    # 1. plot of all data records
    # creates plotting mask ,Z,  where records exist
    Z = np.zeros([4, 78])
    for i in range(78):
        for j in range(78):
            if Ed_NASA.index[j] == i and ~np.isnan(np.array(Ed_NASA['400'])[j]) == 1:
                Z[0,i] = 1
            if Ed_HEREON.index[j] == i and ~np.isnan(np.array(Ed_HEREON['400'])[j]) == 1:
                Z[1,i] = 1
            if Ed_TARTU.index[j] == i and ~np.isnan(np.array(Ed_TARTU['400'])[j]) == 1:
                Z[2,i] = 1
            if Ed_PML.index[j] == i and ~np.isnan(np.array(Ed_PML['400'])[j]) == 1:
                Z[3,i] = 1
    
    # figure for stations
    plt.figure(figsize=(12,8))
    plt.rc('font', size=18)       
    cmap1 = matplotlib.cm.get_cmap('coolwarm')
    plt.pcolor(Z,edgecolors='k',cmap=cmap1)
    plt.title('AAOT station mask')
    
    plt.xlabel('Station number')
    ylabels = ['NASA','HEREON', 'TARTU','PML']
    plt.gca().set_yticklabels(ylabels)
    plt.gca().set_yticks([0.5,1.5,2.5,3.5])
    
    red_patch = mpatches.Patch(color=cmap1(1.001), label='Data')
    blue_patch = mpatches.Patch(color=cmap1(0), label='No data')
    plt.legend(handles=[red_patch,blue_patch],loc=3)
    filename  =  path_output +  'stationmask_4teams.png'
    plt.savefig(filename)
    
    # 2. plot of QC
    # creates plotting mask ,Z,  where records exist
    Z = np.zeros([2, 78])
    Z[0,:] = Q_mask['QC_CI_3']
    Z[1,:] = Q_mask['QC_AOC_3']
  
    # figure for stations
    plt.figure(figsize=(14,6))
    plt.rc('font', size=18)      
    plt.pcolor(Z,edgecolors='k',cmap='coolwarm')
    plt.title('Quality control mask')
    ylabels = ['2. Env.' + '\n' + '& CV[$R_{rs}$]', '1. AOC:'  +'\n' +  'SEAPRISM']
    plt.xlabel('Station number')
    plt.gca().set_yticklabels(ylabels)
    plt.gca().set_yticks([0.5,1.5])
    
    red_patch = mpatches.Patch(color=cmap1(1.001), label='Pass QC')
    blue_patch = mpatches.Patch(color=cmap1(0), label='Fail QC')
    plt.legend(handles=[red_patch,blue_patch],loc=4)  
    
    filename  =  path_output + 'QCmask.png'
    plt.savefig(filename,dpi=300)
    
    # 3. plot of data wwith AOC QC
    
    # creates plotting mask ,Z, where records exist
    Z = np.zeros([4, 78])
    for i in range(78):
        for j in range(78):
            if Ed_NASA.index[j] == i and ~np.isnan(np.array(Ed_NASA['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[0,i] = 1
            if Ed_HEREON.index[j] == i and ~np.isnan(np.array(Ed_HEREON['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[1,i] = 1
            if Ed_TARTU.index[j] == i and ~np.isnan(np.array(Ed_TARTU['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[2,i] = 1
            if Ed_PML.index[j] == i and ~np.isnan(np.array(Ed_PML['400'])[j]) == 1 and Q_mask['QC_AOC_3'][j] == 1:
                Z[3,i] = 1
    
    # figure for stations
    plt.figure(figsize=(12,8))
    plt.rc('font', size=18)         
    plt.pcolor(Z,edgecolors='k',cmap='coolwarm')
    plt.title('AAOT station mask:' + '\n' +  'AOC-SEAPRISM QC and $\geq$ 3 reference systems with data')
    plt.xlabel('Station number')
    ylabels = ['NASA','HEREON', 'TARTU','PML']
    plt.gca().set_yticklabels(ylabels)
    plt.gca().set_yticks([0.5,1.5,2.5,3.5])
    
    red_patch = mpatches.Patch(color=cmap1(1.001), label='Data')
    blue_patch = mpatches.Patch(color=cmap1(0), label='No data')
    plt.legend(handles=[red_patch,blue_patch],loc=3)
    
    filename  =  path_output  + 'stationmask_QC1_4teams.png'
    plt.savefig(filename)
    
    # 3. plot of data with Cloudiness Index QC

    # creates plotting mask ,Z,  where records exist
    Z = np.zeros([4, 78])
    for i in range(78):
        for j in range(78):
            if Ed_NASA.index[j] == i and ~np.isnan(np.array(Ed_NASA['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[0,i] = 1
            if Ed_HEREON.index[j] == i and ~np.isnan(np.array(Ed_HEREON['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[1,i] = 1
            if Ed_TARTU.index[j] == i and ~np.isnan(np.array(Ed_TARTU['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[2,i] = 1
            if Ed_PML.index[j] == i and ~np.isnan(np.array(Ed_PML['400'])[j]) == 1 and Q_mask['QC_CI_3'][j] == 1:
                Z[3,i] = 1
    
    # figure for stations
    plt.figure(figsize=(12,8))
    plt.rc('font', size=18)         
    plt.pcolor(Z,edgecolors='k',cmap='coolwarm')
    plt.title('AAOT station mask:' + '\n' + 'C.I and CV[$R_{rs}$] QC and $\geq$ 3 reference systems with data')
    
    plt.xlabel('Station number')
    ylabels = ['NASA','HEREON', 'TARTU','PML']
    plt.gca().set_yticklabels(ylabels)
    plt.gca().set_yticks([0.5,1.5,2.5,3.5])
    
    red_patch = mpatches.Patch(color=cmap1(1.001), label='Data')
    blue_patch = mpatches.Patch(color=cmap1(0), label='No data')
    plt.legend(handles=[red_patch,blue_patch],loc=3)
    
    filename  =  path_output  + 'stationmask_QC2_4teams.png'
    plt.savefig(filename)
    
    return



def filter_by_azimuth(df_PML, df_NASA, df_RBINS, df_CNR, tol=1):
    '''function to filter NASA, RBINS, CNR via azimuth'''   
    
    delta_phi_NASA = abs(df_NASA['azimuth']) - abs(df_PML['azimuth'])
    delta_phi_RBINS = abs(df_RBINS['azimuth']) - abs(df_PML['azimuth'])
    delta_phi_CNR = abs(df_CNR['azimuth']) - abs(df_PML['azimuth'])
   
    # delta_phi_NASA = df_NASA['azimuth'] - df_PML['azimuth']
    # delta_phi_RBINS = df_RBINS['azimuth'] - df_PML['azimuth']
    
    for i in range(78):
       if abs(delta_phi_NASA[i]) > tol:
         df_NASA.iloc[i] = np.nan
       if abs(delta_phi_RBINS[i]) > tol:
         df_RBINS.iloc[i] = np.nan
       if abs(delta_phi_CNR[i]) > tol:
         df_CNR.iloc[i] = np.nan
         
    return df_NASA, df_RBINS, df_CNR


def filter_by_azimuth_CP(df_PML, df_NASA, df_NASA_CP, tol=1):
    '''function to filter just NASA by azimuth - used for CP
    IP data frame is used for azimuth'''   
    
    delta_phi_NASA = abs(df_NASA['azimuth']) - abs(df_PML['azimuth'])

    # delta_phi_NASA = df_NASA['azimuth'] - df_PML['azimuth']
    
    for i in range(78):
       if abs(delta_phi_NASA[i]) > tol:
         df_NASA_CP.iloc[i] = np.nan
         
    return df_NASA_CP


def azimuth_plot(Ed_PML, Ed_NASA, Ed_RBINS, Ed_CNR, path_output):
    
    plt.figure(figsize=(20,16))
    plt.rc('font', size=18) 
  #  plt.suptitle('Relative azimuth for NASA, RBINS, CNR, compared with fixed systems')
  
    plt.title('NASA')
    plt.subplot(2,2,1)
    plt.plot(Ed_PML.index, abs(Ed_PML['azimuth']), label= 'PML, TARTU, HEREON',color='orange',linewidth = 3)
    plt.plot(Ed_NASA.index, abs(Ed_NASA['azimuth']), label= 'NASA', linestyle='dashed',color='blue',linewidth = 3)
    plt.xlabel('Station number')
    plt.ylabel('|$\phi$| (deg)')
    plt.legend()
    plt.title('RBINS')
    plt.subplot(2,2,2)
    plt.plot(Ed_PML.index, abs(Ed_PML['azimuth']), label= 'PML, TARTU, HEREON',color='orange',linewidth = 3)
    plt.plot(Ed_RBINS.index, abs(Ed_RBINS['azimuth']), label= 'RBINS', linestyle='dashed', color='black',linewidth =3)
    plt.xlabel('Station number')
    plt.ylabel('|$\phi$| (deg)')
    plt.legend()
    
    plt.title('CNR')
    plt.subplot(2,2,3)
    plt.plot(Ed_PML.index, abs(Ed_PML['azimuth']), label= 'PML, TARTU, HEREON',color='orange',linewidth = 3)
    plt.plot(Ed_CNR.index, abs(Ed_CNR['azimuth']), label= 'CNR', linestyle='dashed',color='green',linewidth = 3)
    plt.xlabel('Station number')
    plt.ylabel('|$\phi$| (deg)')
    plt.legend()
    
    plt.tight_layout(pad=2.4)
    
    #plt.figure()
    #plt.title('Relative azimuth')
    #plt.plot(Ed_PML.index, Ed_PML['azimuth'], label= 'PML, TARTU, HEREON',color='orange',  linewidth = 3)
    #plt.plot(Ed_NASA.index, Ed_NASA['azimuth'], label= 'NASA', linestyle='dashed',color='blue', linewidth = 3)
    #plt.plot(Ed_RBINS.index, Ed_RBINS['azimuth'], label= 'RBINS', linestyle='dashed', color='red', linewidth = 3)
    #plt.xlabel('Station number')
    #plt.ylabel('$\phi$ (deg)')
    #plt.legend()
        
    plt.tight_layout()
    filename  =  path_output +  '/azimuth.png'
    plt.savefig(filename)

    return


