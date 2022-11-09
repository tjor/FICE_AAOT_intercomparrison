#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 09:58:07 2022

Script that contains file reader functions for AAOT intercomparrision excercise
CP = `Community Processed' 

Th was previously in intercomparrison master

@author: tjor
"""


from csv import reader
import numpy as np
import pandas as pd
import datetime 
import os

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.dates as mdates