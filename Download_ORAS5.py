#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  7 23:51:46 2023

@author: behrense
"""




import cdsapi
from datetime import datetime
from dateutil.relativedelta import relativedelta
import pandas as pd
from dateutil.rrule import rrule, MONTHLY

c = cdsapi.Client()

start_date = datetime(2015, 1, 1,0,0)
end_date = datetime(2024, 1, 1,0,0)

data_dir='/home/behrense/_niwa02764n/ORAS5_download/'
run_dir='/home/behrense/analyse/ORAS5/'

#%%check last entry
df = pd.read_csv(run_dir+'AntClimNow_Ocean_ACI.csv')
last_date=df[-1:]['Date'].values[0]
last_date=datetime.strptime(last_date,'%Y-%m-%d')

if last_date <= end_date:
    dates = [dt for dt in rrule(MONTHLY, dtstart=last_date, until=end_date)]  
    for da in dates:
        print(da.year,str(da.month).zfill(2))
        #download_ORAS5(da.year)
    
   #return
 

#%%

#%%

def download_ORAS5(year,month):
    var_to_download=['zonal_velocity','meridional_velocity']
    for var in  var_to_download: 
       # start_year =start_year+relativedelta(years=+1) 
        
        if year < 2015:
            pro_typ='consolidated'
        else:
            pro_typ='consolidated'
            
        file=data_dir+'/ORAS5_1m_'+var+'_'+str(year)+month+'.zip'
        c.retrieve(
            'reanalysis-oras5',
        {
        'vertical_resolution':'all_levels',
        'variable': var,
        'product_type': pro_typ,
        'year': year,
        'month': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                    ],

                'format': 'zip',
            },
            file)
        
   
    
