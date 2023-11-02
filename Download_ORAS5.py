#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  7 23:51:46 2023

@author: behrense
"""




import cdsapi
from datetime import datetime
 
from dateutil.relativedelta import relativedelta

c = cdsapi.Client()

start_year = datetime(2015, 1, 1,0,0)
end_year = datetime(2024, 1, 1,0,0)
#delta = timedelta(month=1)




#%%
var_to_download=['zonal_velocity']
#%%

while start_year <= end_year:
   
    for var in  var_to_download: 
       # start_year =start_year+relativedelta(years=+1) 
        
        
        file='/home/behrense/_niwa02764n/ORAS5_download/ORAS5_1m_'+var+'_'+str(start_year.year)+'.zip'
        c.retrieve(
            'reanalysis-oras5',
        {
        'vertical_resolution':'all_levels',
        'variable': var,
        'product_type': 'operational',
        'year': str(start_year.year),
        'month': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                    ],

                'format': 'zip',
            },
            file)
        
    start_year =start_year+relativedelta(years=+1)
    
