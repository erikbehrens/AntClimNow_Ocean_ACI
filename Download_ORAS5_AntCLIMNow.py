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
import zipfile
from netCDF4 import Dataset 
import numpy as np
import copy
import matplotlib.pyplot as plt
c = cdsapi.Client()

start_date = datetime(2015, 1, 1,0,0)
end_date = datetime(2023, 10, 1,0,0)

data_dir='/home/behrense/_niwa02764n/ORAS5_download/'
run_dir='/home/behrense/analyse/ORAS5/'
mesh_dir='/home/behrense/_niwa02764p/MASK025/'


#%% downloads monthly zip files and unpacks them into netcdf

def download_ORAS5(year,month):
    var_to_download=['zonal_velocity','meridional_velocity']
    for var in  var_to_download: 
       # start_year =start_year+relativedelta(years=+1) 
        
        if year < 2015:
            pro_typ='consolidated'
        else:
            pro_typ='operational'
            
        file=data_dir+'/ORAS5_1m_'+var+'_'+str(year)+month+'.zip'
        c.retrieve(
            'reanalysis-oras5',
        {
        'vertical_resolution':'all_levels',
        'variable': var,
        'product_type': pro_typ,
        'year': str(year),
        'month': [ month ],

                'format': 'zip',
            },
            file)
        with zipfile.ZipFile(file,"r") as zip_ref:
            zip_ref.extractall(data_dir)





#%% calculates barotropic transport for ORAS5


def calculate_barotropic_transport(year,month):
    if year < 2015:
   
        pro_typ='CONS'
    else:
        pro_typ='OPER'
    u_file='vozocrtx_control_monthly_highres_3D_'+str(year)+month+'_'+pro_typ+'_v0.1.nc'
    
    cf_ufil = Dataset(data_dir + u_file)
    print(cf_ufil)
    npiglo  = len(cf_ufil.dimensions["x"])
    npjglo  = len(cf_ufil.dimensions["y"])
    npk     = len(cf_ufil.dimensions["depthu"])

    zu      = cf_ufil.variables['vozocrtx'][0,:,:,:]
    zu[zu>100]=0
    print('hier1') 
    cf_ufil.close()

    cn_mask = Dataset(mesh_dir + "/mesh_mask.nc")
    glamf   = cn_mask.variables["glamf"][0, :, :]
    gphif   = cn_mask.variables["gphif"][0, :, :]
    e2u     = cn_mask.variables["e2u"][0, :, :]
    e3u     = cn_mask.variables["e3t_0"][0,  :]
    umask=cn_mask.variables['umask'][0,:,:,:]
    
    print('hier')
    cn_mask.close()
    # opt_dic["iiref"] = npiglo
    # opt_dic["ijref"] = npjglo
    
    dtrpu = np.zeros(e2u.shape)
    for jk in range(npk):
        dtrpu += zu[jk, :, :] * e2u[:, :] * e3u[jk]*umask[jk,:,:]

    # do meridional integration
    dpsiu = np.zeros(e2u.shape)
    for jj in range(1, npjglo):
        dpsiu[jj, :] = dpsiu[jj-1, :] - dtrpu[jj, :]

    dpsi = copy.deepcopy(dpsiu)
    dpsi=dpsi-dpsi[430,1230]
    return dpsi
    
test=calculate_barotropic_transport(2023,'10')

#%%
plt.close('all')
fig = plt.figure(figsize=(16,8))
plt.pcolormesh(test[:,:])
plt.colorbar()

#%%check last entry if new data needs processing



def get_last_entry():
    df = pd.read_csv(run_dir+'AntClimNow_Ocean_ACI.csv')
    last_date=df[-1:]['Date'].values[0]
    last_date=datetime.strptime(last_date,'%Y-%m-%d')
    return last_date

last_date=get_last_entry()
if last_date <= end_date:
    
    dates = [dt for dt in rrule(MONTHLY, dtstart=last_date, until=end_date)]  
    for da in dates[0]: #only donwload next month
        print(da.year,str(da.month).zfill(2))
        download_ORAS5(da.year,str(da.month).zfill(2))
    calculate_barotropic_transport((da.year,str(da.month).zfill(2)))
   #return
 

#%%

        
   
    
