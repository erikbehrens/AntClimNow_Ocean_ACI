#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  7 23:51:46 2023

@author: behrense
"""




import cdsapi
from datetime import datetime
import pandas as pd
from dateutil.rrule import rrule, MONTHLY
import zipfile
from netCDF4 import Dataset 
import numpy as np
import glob
c = cdsapi.Client()
from dateutil.relativedelta import *



#%% downloads monthly zip files and unpacks them into netcdf

def download_ORAS5(year,month):
    var_to_download=['zonal_velocity','meridional_velocity','potential_temperature','salinity']
    varnames=['vozocrtx','vomecrty','votemper','vosaline']
   # var_to_download=['potential_temperature','salinity']
    
    
    
    for var,names in  zip(var_to_download,varnames): 
        if len(glob.glob(data_dir+'*'+names+'*_'+str(year)+month+'*.nc'))==0: # file needs downloading  
            print(var)
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
def nearest_point(lon_model,lat_model,lon_target,lat_target):
    ''' returns the nearest indices (i,j) from given lon lat in nav_lon nav_lat array '''
    import scipy.spatial
    import numpy as np
    points=np.vstack([lon_target,lat_target]).T
    comb_array=np.dstack([lon_model.ravel(),lat_model.ravel()])[0]
    mytree = scipy.spatial.cKDTree(comb_array)
    dist, indexes = mytree.query(points)
    j,i=np.unravel_index(indexes,lon_model.shape)
    return i[0],j[0]


# compute barotropic streamfunction and calculte transport for Ross and Weddell Gyre and Drake transports
def get_barotropic_transport(year,month):

    if year < 2015:
   
        pro_typ='CONS'
    else:
        pro_typ='OPER'
    #read zonal velocity data    
    u_file='vozocrtx_control_monthly_highres_3D_'+str(year)+month+'_'+pro_typ+'_v0.1.nc'
    data = Dataset(data_dir + u_file)
    lon   =  data.variables['nav_lon'][:,:]
    lon[lon<0]+=360
    lat   =  data.variables['nav_lat'][:,:]
    npjglo  = len(data.dimensions["y"])
    npk     = len(data.dimensions["depthu"])
    vozocrtx      = data.variables['vozocrtx'][0,:,:,:]
    vozocrtx[vozocrtx>100]=0
    #read mesh information    
    cn_mask = Dataset(mesh_dir + "/mesh_mask.nc")
    e2u     = cn_mask.variables["e2u"][0, :, :]
    e3u     = cn_mask.variables["e3t_0"][0,  :]
    umask=cn_mask.variables['umask'][0,:,:,:]
    
    #integrate zonal transports vertically  
    dtrpu = np.zeros(e2u.shape)
    for jk in range(npk):
        dtrpu += vozocrtx[jk, :, :] * e2u[:, :] * e3u[jk]*umask[jk,:,:]

    # do meridional integration
    streamfunction = np.zeros(e2u.shape)
    for jj in range(1, npjglo):
        streamfunction[jj, :] = streamfunction[jj-1, :] - dtrpu[jj, :]

    # set zeros value for Africa
    streamfunction=streamfunction-streamfunction[430,1230]

    
    #definition for Ross Gyre
    xx1,yy1=nearest_point(lon, lat, 200,-70)
    xx2,yy2=nearest_point(lon, lat, 220,-65)
    Ross_gyre=(np.mean(np.mean(streamfunction[yy1:yy2,xx1:xx2],1),0)-streamfunction[27,xx2])/1e6
    
    #definition for Ross Gyre
    xx1,yy1=nearest_point(lon, lat, 291.7,-55.53)
    xx2,yy2=nearest_point(lon, lat, 299.3,-64.05)
    ACC_transport=(streamfunction[yy2,xx2]-streamfunction[yy1,xx1])/1e6
    
    #definition for Weddell Gyre
    xx1,yy1=nearest_point(lon, lat, 350,-62.5)
    xx2,yy2=nearest_point(lon, lat, 10,-57.5)
  
    Weddell_gyre=(np.mean(np.mean(streamfunction[yy1:yy2,xx1:xx2],1),0)-streamfunction[101,xx2])/1e6
    return Ross_gyre,Weddell_gyre,ACC_transport
    
# Ross_gyre,Weddell_gyre,ACC_transport=get_barotropic_transport(2023,'10')

#%%
def get_AABW_HSSW_properties(year,month):
     import eos
   
     if year < 2015:
        pro_typ='CONS'
     else:
        pro_typ='OPER'
        
     #read temperature and salinity data
     t_file='votemper_control_monthly_highres_3D_'+str(year)+month+'_'+pro_typ+'_v0.1.nc'
     data = Dataset(data_dir + t_file)
     lon   =  data.variables['nav_lon'][:,:]
     lon[lon<0]+=360
     lat   =  data.variables['nav_lat'][:,:]
     npk     = len(data.dimensions["deptht"])
     votemper= data.variables['votemper'][0,:,:,:]
     s_file='vosaline_control_monthly_highres_3D_'+str(year)+month+'_'+pro_typ+'_v0.1.nc'
     data = Dataset(data_dir +s_file)
     vosaline= data.variables['vosaline'][0,:,:,:]
     # compute sigma2
     sigma2 = eos.sigmai_dep(votemper[:, :, :], vosaline[:, :, :], 2000)
     # get the mesh information
     cn_mask = Dataset(mesh_dir + "/mesh_mask.nc")
     e1t     = cn_mask.variables["e1t"][0, :, :]
     tmask   = cn_mask.variables["tmask"][0,:, :, :]
     e2t     = cn_mask.variables["e2t"][0, :, :]
     e3t     = cn_mask.variables["e3t_0"][0, :]
     # volume fore each grid cell
     volume=np.zeros(votemper.shape)
     for k in range(npk):
         
         volume[k,:,:]=e1t*e2t*e3t[k]*tmask[k,:,:]
     # integrate AABW volume where sigma is >37 up to the equator
     AABW_volume=np.sum(np.where(sigma2[:,:550,:]>37,volume[:,:550,:],0))/1e9 #km3
     AABW_temperature=np.nanmedian(np.where(sigma2[:,:550,:]>37,votemper[:,:550,:],np.nan))
     AABW_salinity=np.nanmedian(np.where(sigma2[:,:550,:]>37,vosaline[:,:550,:],np.nan))
     # integrate volume of HSSW in Ross Sea sector where sigma >37.05 and get median temperature and salinity
     sigma_threshold=37.05
     depth_level=45 #upper 1000m
     xx1,yy1=nearest_point(lon, lat, 160,-70) # only need to the northern end
     xx2,yy2=nearest_point(lon, lat, 200,-70)
     HSSW_Ross_volume=np.sum(np.where(sigma2[:depth_level,:yy1,xx1:xx2]>sigma_threshold,volume[:depth_level,:yy1,xx1:xx2],0))/1e9 #upper 1000m in km
     HSSW_Ross_temperature=np.nanmedian(np.where(sigma2[:depth_level,:yy1,xx1:xx2]>sigma_threshold,votemper[:depth_level,:yy1,xx1:xx2],np.nan))
     HSSW_Ross_salinity=np.nanmedian(np.where(sigma2[:depth_level,:yy1,xx1:xx2]>sigma_threshold,vosaline[:depth_level,:yy1,xx1:xx2],np.nan))
     # integrate volume of HSSW in Weddell Sea sector where sigma >37.05 and get median temperature and salinity
     xx1,yy1=nearest_point(lon, lat, 295,-70) # only need to the northern end
     xx2,yy2=nearest_point(lon, lat, 335,-70)
     HSSW_Weddell_volume=np.sum(np.where(sigma2[:depth_level,:yy1,xx1:xx2]>sigma_threshold,volume[:depth_level,:yy1,xx1:xx2],0))/1e9 #upper 1000m in km
     HSSW_Weddell_temperature=np.nanmedian(np.where(sigma2[:depth_level,:yy1,xx1:xx2]>sigma_threshold,votemper[:depth_level,:yy1,xx1:xx2],np.nan))
     HSSW_Weddell_salinity=np.nanmedian(np.where(sigma2[:depth_level,:yy1,xx1:xx2]>sigma_threshold,vosaline[:depth_level,:yy1,xx1:xx2],np.nan))
     return AABW_volume, AABW_temperature,AABW_salinity,HSSW_Ross_volume,HSSW_Ross_temperature,HSSW_Ross_salinity,HSSW_Weddell_volume,HSSW_Weddell_temperature,HSSW_Weddell_salinity
       
#%%


def get_sigma_overturning(year,month):    
    import eos
    if year < 2015:
        pro_typ='CONS'
    else:
        pro_typ='OPER'
    #read temperature and salinity data to compute sigma2    
    t_file='votemper_control_monthly_highres_3D_'+str(year)+month+'_'+pro_typ+'_v0.1.nc'
    s_file='vosaline_control_monthly_highres_3D_'+str(year)+month+'_'+pro_typ+'_v0.1.nc'
    data = Dataset(data_dir + t_file)
    lon   =  data.variables['nav_lon'][:,:]
    lon[lon<0]+=360
    lat   =  data.variables['nav_lat'][:,:]
    npiglo  = len(data.dimensions["x"])
    npjglo  = len(data.dimensions["y"])
    votemper= data.variables['votemper'][0,:,:,:]
    data = Dataset(data_dir +s_file)
    vosaline= data.variables['vosaline'][0,:,:,:]
    
    
    #read meridional velocity data to compute sigma2    
    v_file='vomecrty_control_monthly_highres_3D_'+str(year)+month+'_'+pro_typ+'_v0.1.nc'
    data = Dataset(data_dir +v_file)
    vomecrty= data.variables['vomecrty'][0,:,:,:]
    #read mesh data
    cn_mask = Dataset(mesh_dir + "mesh_mask.nc")
    e1v     = cn_mask.variables["e1v"][0, :, :]
    e3v     = cn_mask.variables["e3t_0"][0, :]
    nbins, sigmin, sigstp =158,30,0.05
    
    # make the density vector
    sigma = sigmin + (np.linspace(1, nbins, nbins) - 0.5) * sigstp
    
    # calculate sigma2
    sigma2 = eos.sigmai_dep(votemper[:, :, :], vosaline[:, :, :], 2000)
    #  inflate matrix
    lat=np.repeat(lat[np.newaxis,:,:,],75,axis=0)
    e1v=np.repeat(e1v[np.newaxis,:,:,],75,axis=0)
    e3v=np.repeat(e3v[:,np.newaxis],npjglo,axis=1)
    e3v=np.repeat(e3v[:,:,np.newaxis],npiglo,axis=2)
    #  use sigma as mask so vector have the same length
    lat = np.ma.masked_where(np.ma.getmask(sigma2), lat)
    vomecrty[vomecrty>10]=0
    vomecrty= np.ma.masked_where(np.ma.getmask(sigma2), vomecrty.data)
    e1v= np.ma.masked_where(np.ma.getmask(sigma2), e1v)
    e3v= np.ma.masked_where(np.ma.getmask(sigma2), e3v)
    
    #latitude vector
    lat_ax=np.arange(-90,90,.25)
    # do 2d historgam with volume transports acuumulated
    hist2d_tmp=np.histogram2d(sigma2.compressed(),lat.compressed(),bins=(sigma,lat_ax),weights=vomecrty.compressed()*e3v.compressed()*e1v.compressed())
    #set zeros  nan
    hist2d=np.where(hist2d_tmp[0]==0,np.nan,hist2d_tmp[0])
    #generate streamfunction
    streamfunction=np.nancumsum(np.flipud(hist2d),axis=0)/1e6 #integrate from dense to light at each latitude -> Sv
    
    AABW_export= streamfunction[136:,210:230].max() #get maximum at 37.5:32.5S and denser than 36.8 which is the AABW export Sv
    AABW_overturning= streamfunction[136:,:220].max() #get maximum at 90S:35S and denser than 36.8 which is  AABW in ORAS5 Sv
    return AABW_export,AABW_overturning




#%%check last entry if new data needs processing



def get_last_entry():
    df = pd.read_csv(run_dir+'AntClimNow_Ocean_ACI_v1.csv')
    last_date=df[-1:]['Date'].values[0]
    last_date=datetime.strptime(last_date,'%Y-%m-%d')
    return last_date

def main():
    global start_date,end_date,data_dir,run_dir,mesh_dir,sigma2
    
    
    
    ### this section need adjusting by user
    start_date = datetime(1960, 1, 1,0,0)
    end_date = datetime(2023, 12, 1,0,0)
    
    data_dir='/home/behrense/_niwa02764n/ORAS5_download/'
    run_dir='/home/behrense/analyse/ORAS5/'
    mesh_dir='/home/behrense/_niwa02764p/MASK025/'
    ### this section above need adjusting by user
    
    #check for csv files
    if len(glob.glob(run_dir+'AntClimNow_Ocean_ACI_v1.csv'))==0: # file needs creating
        idx = pd.date_range("1957-12-01", periods=1, freq="MS")
        init_data=np.zeros((1,14,))
        df = pd.DataFrame(init_data,columns=['Ross_gyre','Weddell_gyre','ACC_transport','AABW_export','AABW_overturning',\
                                   'AABW_volume', 'AABW_temperature','AABW_salinity','HSSW_Ross_volume',\
                                   'HSSW_Ross_temperature','HSSW_Ross_salinity','HSSW_Weddell_volume',\
                                   'HSSW_Weddell_temperature','HSSW_Weddell_salinity'],index= idx)
        df.index.name='Date'
        df.to_csv(run_dir+'AntClimNow_Ocean_ACI_v1.csv')

    last_date=get_last_entry()
    if last_date <= end_date and start_date <=last_date:
        
        dates = [dt for dt in rrule(MONTHLY, dtstart=last_date+relativedelta(months=+1), until=end_date)]  
    elif last_date <= end_date and start_date >=last_date:
        
        dates = [dt for dt in rrule(MONTHLY, dtstart=start_date, until=end_date)]     
    if len(dates) !=0:
        for da in dates: #only donwload next month
            year=da.year
            month=str(da.month).zfill(2)
            print('Processing', year,month,'now. Takes about 3-5 minutes after download.')
            download_ORAS5(da.year,str(da.month).zfill(2))
            idx = pd.date_range(da, periods=1, freq="MS")
            #compute barotropic stremafunction
            Ross_gyre,Weddell_gyre,ACC_transport=get_barotropic_transport(year,month)
            #compute sigma2 overturning 
            AABW_export,AABW_overturning=get_sigma_overturning(year,month)
            #compute sigma2 
            AABW_volume, AABW_temperature,AABW_salinity,HSSW_Ross_volume,HSSW_Ross_temperature,\
            HSSW_Ross_salinity,HSSW_Weddell_volume,HSSW_Weddell_temperature,HSSW_Weddell_salinity=get_AABW_HSSW_properties(year,month)
            data=np.array([[Ross_gyre,Weddell_gyre,ACC_transport,AABW_export,AABW_overturning,\
                            AABW_volume, AABW_temperature,AABW_salinity,HSSW_Ross_volume,\
                            HSSW_Ross_temperature,HSSW_Ross_salinity,\
                            HSSW_Weddell_volume,HSSW_Weddell_temperature,HSSW_Weddell_salinity]])
            print(data)
            df = pd.DataFrame(data,index= idx)
            #df = pd.DataFrame(init_data+2,index= idx)
            df.to_csv(run_dir+'AntClimNow_Ocean_ACI_v1.csv',mode='a', header=False)
    return data       
if __name__ == '__main__':
    data=main()        

 
#%%

# da=Dataset('/home/behrense/_niwa02764p/MASK025//mesh_mask.nc')

# lon=da['nav_lon'][:]
# lon[lon<0]+=360
# lat=da['nav_lat'][:]

# import matplotlib.pyplot as plt

# plt.close('all')
# fig = plt.figure(figsize=(16,8))

# plt.pcolormesh(lon[:550,:1050],lat[:550,:1050],np.max(sigma2[:45,:550,:1050],axis=0),vmin=37.05)
# plt.colorbar()


# # #%%
# AABW_export,AABW_overturning=get_sigma_overturning(2023,'09')
# AABW_volume, AABW_temperature,AABW_salinity,sigma2,HSSW_Ross_volume,HSSW_Ross_temperature,HSSW_Ross_salinity,HSSW_Weddell_volume,HSSW_Weddell_temperature,HSSW_Weddell_salinity=get_AABW_HSSW_properties(2023,'09')
# Ross_gyre,Weddell_gyre,ACC_transport=get_barotropic_transport(2023,'09')
# #%%      
   
    
