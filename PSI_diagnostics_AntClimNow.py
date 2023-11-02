#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 22:52:01 2018

@author: behrense
"""
#Analyse of MHW

from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
from mpl_toolkits.basemap import Basemap
import sys
import numpy as np

sys.path.append('/home/behrense/python_scripts/ebtools/')
import plotting
import calc
import interpolation
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib as mpl
import cmocean

import pandas as pd
#%% read data that is an NZESM 95:14orical simulation
d1=Dataset('/home/behrense/_niwa02764p/MASK025/mesh_mask.nc')
lon=d1.variables['nav_lon'][:]
lon[lon<0]+=360
lat=d1.variables['nav_lat'][:]
tmask=d1['tmask'][:]
#%%
d1=Dataset('/home/behrense/_niwa02764n/ORAS5/psi/vozocrtx_control_monthly_highres_3D_195801_202309_OPER_v0.1.nc')

psi=d1['sobarstf'][:]

#%%
psi=psi*tmask[:,0,:,:]
psi=np.ma.masked_equal(psi, 0)
#%%
plt.close('all')
fig = plt.figure(figsize=(14,6))
plt.contour(lon[:300,1000:1150]-360,lat[:300,1000:1150],np.mean(psi[500:,:300,1000:1150],axis=0)/1e6,levels=55)
plt.contour(lon[:300,1155:],lat[:300,1155:],np.mean(psi[:120,:300,1155:],axis=0)/1e6,levels=55)
plt.colorbar()
plt.plot([200,220,220,200,200],[-70,-70,-65,-65,-70])

plt.plot([-10,10,10,-10,-10],[-62.5,-62.5,-57.5,-57.5,-62.5])

plt.plot([299,291],[-64,-55.5])
#%%

xx1,yy1=interpolation.nearest_point(lon, lat, 200,-70)
xx2,yy2=interpolation.nearest_point(lon, lat, 220,-65)
xx1=xx1[0]
xx2=xx2[0]
yy1=yy1[0]
yy2=yy2[0]

Ross_gyre=(np.mean(np.mean(psi[:,yy1:yy2,xx1:xx2],2),1)-psi[500,27,xx2])/1e6


xx1,yy1=interpolation.nearest_point(lon, lat, 291.7,-55.53)
xx2,yy2=interpolation.nearest_point(lon, lat, 299.3,-64.05)

xx1=xx1[0]
xx2=xx2[0]
yy1=yy1[0]
yy2=yy2[0]
ACC_transport=(psi[:,yy2,xx2]-psi[:,yy1,xx1])/1e6

xx1,yy1=interpolation.nearest_point(lon, lat, 350,-62.5)
xx2,yy2=interpolation.nearest_point(lon, lat, 10,-57.5)
xx1=xx1[0]
xx2=xx2[0]
yy1=yy1[0]
yy2=yy2[0]

Weddell_gyre=(np.mean(np.mean(psi[:,yy1:yy2,xx1:xx2],2),1)-psi[500,101,xx2])/1e6

#%%
plt.close('all')
fig = plt.figure(figsize=(14,6))
plt.plot(Weddell_gyre)

#%%
idx = pd.date_range("1958-01-01", periods=Weddell_gyre.shape[0], freq="MS")
data=np.array([Ross_gyre,Weddell_gyre,ACC_transport])
df = pd.DataFrame(data.T,columns=['Ross_gyre','Weddell_gyre','ACC_transport'],index= idx)


df.to_csv("AntClimNow_Ocean_ACI.csv")