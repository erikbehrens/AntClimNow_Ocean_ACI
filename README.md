# AntClimNoW    
# AntClimNow Ocean ACI

This repository contains the python scripts to calculate following oceanic Antarctic Climate Indicators (ACI):
1. Barotropic Ross Gyre transport [Sv]
2. Barotropic Weddell Gyre transport [Sv]
3. Barotropic Drake Passage transport [Sv]
4. Antarctic Bottom Water volume [km^3]
5. Antartcic Bottom Water temperature [degrees C] (potential temperature)
6. Antarctic Bottom Water salinity [psu]
7. High Salinity Shelf Water Ross Sea volume [km^3]
8. High Salinity Shelf Water Ross Sea temperature [degrees C] (potential temperature)
9. High Salinity Shelf Water Ross Sea salinity [psu]
10. High Salinity Shelf Water Weddell Sea volume [km^3]
11. High Salinity Shelf Water Weddell Sea temperature [degrees C] (potential temperature)
12. High Salinity Shelf Water Weddell Sea salinity [psu]
13. Antarctic Bottom water export northwar across 35S [Sv]
14. Strength of the Bottom water overturning cell [Sv]

Monthly means of this diagnostics are provide in AntClimNow_Ocean_ACI_v1.csv for the period 1960 to 2023.
Data originates from ORAS5 reanalysis datasets. (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-oras5?tab=form)

# What is needed to run the code?
1. A CDS account (https://cds.climate.copernicus.eu/user/register)
2. A CDS python environment ( a working environment.yml file is proivde) conda env create -f environment.yml
3. Three paths wihtin the AntClimNow_Ocean_ACI.py are required to be modified in main()
  3.1     data_dir='/home/ORAS5_download/'  # this is the directory where the ORAS5 data is donwloaded to, on HPCF that might be a temporary directory
  3.2     run_dir='/home/analyse/ORAS5/'    # that is the run directory contain the cloned python scripts
  3.3     mesh_dir='/home//MASK025/'        # this folder points to grid specification of ORAS5 (mesh_mask.nc, file is here https://icdc.cen.uni-hamburg.de/thredds/catalog/ftpthredds/EASYInit/oras5/ORCA025/mesh/catalog.html and in zenodo)
4.  Modify
    start_date = datetime(1960, 1, 1,0,0)
    end_date = datetime(2023, 12, 1,0,0) in main() as required to rerun or update the existing ACI values.
    The code takes the last time stamp (AntClimNow_Ocean_ACI_v1.csv) into account and starts from that point in time and continues until end_date is reached. It simply append
    updated values in AntClimNow_Ocean_ACI_v1.csv. It also checks if data has been already downloaded, which saves a bit of time.
5. run python AntClimNow_Ocean_ACI.py (depending on the download speed it might take around 10-15 minutes to obtain a new monthly values for all ACI)

# Diagnostics
ACI diagnostics follow widely used definitions, but might not satisfy everyone. (more documentation is here https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JC011286)
1. Ross and Weddell Gyre strength is computed by using the barotropic transport between gyre averaged value and the nearest costal value. The averaged region covers the gyre centred and is large enough to not just pick up a single eddy.
2. Drake Passage transport is the transport difference between both closest land point from the barotropic streamfunction
3. AABW bottom water, sigma2 density has been calculated and the volume south of the equaltor with sigma denser of 37 kg/m^3 integrated. Median values for salinity and temperature have been provided too to characterise changes in the water mass. It appears ORAS5 does not produce very dense observed AABW values, therefore sigma2 values have slightly adjusted towards lighter water
4. HSSW water for Ross and Weddell Sea has been computed similar to AABW above, but with sigma2 denser than 37.05 and shallower than 1000m and just of over the formation region in Ross and Weddell Sea
5. Overturning and AABW export, the overturning has been calculted in sigma2 coordinates and water denser than 36.8 across 35S diagnosed. Overtunring strenth is quantified by the strength of streamfunction denser 36.8 between 90S and 35S.


# Funding 
This toolbox was supported by funding through SCAR/AntClimNow https://www.scar.org/science/antclimnow/home
    

