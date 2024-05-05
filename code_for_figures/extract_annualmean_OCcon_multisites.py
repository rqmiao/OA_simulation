#!/usr/bin/env python

# Imports
import gcpy.constants as gcon
import os
import numpy as np
import xarray as xr
import warnings
import pandas as pd
from datetime import datetime

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

def find_files_in_dir(path, substrs, timerange):

    # Initialize
    file_list = []

    # Walk through the given data directory.  Then for each file found,
    # add it to file_list if it matches text in search_list.
    for root, directory, files in os.walk(path):
        for f in files:
            for t in timerange:
                for s in substrs:
                    if s in f and t in f: #  and '2020' in f and not '20201001' in f and not '20201101' in f: #
                        file_list.append(os.path.join(root, f))

    # Return an alphabetically sorted list of files
    file_list.sort()
    print('file found {}'.format(file_list))
    return file_list

def find_value_index(seq, val):
    '''
    Finds the index of a numpy array that is close to a value.
    '''
    r = np.where(np.diff(np.sign(seq - val)) != 0)
    idx = r + (val - seq[r]) / (seq[r + np.ones_like(r)] - seq[r])
    idx = np.append(idx, np.where(seq == val))
    idx = np.sort(idx)
    result = np.round(idx)

    # NOTE: xarray needs integer values, so convert here!
    return int(result[0])

def read_geoschem_data(path, collections, timerange):

    # Get a list of variables that GCPy should not read.
    # These are mostly variables introduced into GCHP with the MAPL v1.0.0
    # update.  These variables contain either repeated or non-standard
    # dimensions that can cause problems in xarray when combining datasets.
    skip_vars = gcon.skip_these_vars

    # Find all files in the given 
    file_list = find_files_in_dir(path, collections, timerange)

    # Return a single xarray Dataset containing data from all files
    # NOTE: Need to add combine="nested" for xarray 0.15 and higher
    v = xr.__version__.split(".")
    #print(v)
    return xr.open_mfdataset(file_list,
                             drop_variables=skip_vars)
    #if int(v[0]) == 0 and int(v[1]) >= 15:
    #    return xr.open_mfdataset(file_list,
    #                             drop_variables=skip_vars,
    #                             combine="nested",
    #                             concat_dim=None)
    #else:
    #    return xr.open_mfdataset(file_list,
    #                             drop_variables=skip_vars)
    
def calculate_timerange(starttime, endtime):
    time_range = pd.date_range(start=starttime, end=endtime, freq='MS')
    time_range = time_range.strftime('%Y%m%d')
    settime = time_range._data
    return settime

def main():
    '''
    Main program.
    '''
    # ----------------------------------------------------------------------
    # Set path where the data files live
    # ----------------------------------------------------------------------
    # Path for OA concentrations
    path_to_data = '/home/mrq/GEOS/13.3.1/updateHONO_sivocSOA_correctwl/rundirs/gc_05x0625_merra2_fullchem_SVPOA_mPOA_ch_meic_fuelbased_wdluo/OutputDir'
    # Path for site information
    input_path = '/home/mrq/pycode/'
    sitedata = pd.read_csv(input_path+'OC_network_sites.csv', usecols=['lon', 'lat'])
    nsites = len(sitedata.index)
    print('site number:', nsites)
    condata = np.zeros((nsites), dtype=float)

    # Set starttime and endtime to calculate time range
    year = 2018
    starttime = '1/1/' + str(year)
    endtime = '12/1/' + str(year)
    settime = calculate_timerange(starttime, endtime)

    # Find the species names
    species_name = ['POA1', 'POA2', 'POA3', 'POA4', 'POA5',
                    'OPOA1', 'OPOA2', 'OPOA3',
                    'ASOA1', 'ASOA2', 'ASOA3', 'ASOA0',
                    'TSOA0', 'TSOA1', 'TSOA2', 'TSOA3',
                    'ISOA0', 'ISOA1', 'ISOA2', 'ISOA3',
                    'IVSOA0', 'IVSOA1', 'IVSOA2', 'IVSOA3', 'IVSOAN',
                    'SOAGX', 'SOAIE', 'LVOCOA', 'IONITA']
    MW = [12.01, 12.01, 12.01, 12.01, 12.01,
          12.01, 12.01, 12.01,
          150/2.1, 150/2.1, 150/2.1, 150/2.1,
          150/2.1, 150/2.1, 150/2.1, 150/2.1,
          150/2.1, 150/2.1, 150/2.1, 150/2.1,
          150/2.1, 150/2.1, 150/2.1, 150/2.1, 150/2.1,
          12.01*2, 12.01*5, 12.01*5, 14.01/2.1]
    ppb_ugm3 = 1e3 * 1e6 / 28.9647
    nspecies = len(species_name)
    print('nspecies:', nspecies)
    
    # ----------------------------------------------------------------------
    # Read and calculate concentration data
    # ----------------------------------------------------------------------
    # Get a list of files in the SpeciesConc and StateMet collections in GEOS-Chem output
    # (YOU CAN EDIT THIS FOR YOUR OWN PARTICULAR APPLICATION!)
    collections = ['SpeciesConc', 'StateMet']

    # Read GEOS-Chem data into an xarray Dataset
    ds = read_geoschem_data(path_to_data, collections, settime)
    dataAirDen = ds['Met_AIRDEN'].isel(lev=0)
    time = ds['time'].values
    lat = ds.lat.values
    lon = ds.lon.values

    # Set xarray for OC concentrations
    data = np.zeros((time.size, lat.size, lon.size), dtype=float)
    con = xr.DataArray(data, coords=[time, lat, lon], dims=['time', 'lat', 'lon'])
    print(time.size, lat.size, lon.size)

    # Read concentration data from the GEOS-Chem xarray Dataset
    for species in range(nspecies):
        vname = 'SpeciesConc_' + species_name[species]
        species_con = ds[vname].isel(lev=0)
        scalar = dataAirDen * MW[species] * ppb_ugm3
        con = con + species_con * scalar

    # Calculate annual mean concentrations
    meancon = con.mean(dim='time')

    # Read concentration for each site
    for indexs in sitedata.index:
        site_coords = sitedata.loc[indexs].values

        # Find the indices corresponding to the site lon and lat
        lon_idx = find_value_index(ds.lon.values, site_coords[0])
        lat_idx = find_value_index(ds.lat.values, site_coords[1])
        condata[indexs] = meancon.isel(lon=lon_idx, lat=lat_idx)

    # Save csv file to disk
    columns_name = ['OC']
    df = pd.DataFrame(condata, columns=columns_name)
    df.to_csv('sites_OCcon.csv')

if __name__ == "__main__":
    main()
