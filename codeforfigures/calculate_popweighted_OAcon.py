#!/usr/bin/env python

# Imports
import gcpy.constants as gcon
import os
import numpy as np
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import xarray as xr
import warnings
import gcpy.plot as gcplot
import cartopy.crs as ccrs
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
                    if s in f and t in f:
                        file_list.append(os.path.join(root, f))

    # Return an alphabetically sorted list of files
    file_list = list(set(file_list))
    file_list.sort()
    print('file found {}'.format(file_list))
    return file_list

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
    # Path for China mask
    maskpath = '/home/mrq/pycode/mask/china_map/'
    chnmaskfile = maskpath + 'china.mask.generic.05x0625.nc'
    chnmask = xr.open_dataset(chnmaskfile)
    # Path for population data
    popfile = '/home/mrq/PopulationGrid/GPW_v4/gpwv4_population_2015c_05x0625.nc'
    popds = xr.open_dataset(popfile)

    # ----------------------------------------------------------------------
    # Read and calculate concentration data
    # ----------------------------------------------------------------------
    # Get a list of files in the SpeciesConc and StateMet collections in GEOS-Chem output
    # (YOU CAN EDIT THIS FOR YOUR OWN PARTICULAR APPLICATION!)
    collections = ['SpeciesConc','StateMet']

    # Set starttime and endtime to calculate time range
    year = 2017
    starttime = '1/1/'+str(year)
    endtime = '12/1/'+str(year)
    settime = calculate_timerange(starttime, endtime)
    season = 'annual'

    # Read GEOS-Chem data into an xarray Dataset
    ds = read_geoschem_data(path_to_data, collections, settime)

    # Find the species names
    species_name = ['POA1', 'POA2', 'POA3', 'POA4', 'POA5',
                    'OPOA1', 'OPOA2', 'OPOA3',
                    'ASOA1', 'ASOA2', 'ASOA3', 'ASOA0',
                    'TSOA0', 'TSOA1', 'TSOA2', 'TSOA3',
                    'ISOA0', 'ISOA1', 'ISOA2', 'ISOA3',
                    'IVSOA0', 'IVSOA1', 'IVSOA2', 'IVSOA3', 'IVSOAN',
                    'SOAGX', 'SOAIE', 'LVOCOA', 'IONITA']
    MW = [16.814, 16.814, 16.814, 16.814, 16.814,
          25.221, 25.221, 25.221,
          150, 150, 150, 150,
          150, 150, 150, 150,
          150, 150, 150, 150,
          150, 150, 150, 150, 150,
          58.04, 118.15, 154.19, 14.01]
    POA_species = ['POA1', 'POA2', 'POA3', 'POA4', 'POA5' ]
    OPOA_species = ['OPOA1', 'OPOA2', 'OPOA3']
    ASOA_species = ['ASOA1', 'ASOA2', 'ASOA3', 'ASOA0']
    IVSOA_species = ['IVSOA0', 'IVSOA1', 'IVSOA2', 'IVSOA3', 'IVSOAN']
    BSOA_species = ['TSOA0', 'TSOA1', 'TSOA2', 'TSOA3', 'ISOA0', 'ISOA1', 'ISOA2', 'ISOA3']
    aqSOA_species = ['SOAGX', 'SOAIE', 'LVOCOA', 'IONITA']

    ppb_ugm3 = 1e3 * 1e6 / 28.9647
    nspecies = len(species_name)
    print('nspecies:', nspecies)

    # Read air density data from the GEOS-Chem xarray Dataset
    dataAirDen = ds['Met_AIRDEN'].isel(lev=0)

    # Set xarray for OA concentrations
    time = ds['time'].values
    lat = ds.lat.values
    lon = ds.lon.values
    data = np.zeros((time.size, lat.size, lon.size), dtype=float)
    OAcon = xr.DataArray(data, coords=[time, lat, lon], dims=['time', 'lat', 'lon'])
    POAcon = xr.DataArray(data, coords=[time, lat, lon], dims=['time', 'lat', 'lon'])
    OPOAcon = xr.DataArray(data, coords=[time, lat, lon], dims=['time', 'lat', 'lon'])
    ASOAcon = xr.DataArray(data, coords=[time, lat, lon], dims=['time', 'lat', 'lon'])
    BSOAcon = xr.DataArray(data, coords=[time, lat, lon], dims=['time', 'lat', 'lon'])
    IVSOAcon = xr.DataArray(data, coords=[time, lat, lon], dims=['time', 'lat', 'lon'])
    aqSOAcon = xr.DataArray(data, coords=[time, lat, lon], dims=['time', 'lat', 'lon'])
    print(time.size, lat.size, lon.size)

    # Read concentration data from the GEOS-Chem xarray Dataset
    for species in range(nspecies):
        vname = 'SpeciesConc_' + species_name[species]
        species_con = ds[vname].isel(lev=0)
        OAcon = OAcon + species_con * dataAirDen * MW[species] * ppb_ugm3
        if species_name[species] in POA_species:
            POAcon = POAcon + species_con * dataAirDen * MW[species] * ppb_ugm3
        if species_name[species] in OPOA_species:
            OPOAcon = OPOAcon + species_con * dataAirDen * MW[species] * ppb_ugm3
        if species_name[species] in ASOA_species:
            ASOAcon = ASOAcon + species_con * dataAirDen * MW[species] * ppb_ugm3
        if species_name[species] in BSOA_species:
            BSOAcon = BSOAcon + species_con * dataAirDen * MW[species] * ppb_ugm3
        if species_name[species] in IVSOA_species:
            IVSOAcon = IVSOAcon + species_con * dataAirDen * MW[species] * ppb_ugm3
        if species_name[species] in aqSOA_species:
            aqSOAcon = aqSOAcon + species_con * dataAirDen * MW[species] * ppb_ugm3
    SOAcon = OPOAcon + ASOAcon + BSOAcon + IVSOAcon + aqSOAcon

    # Choose components and regions to output population-weighted concentrations
    plot = ['POA', 'SOA']
    region = ['chn']
    for rg in region:
        rgmask = eval(rg + 'mask')
        rgpop = popds['population'] * rgmask['mask']
        for pl in range(len(plot)):
            spcon = eval(plot[pl] + 'con')
            rgconpop = spcon.mean(dim='time') * rgpop
            print(rgconpop.sum().values / rgpop.sum().values)

if __name__ == "__main__":
    main()
