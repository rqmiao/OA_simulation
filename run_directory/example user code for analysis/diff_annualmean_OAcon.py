#!/usr/bin/env python

# example for drawing Figure 2c by Chen et al. (2024).

# Imports
import gcpy.constants as gcon
import os
import numpy as np
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
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
    
def calculate_timerange(starttime, endtime):
    time_range = pd.date_range(start=starttime, end=endtime, freq='MS')
    time_range = time_range.strftime('%Y%m%d')
    settime = time_range._data
    return settime

def main():
    '''
    Main program.
    '''
    # Path where the data files live
    # (YOU MUST EDIT THIS FOR YUR OWN PARTICULAR APPLICATION!)
    path_to_data1 = '/home/mrq/GEOS/13.3.1/updateHONO_sivocSOA_correctwl/rundirs/gc_05x0625_merra2_fullchem_SVPOA_mPOA_ch_meic_fuelbased_wdluo/OutputDir'
    path_to_data2 = path_to_data1
    maskpath = '/home/mrq/pycode/mask/china_map/'
    chnmaskfile = maskpath + 'china.mask.generic.05x0625.nc'
    chnmask = xr.open_dataset(chnmaskfile)
    xmask = xr.where(chnmask['mask'] == 0.0, np.NaN, 1.0)

    # Get a list of files in the ConcAboveSfc and SpeciesConc collections
    # (YOU CAN EDIT THIS FOR YOUR OWN PARTICULAR APPLICATION!)
    collections = ['SpeciesConc','StateMet']

    # Set starttime and endtime to calculate time range
    starttime1 = '1/1/2013'
    endtime1 = '12/1/2013'
    starttime2 = '1/1/2020'
    endtime2 = '12/1/2020'
    settime1 = calculate_timerange(starttime1, endtime1)
    settime2 = calculate_timerange(starttime2, endtime2)
    season = 'annual'

    # Read GEOS-Chem data into an xarray Dataset
    ds_ref = read_geoschem_data(path_to_data1, collections, settime1)
    ds_new = read_geoschem_data(path_to_data2, collections, settime2)
    dataAirDen_ref = ds_ref['Met_AIRDEN'].isel(lev=0)
    dataAirDen_new = ds_new['Met_AIRDEN'].isel(lev=0)

    # Define OA species in GC 13.3.1 and their MW
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
    BSOA_species = ['TSOA0', 'TSOA1', 'TSOA2', 'TSOA3','ISOA0', 'ISOA1', 'ISOA2', 'ISOA3']
    aqSOA_species = ['SOAGX', 'SOAIE', 'LVOCOA', 'IONITA']

    ppb_ugm3 = 1e3 * 1e6 / 28.9647
    nspecies = len(species_name)
    print('nspecies:', nspecies)

    time_ref = ds_ref['time'].values
    lat_ref = ds_ref.lat.values
    lon_ref = ds_ref.lon.values
    time_new = ds_new['time'].values
    lat_new = ds_new.lat.values
    lon_new = ds_new.lon.values
    data_ref = np.zeros((time_ref.size, lat_ref.size, lon_ref.size), dtype=float)
    data_new = np.zeros((time_new.size, lat_new.size, lon_new.size), dtype=float)
    OAcon_ref = xr.DataArray(data_ref, coords=[time_ref, lat_ref, lon_ref], dims=['time', 'lat', 'lon'])
    OAcon_new = xr.DataArray(data_new, coords=[time_new, lat_new, lon_new], dims=['time', 'lat', 'lon'])
    POAcon_ref = xr.DataArray(data_ref, coords=[time_ref, lat_ref, lon_ref], dims=['time', 'lat', 'lon'])
    POAcon_new = xr.DataArray(data_new, coords=[time_new, lat_new, lon_new], dims=['time', 'lat', 'lon'])
    OPOAcon_ref = xr.DataArray(data_ref, coords=[time_ref, lat_ref, lon_ref], dims=['time', 'lat', 'lon'])
    OPOAcon_new = xr.DataArray(data_new, coords=[time_new, lat_new, lon_new], dims=['time', 'lat', 'lon'])
    ASOAcon_ref = xr.DataArray(data_ref, coords=[time_ref, lat_ref, lon_ref], dims=['time', 'lat', 'lon'])
    ASOAcon_new = xr.DataArray(data_new, coords=[time_new, lat_new, lon_new], dims=['time', 'lat', 'lon'])
    BSOAcon_ref = xr.DataArray(data_ref, coords=[time_ref, lat_ref, lon_ref], dims=['time', 'lat', 'lon'])
    BSOAcon_new = xr.DataArray(data_new, coords=[time_new, lat_new, lon_new], dims=['time', 'lat', 'lon'])
    IVSOAcon_ref = xr.DataArray(data_ref, coords=[time_ref, lat_ref, lon_ref], dims=['time', 'lat', 'lon'])
    IVSOAcon_new = xr.DataArray(data_new, coords=[time_new, lat_new, lon_new], dims=['time', 'lat', 'lon'])
    aqSOAcon_ref = xr.DataArray(data_ref, coords=[time_ref, lat_ref, lon_ref], dims=['time', 'lat', 'lon'])
    aqSOAcon_new = xr.DataArray(data_new, coords=[time_new, lat_new, lon_new], dims=['time', 'lat', 'lon'])

    # calculate POA and SOA concentrations
    for species in range(nspecies):
        vname = 'SpeciesConc_' + species_name[species]
        species_con_ref = ds_ref[vname].isel(lev=0) * dataAirDen_ref * MW[species] * ppb_ugm3
        species_con_new = ds_new[vname].isel(lev=0) * dataAirDen_new * MW[species] * ppb_ugm3

        OAcon_ref = OAcon_ref + species_con_ref
        OAcon_new = OAcon_new + species_con_new
        if species_name[species] in POA_species:
            POAcon_ref = POAcon_ref + species_con_ref
            POAcon_new = POAcon_new + species_con_new
        if species_name[species] in OPOA_species:
            OPOAcon_ref = OPOAcon_ref + species_con_ref
            OPOAcon_new = OPOAcon_new + species_con_new
        if species_name[species] in ASOA_species:
            ASOAcon_ref = ASOAcon_ref + species_con_ref
            ASOAcon_new = ASOAcon_new + species_con_new
        if species_name[species] in BSOA_species:
            BSOAcon_ref = BSOAcon_ref + species_con_ref
            BSOAcon_new = BSOAcon_new + species_con_new
        if species_name[species] in IVSOA_species:
            IVSOAcon_ref = IVSOAcon_ref + species_con_ref
            IVSOAcon_new = IVSOAcon_new + species_con_new
        if species_name[species] in aqSOA_species:
            aqSOAcon_ref = aqSOAcon_ref + species_con_ref
            aqSOAcon_new = aqSOAcon_new + species_con_new
    SOAcon_ref = OPOAcon_ref + ASOAcon_ref + BSOAcon_ref + IVSOAcon_ref + aqSOAcon_ref
    SOAcon_new = OPOAcon_new + ASOAcon_new + BSOAcon_new + IVSOAcon_new + aqSOAcon_new

    # Set China region
    extent = [73.75, 135.0, 18.0, 53.5]

    # Set figure configulation
    plot = ['POA', 'SOA']
    diffvmax = [12.0, 5.0]
    nplot = len(plot)
    unit ='ug/'+u'm\u00B3'

    for pl in range(len(plot)):
        # Create a PDF file of the plots
        pdf = PdfPages( 'diff_' + plot[pl] + '.pdf')
        proj = ccrs.PlateCarree()
        figs = plt.figure(figsize=[12, 6])

        # read POA or SOA concentration and calculate the difference between 2020 and 2013
        spcon_new = eval(plot[pl] + 'con_new') * xmask
        spcon_ref = eval(plot[pl] + 'con_ref') * xmask
        diffcon = spcon_new.mean(dim='time') - spcon_ref.mean(dim='time')
        
        # draw a figure for the difference
        gcplot.single_panel(diffcon, vmin=-diffvmax[pl], vmax=diffvmax[pl], use_cmap_RdBu=True,
                            extent=extent, china_map=True, coastlines=False, unit=unit)
        plt.show()

        # Save this page to PDF
        pdf.savefig(figs)
        plt.close(figs)
        # Save the PDF file to disk
        pdf.close()

if __name__ == "__main__":
    main()
