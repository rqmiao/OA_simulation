#!/usr/bin/env python

# Imports
import os
import numpy as np
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import warnings
import gcpy.plot as gcplot
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import cartopy.io.shapereader as shapereader
import pandas as pd
from datetime import datetime

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# ----------------------------------------------------------------------
# Read information of sites
# ----------------------------------------------------------------------
input_path = '/home/mrq/pycode/'
sitedata = pd.read_csv(input_path + 'OC_network_sites.csv', usecols=['lon', 'lat'])
sitedata2 = pd.read_csv(input_path + 'OC_research_sites.csv', usecols=['lon', 'lat'])
nsites = len(sitedata.index)
nsites2 = len(sitedata2.index)

# ----------------------------------------------------------------------
# Create a PDF file of the plots
# ----------------------------------------------------------------------
# Define a PDF object so that we can save the plots to PDF
pdf = PdfPages('OC_site_location.pdf')
proj = ccrs.PlateCarree()
fig = plt.figure(figsize=[12, 6])

# ----------------------------------------------------------------------
# Read and display map of China
# ----------------------------------------------------------------------
# Read shapefile for China map
shp_path = '/home/mrq/pycode/mask/china_map/shapefiles/'
ax = plt.axes(projection=proj)
ax.set_extent([107, 125, 27, 43])  #NCP+YRD
provinces = cfeat.ShapelyFeature(
    shapereader.Reader(shp_path + 'China_provinces.shp').geometries(),
    proj, edgecolor='k', facecolor='none')
ax.add_feature(provinces, linewidth=0.4, zorder=2)
islands = cfeat.ShapelyFeature(
    shapereader.Reader(shp_path + 'bou2_4l.shp').geometries(),
    proj, edgecolor='k', facecolor='none')
ax.add_feature(islands, linewidth=0.4, zorder=2)
dashline = cfeat.ShapelyFeature(
    shapereader.Reader(shp_path + 'China_10-dash_line.shp').geometries(),
    proj, edgecolor='k', facecolor='none')
ax.add_feature(dashline, linewidth=0.4, zorder=2)
ax.add_feature(cfeat.LAND, facecolor='#FFFFF7')
ax.add_feature(cfeat.OCEAN, facecolor='#D4EBFE')

# Display subplot for South China Sea
fig = plt.gcf()
left, bottom, width, height = 0.75, 0.15, 0.08, 0.10
ax2 = fig.add_axes([left, bottom, width, height],projection=proj)
ax2.add_feature(provinces, linewidth=0.4, zorder=2)
ax2.add_feature(islands, linewidth=0.4, zorder=2)
ax2.add_feature(dashline, linewidth=0.4, zorder=2)
ax2.set_extent([105, 125, 4, 22])
ax2.add_feature(cfeat.LAND, facecolor='#FFFFF7')
ax2.add_feature(cfeat.OCEAN, facecolor='#D4EBFE')

# ----------------------------------------------------------------------
# Display sites by monitoring network and research sites
# ----------------------------------------------------------------------
for indexs in sitedata.index:
    site_coords = sitedata.loc[indexs].values
    # Find the indices corresponding to the site lon and lat
    lon = site_coords[1]
    lat = site_coords[0]
    ax.plot(lon, lat, 'h', color='#AFABAB', transform=ccrs.Geodetic(), ms=10.0, label='OC measurement sites')
for indexs in sitedata2.index:
    site_coords = sitedata2.loc[indexs].values
    # Find the indices corresponding to the site lon and lat
    lon = site_coords[1]
    lat = site_coords[0]
    ax.plot(lon, lat, 'h', color='#ffb299', transform=ccrs.Geodetic(), ms=10.0, label='OC measurement sites')

plt.show()

# -----------------------------
# Save this page to PDF
# -----------------------------
pdf.savefig(fig)

# ----------------------------------------------------------------------
# Save the PDF file to disk
# ----------------------------------------------------------------------
pdf.close()
