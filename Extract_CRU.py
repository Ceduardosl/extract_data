#%%
import rioxarray as rxr
import xarray as xr
import numpy as np
import pandas as pd
import os
import geopandas as gpd
import matplotlib.pyplot as plt
import pyproj
from glob import glob
import netCDF4 as nc
import scipy as sc
from shapely.geometry import point
#%%
basins = ["Redencao_WGS84"]

list_nc = glob("{}/Data/*.nc".format(os.getcwd()))

for basin in basins:

    dir_shp = "{}/Shapes/{}.shp".format(os.getcwd(), basin)
    shp = gpd.read_file(dir_shp)

    for path_nc in list_nc:

        ds = rxr.open_rasterio(path_nc, engine = "netcdf4", masked=True)
        ds.rio.write_crs("EPSG:4326", inplace = True)
        name = path_nc.split(".")[5]

        ds_clip = ds.rio.clip(shp['geometry'], all_touched = True)
        df = ds_clip.to_dataframe().reset_index()
        #df.drop(["spatial_ref", "stn"], axis = 1, inplace = True)

        df.to_csv("{}/Extracted/CRU_{}.csv".format(os.getcwd(), name), sep = ";", index = False)
#%%