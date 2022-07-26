__author__ = "Carlos Eduardo Sousa Lima"
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
import rtree
import pygeos
#%%
dir_ncdf = glob("{}/ncdf/*.nc".format(os.getcwd()))
dir_shp = "{}/shapes/09_oros_CE_SRTM2000_20181024.shp".format(os.getcwd())
#%%
ds = xr.open_dataset(dir_ncdf[0])
ds.coords["lon"] = ds.coords["lon"] - 180 #verificar se lon ta de 0 até 360, caso sim, executar
#%%
#alterar o nome das coords, caso necessário
coords = np.meshgrid(ds.variables["lon"].values, ds.variables["lat"].values)
lon = np.ravel(coords[0])
lat = np.ravel(coords[1])
df_coords = pd.DataFrame({"lon": lon, "lat": lat})
# %%
shp = gpd.read_file(dir_shp)
gdf_grid = gpd.GeoDataFrame(
    df_coords,
    geometry = gpd.points_from_xy(df_coords["lon"], df_coords["lat"]),
    crs = shp.crs)
# %%
gdf_ins = gpd.sjoin(gdf_grid, shp, how="left")
gdf_ins = gdf_ins.loc[~np.isnan(gdf_ins['index_right'])]
# %%
ts_df = pd.DataFrame({"time": ds["time"].values})
for i,j in zip(gdf_ins["lon"], gdf_ins["lat"]):
    ts = ds[ds.variable_id].sel(lon = i, lat = j)
    ts = ts.to_dataframe().reset_index()
    ### Organizar a saída
#%%
extents = shp.geometry.total_bounds
long_esq, lat_inf, long_dir, lat_sup = extents
fig, ax = plt.subplots()
#shp.plot(ax = ax)
gdf_grid.plot(ax = ax)
ax.set_xlim(long_esq, long_dir)
ax.set_ylim(lat_inf, lat_sup)
# %%
for i,j in zip(gdf_ins["lon"], gdf_ins["lat"]):
    print(i)
    print(j)
    print("_")

# %%
