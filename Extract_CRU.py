#%%
import pandas as pd
import geopandas as gpd
import os
import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
from glob import glob

list_nc = glob("CRU/*.nc")

for dir_nc in list_nc:
    
    var = dir_nc.split(".")[-3]
    nc_data = xr.open_dataset(dir_nc)

    lon, lat = np.meshgrid(nc_data["lon"], nc_data["lat"])
    df_point = pd.DataFrame({"lon": lon.flatten(), "lat": lat.flatten()})
    shp_point = gpd.GeoDataFrame(df_point, geometry = gpd.points_from_xy(df_point["lon"], df_point["lat"], crs="EPSG:4326"))

    basins = glob("shapes/*.shp")

    for basin in basins:
        shp_basin = gpd.read_file(basin)

        gdf_result = gpd.sjoin(shp_point, shp_basin, how = "left")
        gdf_result = gdf_result.loc[~np.isnan(gdf_result["index_right"])]
        ins_point = gdf_result[["lon", "lat"]].values[:]
        ins_df = pd.DataFrame(ins_point, columns = ["lon", "lat"])

        var_ins = []

        for i in range(len(ins_df)):
            var_ins.append(nc_data[var].sel(lat = ins_df["lat"][i], lon = ins_df["lon"][i]).values)

        var_df = ins_df[["lon", "lat"]].join(pd.DataFrame(var_ins, columns = nc_data["time"]))
        var_df.to_csv("Extracted/{}_{}.csv".format(basin.split("\\")[1].split(".")[0], var), sep = ";", index = None, header = True)
#%%