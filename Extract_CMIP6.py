#%%
import xarray as xr
import numpy as np
import pandas as pd 
import geopandas as gpd
import os
from glob import glob
import matplotlib.pyplot as plt
#%%
use_buffer = True
models = ["BCC-CSM2-MR", "CanESM5", "MIROC6", "MRI-ESM2-0"]
shp_base = gpd.read_file("{}/Shape_CE/Ceara_WGS84.shp".format(os.getcwd()))

for model in models:

    dir_ncdf = "{}/Climate_Data/{}".format(os.getcwd(), model)

    list_ncdf = glob("{}/*.nc".format(dir_ncdf))

    for path_ncdf in list_ncdf:

        scenario = path_ncdf.split("\\")[-1].split("_")[3]

        nc_data = xr.open_dataset(path_ncdf)

        nc_data.coords["lon"] = (nc_data.coords["lon"] + 180) % 360 - 180
        nc_data = nc_data.sortby("lon")

        if nc_data.nominal_resolution == "100 km":
            shp_buffer = gpd.GeoDataFrame(geometry = shp_base.buffer(1))
        if nc_data.nominal_resolution == "250 km":
            shp_buffer = gpd.GeoDataFrame(geometry = shp_base.buffer(1))
        if nc_data.nominal_resolution == "500 km":
            shp_buffer = gpd.GeoDataFrame(geometry = shp_base.buffer(2))
        
        if use_buffer:
            shp = shp_buffer
        else:
            shp = shp_base
            
        coords = np.meshgrid(nc_data.variables["lon"].values, nc_data.variables["lat"].values)
        lon = np.ravel(coords[0])
        lat = np.ravel(coords[1])
        df_coords = pd.DataFrame({"lon": lon, "lat": lat})

        gdf_grid = gpd.GeoDataFrame(
            df_coords,
            geometry = gpd.points_from_xy(df_coords.lon, df_coords.lat),
            crs = shp.crs
        )
        gdf_ins = gpd.sjoin(gdf_grid, shp, how = "left")
        gdf_ins = gdf_ins.loc[~np.isnan(gdf_ins.index_right)]
        coords_ins = gdf_ins[["lon", "lat"]]

        index_list = pd.Index(nc_data.time.values)
        index_list = index_list.insert(0, "lat")
        index_list = index_list.insert(0, "lon")
        ts_df = pd.DataFrame(index = index_list)

        count = 0
        for i, j in zip(coords_ins.lon, coords_ins.lat):
            ts = nc_data[nc_data.variable_id].sel(lon = i, lat = j).to_dataframe()[nc_data.variable_id]
            if (nc_data.variable_id == "tasmax") or (nc_data.variable_id == "tasmin"):
                ts = ts - 273.156
            if (nc_data.variable_id == "pr"):
                ts = ts*86400*30
            ts = pd.Series([i,j], index = ["lon", "lat"]).append(ts)
            ts_df.insert(len(ts_df.columns), count, ts)
            count += 1
        ts_df = ts_df.T
        ts_df.to_csv("{}/Climate_Data/{}/{}_{}_{}.csv".format(os.getcwd(), model, model, nc_data.variable_id, scenario), index = None, header = True)

        fig, ax = plt.subplots()
        shp.plot(ax = ax, edgecolor = "black", facecolor = "none")
        ax.scatter(coords_ins.lon, coords_ins.lat)
        ax.set_ylim(shp.geometry.bounds["miny"].values, shp.geometry.bounds["maxy"].values)
        ax.set_xlim(shp.geometry.bounds["minx"].values, shp.geometry.bounds["maxx"].values)
        ax.set_title("{} - {} - {}".format(model, nc_data.variable_id, scenario))
        ax.set_ylabel("latitude (°)")
        ax.set_xlabel("longitude (°)")
        fig.savefig("{}/Climate_Data/{}/{}_{}_{}.png".format(os.getcwd(), model, model, nc_data.variable_id, scenario),
                    dpi = 600, bbox_inches = "tight", facecolor = "w")

print("############# Finalizado!!!!! ##################")
#%%