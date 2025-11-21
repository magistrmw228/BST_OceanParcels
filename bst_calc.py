from datetime import timedelta as delta
from glob import glob
import os
from datetime import timedelta
import xarray as xr
import numpy as np
from scipy.spatial import KDTree
from haversine import haversine, Unit
from time import time as clock
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import parcels.rng as ParcelsRandom
from parcels import (
    AdvectionRK4_3D,
    FieldSet,
    JITParticle,
    ParticleSet,
    StatusCode,
)
import random

## The vertical velocity [m/s] can be changed here (negative for ascent)
# vertical_speed = - 0.01 # This can be changed only inside Buoyancy kernel
multi = 1 # how many times are there more particles with the same trajectories to calculate statistics

BASE_DIR = os.path.dirname(os.path.abspath(__file__))  # the directory of the current script
input_folder = os.path.join(BASE_DIR, 'Input')
output_folder = os.path.join(BASE_DIR, 'Output')
mesh_mask = os.path.join(input_folder, 'mesh_mask_new.nc')

ufiles = sorted(glob(f"{input_folder}/BAMS*_U_*.nc"))
vfiles = sorted(glob(f"{input_folder}/BAMS*_V_*.nc"))
wfiles = sorted(glob(f"{input_folder}/BAMS*_W_*.nc"))
tfiles = sorted(glob(f"{input_folder}/BAMS*_T_*.nc"))
mesh = xr.open_dataset(mesh_mask)

dt = 1 # calculations timestep in hours
dt_repeat = 24 # period of particle releases from the same locations
dt_output = 3 # output timestep in hours

## select geographic region of particle releases 
# (release area depth, latitude and longitude boundaries)
depth_start_min = 40
depth_start_max = 100
lat_start_min = 44
lat_start_max = 60
lon_start_min = 33
lon_start_max = 38


### Preprocessing start below ####
# Fixing the initial time
start_time = clock()
# create a folder for the output data, if it doesn't exist yet
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
## limiting the release area
gdepw_1d = mesh.gdepw_1d[:]
gdepw_1d = np.array(gdepw_1d[0])
mbathy = np.array(mesh.mbathy[:][0])
nav_lat = np.array(mesh.nav_lat[:])
nav_lon = np.array(mesh.nav_lon[:])
for j in range(mbathy.shape[0]):
    for i in range(mbathy.shape[1]):
        mbathy[j,i] = gdepw_1d[mbathy[j,i]]
index = (mbathy >= depth_start_min) & (mbathy <= depth_start_max) & (nav_lon >= lon_start_min) & (nav_lon <= lon_start_max) & (nav_lat >= lat_start_min) & (nav_lat <= lat_start_max)
lat_bottom = nav_lat[index]
lat_bottom = nav_lat[index]
lon_bottom = nav_lon[index]
depth_bottom = mbathy[index]
depth_bottom = np.astype(depth_bottom, np.float64)

## loading NEMO mesh
fmask = np.array(mesh.fmask[:])
fmask = np.squeeze(fmask)
fmask = fmask[0, :, :] 
land_mask_F = np.where(fmask == 0) 
glamt = np.array(mesh.glamt[:])
glamt = np.squeeze(glamt)
glamt = glamt[land_mask_F]
gphit = np.array(mesh.gphit[:])
gphit = np.squeeze(gphit)
gphit = gphit[land_mask_F]
glamu = np.array(mesh.glamu[:])
glamu = np.squeeze(glamu)
glamu = glamu[land_mask_F]
gphiu = np.array(mesh.gphiu[:])
gphiu = np.squeeze(gphiu)
gphiu = gphiu[land_mask_F]
glamv = np.array(mesh.glamv[:])
glamv = np.squeeze(glamv)
glamv = glamv[land_mask_F]
gphiv = np.array(mesh.gphiv[:])
gphiv = np.squeeze(gphiv)
gphiv = gphiv[land_mask_F]
lat_coast = np.vstack([gphit, gphiu, gphiv])
lon_coast = np.vstack([glamt, glamu, glamv])
land_mesh = np.column_stack((lat_coast.ravel(), lon_coast.ravel()))
# Creating a KD-tree to quickly find the nearest node
tree_mesh = KDTree(land_mesh)

# uploading data on distances to the coastline
GSHHG = xr.open_dataset(input_folder + "/dist_to_GSHHG_v2.3.7_1m.nc")
lon_min = 31
lon_max = 39
lat_min = 42
lat_max = 47
land_lat = np.array(GSHHG['lat'][:])
land_lon = np.array(GSHHG['lon'][:])
dist = np.array(GSHHG['dist'][:])
land_lon, land_lat = np.meshgrid(land_lon, land_lat)
dist = dist.ravel()
land_points = np.column_stack((land_lat.ravel(), land_lon.ravel()))
land_mask = np.where((land_points[:,0] >= lat_min) &
                                (land_points[:,0] <= lat_max) &
                                (land_points[:,1] >= lon_min)
                                & (land_points[:,1] <= lon_max))
land_points = land_points[land_mask]
dist = dist[land_mask]
tree_gshhg = KDTree(land_points)

## apply "beaching" threshold
for parcel in range(len(depth_bottom)):
    _, idx = tree_mesh.query((lat_bottom[parcel], lon_bottom[parcel]))
    nearest_node = land_mesh[idx]
    distance_mesh = haversine((lat_bottom[parcel], lon_bottom[parcel]), nearest_node, unit=Unit.KILOMETERS)
    _, idx = tree_gshhg.query((lat_bottom[parcel], lon_bottom[parcel]))
    distance_gsshg = -1 * dist[idx]
    distance = np.nanmin([distance_mesh, distance_gsshg])
    if distance < 0.5:
        lat_bottom[parcel] = np.nan
        lon_bottom[parcel] = np.nan
        depth_bottom[parcel] = np.nan
        
lat_bottom = lat_bottom[~np.isnan(lat_bottom)]
lon_bottom = lon_bottom[~np.isnan(lon_bottom)]
depth_bottom = depth_bottom[~np.isnan(depth_bottom)]

## loading NEMO files
filenames = {
    "U": {"lon": mesh_mask, "lat": mesh_mask, "depth": wfiles[0], "data": ufiles},
    "V": {"lon": mesh_mask, "lat": mesh_mask, "depth": wfiles[0], "data": vfiles},
    "W": {"lon": mesh_mask, "lat": mesh_mask, "depth": wfiles[0], "data": wfiles},
    "T": {"lon": mesh_mask, "lat": mesh_mask, "depth": wfiles[0], "data": tfiles},
}
variables = {"U": "uo", "V": "vo", "W": "wo", "T": "thetao"}

# Note that all variables need the same dimensions in a C-Grid
c_grid_dimensions = {
    "lon": "glamf",
    "lat": "gphif",
    "depth": "depthw",
    "time": "time_counter",
}
dimensions = {
    "U": c_grid_dimensions,
    "V": c_grid_dimensions,
    "W": c_grid_dimensions,
    "T": c_grid_dimensions,
}

fieldset = FieldSet.from_nemo(filenames, variables, dimensions)

## Set kernels to use
def CheckOutOfBounds(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
        particle.delete()
        
### Kernel for recording particle temperature
Kalkan = JITParticle.add_variable("temperature")
def SampleT(particle, fieldset, time):
    particle.temperature = fieldset.T[time, particle.depth, particle.lat, particle.lon]

# Kernel for positive buoyancy
def Buoyancy(particle, fieldset, time):
    # Case 1: the particle crossed the 0 m surface - we return it to 0
    if particle.state == StatusCode.ErrorThroughSurface:
        particle.depth = 0.0
        particle.state = StatusCode.Success
    else:
    # Case 2: The particle is below the surface layer of 0-2 m, - it is necessary to give it buoyancy
        if (particle.depth + particle_ddepth >= 1.0):
            # vertical_speed = - ParcelsRandom.uniform(0.0005, 0.001)
            vertical_speed = -0.01
            upper_depth = ParcelsRandom.uniform(0.0, 2.0)
            particle_ddepth = particle_ddepth + (vertical_speed * particle.dt)
            # If a particle crosses a randomly selected depth in the surface layer,
            # then we launch it to that depth.
            if particle.depth + particle_ddepth < upper_depth:
                particle_ddepth = upper_depth - particle.depth
    # Case 3: A particle in a 0-2 m layer, - it is necessary to remove the buoyancy
        elif particle.depth < 2.0:
            particle_ddepth = particle_ddepth
            
## Connecting the used data and kernels
pset = ParticleSet.from_list(
    fieldset=fieldset,
    pclass=Kalkan,
    depth=depth_bottom,
    lon=lon_bottom,
    lat=lat_bottom,
    repeatdt = delta(hours=dt_repeat)
)

kernels = pset.Kernel([AdvectionRK4_3D,
                    CheckOutOfBounds, 
                    Buoyancy, 
                    SampleT])

output_file = pset.ParticleFile(name=output_folder + '/Parcels_2001.zarr', 
                                outputdt=timedelta(hours=dt_output), 
                                chunks=(len(depth_bottom)*10, 10),
)

## Launching the model
# Starting the spawning period
pset.execute(kernels, 
            runtime=delta(days=75), 
            dt=delta(hours=dt),
            output_file=output_file,
)
# finish spawning
pset.repeatdt = None

# Continue tracking for another 14 days
pset.execute(kernels, 
            runtime=delta(days=14), 
            dt=delta(hours=dt),
            output_file=output_file,
)
print(f'Finish calculating trajectories in {(clock() - start_time) / 60:.1f} min. from start')


##### POSTPROCESSING ######
# Open the file and load variables
ds = xr.open_dataset(output_folder + '/Parcels_2001.zarr', 
                    engine="zarr")
lon = np.array(ds["lon"].values)
lat = np.array(ds["lat"].values)
z = np.array(ds["z"].values)
temp = np.array(ds["temperature"].values)
time = np.array(ds["time"].values)
if 'ds' in globals():
    del ds
    
# We make sure that the particles that have not yet started have NaN coordinates and time.
step = len(depth_bottom)
start = step
end = 0
time_step = int(dt_repeat/dt_output)
shear = 0
while end < lon.shape[0]:
    end = start + step
    shear = shear + time_step
    lon[start:end, shear:] = lon[start:end, 0:-shear]
    lon[start:end, :shear] = np.nan
    lat[start:end, shear:] = lat[start:end, 0:-shear]
    lat[start:end, :shear] = np.nan
    z[start:end, shear:] = z[start:end, 0:-shear]
    z[start:end, :shear] = np.nan
    temp[start:end, shear:] = temp[start:end, 0:-shear]
    temp[start:end, :shear] = np.nan
    time[start:end, shear:] = time[start:end, 0:-shear]
    time[start:end, :shear] = np.datetime64('NaT')
    start = end
    
# Spawning only in the layer from 7.5 to 12 degrees Celsius
n_threads = 10 # Number of threads

def check_row(row_index, lon, temp):
    not_nan_indices = np.where(~np.isnan(lon[row_index]))[0]
    if len(not_nan_indices) > 0:
        first = not_nan_indices[0]
        if (temp[row_index, first] < 7.5) or (temp[row_index, first] > 12):
            return False
    return True

def process_chunk(start, end, lon, temp):
    return [check_row(i, lon, temp) for i in range(start, end)]

rows_per_thread = temp.shape[0] // n_threads

mask = np.ones(temp.shape[0], dtype=bool)

with ThreadPoolExecutor(max_workers=n_threads) as executor:
    futures = []
    for i in range(n_threads):
        start = i * rows_per_thread
        end = (i + 1) * rows_per_thread if i < n_threads - 1 else temp.shape[0]
        futures.append(executor.submit(process_chunk, start, end, lon, temp))
    results = [future.result() for future in futures]
    
mask = np.concatenate(results)
lat = lat[mask]
lon = lon[mask]
z = z[mask]
time = time[mask]
temp = temp[mask]
on_land = np.zeros(lat.shape)


for row in range(temp.shape[0]):
    zero_index = np.argmax(temp[row] == 0)
    if temp[row, zero_index] == 0:
        temp[row, zero_index:] = np.nan
        lat[row, zero_index+1:] = np.nan
        lon[row, zero_index+1:] = np.nan
        z[row, zero_index+1:] = np.nan
        time[row, zero_index+1:] =  np.datetime64('NaT')
        on_land[row, zero_index:] = True
       
for row_index in range(lat.shape[0]):
    not_nan_indices = np.where(~np.isnan(lon[row_index]))[0]
    if len(not_nan_indices) > 0:
        first = not_nan_indices[0]
        last = not_nan_indices[-1]
        
        for timestep in range(first, last + 1):
            _, idx = tree_gshhg.query((lat[row_index, timestep], lon[row_index, timestep]))
            distance = dist[idx] * -1
            
            if distance < 0.5:
                lat[row_index, timestep+1:] = np.nan
                lon[row_index, timestep+1:] = np.nan
                z[row_index, timestep+1:] = np.nan
                temp[row_index, timestep:] = np.nan
                time[row_index, timestep+1:] = np.datetime64('NaT')
                on_land[row_index, timestep:] = True
                break
vars_to_delete = ['fmask', 'land_mask_F', 'glamt', 'gphit', 'glamu', 'gphiu', 'glamv', 'gphiv', 
    'lat_coast', 'lon_coast', 'distance_mesh', 'distance_gsshg', 'tree_gshhg', 'tree_mesh']
for var in vars_to_delete:
    if var in globals():
        del globals()[var]
        
# Now we have trajectories for those parcels that started in good temperature spawning conditions
## Write to NetCDF file
ds = xr.Dataset(
    {
        "temperature": (["parcels", "timesteps"], temp),
        "latitude": (["parcels", "timesteps"], lat),  
        "longitude": (["parcels", "timesteps"], lon),  
        "depth": (["parcels", "timesteps"], z),      
        "time": (["parcels", "timesteps"], time),      
        "on_land": (["parcels", "timesteps"], on_land),      
        
    },
    coords={
        "timesteps": np.arange(1, time.shape[1]+1, 1),          
        "parcels": np.arange(1, time.shape[0]+1, 1)
    }
)
ds.temperature.attrs["units"] = "degrees Celsius"
ds.latitude.attrs["units"] = "degrees_north"
ds.longitude.attrs["units"] = "degrees_east"
ds.depth.attrs["units"] = "meters"
ds.attrs = {
    "title": "OceanParcels simulation",
    "description": "This dataset contains information about trajectories of passive Black Sea Kalkan fish eggs released in NEMO fields (MHI RAS configuration) using OceanParcels",
    "creator": "Dmitriy Krasilnikov",
    "date_created": pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
    "contact": "krasilnikov.dmitrij@inbox.ru",
    "institution": "MHI RAS",
    "source": "Generated for science purposes"
}

ds.to_netcdf(output_folder + '/Trajectories_2001.nc')
    
# Delete "ds" to free up memory
if 'ds' in globals():
    del ds
    
## The block of development from the egg to the hatched larva under the influence of variable temperatures is presented below.
# Experimental data on fish eggs to hatched larcae development rate
x = np.array([13, 15, 17, 19, 21])
y = np.array([146.4, 115.2, 96.7, 85.6, 78.3])
# Approximation by 3d-order polynomial function
order = 3
coefficients = np.polyfit(x, y, order)
poly_func = np.poly1d(coefficients)
# Function used later to calculate probability of egg to die on each time step
def P_temp(p_zero):
    return random.choices([0, 1], weights=[p_zero, 1 - p_zero])[0]   

# Single particle processing function
def process_parcel(parcel): # 1
    original_index = parcel % N
    not_nan_indices = np.where(~np.isnan(lon[original_index]))[0]
    
    if len(not_nan_indices) > 0: 
        first = not_nan_indices[0]
        last = not_nan_indices[-1]
        
        stage = 0  # Stage of development variable
        is_recorded = False
        # gastrulation starts when "stage" = 0.2 (20%)
        gastrulation_started = False
        final_data = None
        is_dead = False
        why_dead = None
        for timestep in range(first, last + 1, 1):
            if ((stage < 1) and (not final_data)):
                growth = 1 / (poly_func(temp[original_index, timestep])) * dt_output
                stage += growth
                # if this is the last step, and the egg has not reached its full development (stage = 1), then death
                if ((stage < 1) and (timestep == last)):
                    is_dead = True 
                    if on_land[original_index, timestep+1]:
                        why_dead = 'on_land'
                    else:
                        why_dead = 'simulation_end'
                # if the temperature is < 10 C at the beginning of gastrulation, then death
                elif ((stage > 0.2) and (temp[original_index, timestep] < 10.0)):
                    is_dead = True
                    why_dead = 'low_temp_gastrulation'   
                elif ((stage < 1) and (temp[original_index, timestep] >= 10.0) and (temp[original_index, timestep] <= 16.0) 
                    and (timestep < last) and (stage > 0.2)):
                    # Survival probability equation for the range of 10-16 C
                    is_dead = P_temp(0.1417 * temp[original_index, timestep] - 1.3667)  
                elif ((stage < 1) and (temp[original_index, timestep] > 16.0) and (temp[original_index, timestep] <= 18.0)
                    and (timestep < last) and (stage > 0.2)):
                    # The probability of survival equation (const) for the range of 16-18 C 
                    is_dead = P_temp(0.9)
                elif ((stage < 1) and (temp[original_index, timestep] > 18.0) and (temp[original_index, timestep] <= 21.0)
                    and (timestep < last) and (stage > 0.2)):
                    # Survival probability equation for the range of 18-22 C
                    is_dead = P_temp(-0.0667 * temp[original_index, timestep] + 2.1)
                else:
                    is_dead = False
                
                if is_dead:
                    if not why_dead:
                        why_dead = 'low_temp_after_gastrulation'
                        gastrulation_started = True
                                                    
                    final_data = [lat[original_index, first], lon[original_index, first], time[original_index, first], 
                                z[original_index, first], temp[original_index, first], 
                                lat[original_index, timestep], lon[original_index, timestep], time[original_index, timestep], 
                                z[original_index, timestep], temp[original_index, timestep], int(gastrulation_started), np.min([np.round(stage, 2), 0.99]), why_dead]
                    is_recorded = True    
            else:
                if not is_recorded:
                    gastrulation_started = True
                    final_data = [lat[original_index, first], lon[original_index, first], time[original_index, first], 
                                z[original_index, first], temp[original_index, first], 
                                lat[original_index, timestep], lon[original_index, timestep], time[original_index, timestep], 
                                z[original_index, timestep], temp[original_index, timestep], int(gastrulation_started), 1, why_dead]
                    is_recorded = True
                break
        
        if final_data:
            return final_data
        else:
            return None

N = temp.shape[0]

## Running parallel processing using threads
with ThreadPoolExecutor(max_workers=14) as executor:
    results = list(executor.map(process_parcel, range(N * multi)))
results = [r for r in results if r is not None]


## Converting the results to a DataFrame and then saving as NetCDF
columns = ["lat_start", "lon_start", "time_start", "depth_start", "temp_start", 
        "lat_finish", "lon_finish", "time_finish", "depth_finish", "temp_finish", "gastrulation", "stage", "why_dead"]
df = pd.DataFrame(results, columns=columns)
ds = xr.Dataset.from_dataframe(df)
ds.to_netcdf(output_folder + '/Start_and_Finish_2001.nc')
    
print(f'Finish postprocessing in {(clock() - start_time) / 60:.1f} min. from start')
print(f"Data stored in {output_folder}")
print('Trajectories_2001.nc contains trajectories for those parcels that started in good temperature spawning conditions')
print('Start_and_Finish_2001.nc contains final table')

