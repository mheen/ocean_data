from utilities import get_ncfiles_in_dir
from modeldata import ModelData, Dimension, Quantity3D, Quantity4D
from modeldata import from_local_file as modeldata_from_local_file
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import os
import log

def all_files_in_dir_horizontally(input_dir : str, output_dir : str,
                                  dx=1, x0=0.5, filename_format='%Y%m', method='linear',
                                  log_file='edt/interpolate_horizontally.log'):
    lon_interp = np.arange(-180.+x0,180+x0,dx)
    lat_interp = np.arange(-80+x0,90+x0,dx)
    filenames = get_ncfiles_in_dir(input_dir)
    for filename in filenames:
        modeldata = modeldata_from_local_file(input_dir+filename)
        modeldata.sort_lat_ascending()
        output_path = modeldata.get_output_path(output_dir)
        if os.path.exists(output_path):
            log.info(log_file,f'File already exists, skipping: {output_path}')
            continue
        log.info(log_file,f'Horizontally interpolating to {str(dx)} resolution, saving file: {output_path}')
        modeldata_interp = horizontally(modeldata,lat_interp,lon_interp,method=method)
        modeldata_interp.write_to_netcdf(output_dir,filename_format=filename_format)

def horizontally(modeldata : ModelData, lat_interp : np.ndarray, lon_interp : np.ndarray,
                 method='linear') -> ModelData:
    lat_interp_dimension = Dimension('lat',lat_interp,'degrees_north')
    lon_interp_dimension = Dimension('lon',lon_interp,'degrees_east')
    modeldata_interp = ModelData(modeldata.time,modeldata.depth,lat_interp_dimension,lon_interp_dimension)
    lat_interp2d,lon_interp2d = np.meshgrid(lat_interp,lon_interp)
    lat_interp2d = lat_interp2d.transpose()
    lon_interp2d = lon_interp2d.transpose()
    variable_names = modeldata.get_variable_names()
    if modeldata.lon.values[0] > -179.5:
        lon = np.insert(modeldata.lon.values,0,-modeldata.lon.values[-1])        
        boundary = 'periodic'
    else:
        lon = modeldata.lon.values
        boundary = None
    for variable_name in variable_names:
        variable = getattr(modeldata,variable_name)
        variable_interp = _horizontally_interpolate_quantity(modeldata.lat.values,lon,
                                           lat_interp2d,lon_interp2d,variable,method=method,boundary=boundary)
        modeldata_interp.fill_variable(variable_name,variable_interp)
    return modeldata_interp

def _horizontally_interpolate_quantity(lat, lon, lat2, lon2, quantity, method='linear',boundary=None):
    if quantity is None:
        return None
    if boundary == 'periodic':
        values = _add_periodic_boundary_to_values(quantity.values,quantity.dimensions)
    else:
        values = quantity.values
    quantity_shape = values.shape
    if len(quantity.dimensions) == 3:
        interp_values = np.empty((quantity_shape[0],lat2.shape[0],lon2.shape[1]))*np.nan
        for t in range(interp_values.shape[0]):
            f_interp = RegularGridInterpolator((lat,lon),values[t,:,:],method=method)
            interp_values[t,:,:] = f_interp((lat2,lon2))
        return Quantity3D(quantity.name,interp_values,quantity.units)
    if len(quantity.dimensions) == 4:
        interp_values = np.empty((quantity_shape[0],quantity_shape[1],lat2.shape[0],lon2.shape[1]))*np.nan
        for t in range(interp_values.shape[0]):
            for d in range(interp_values.shape[1]):
                f_interp = RegularGridInterpolator((lat,lon),values[t,d,:,:],method=method)
                interp_values[t,d,:,:] = f_interp((lat2,lon2))
        return Quantity4D(quantity.name,interp_values,quantity.units)

def _add_periodic_boundary_to_values(values,dimensions):
    if len(dimensions) == 3:
        insert_values = values[:,:,-1]
        values_periodic = np.insert(values,0,insert_values,axis=2)
    if len(dimensions) == 4:
        insert_values = values[:,:,:,-1]
        values_periodic = np.insert(values,0,insert_values,axis=3)
    return values_periodic    
