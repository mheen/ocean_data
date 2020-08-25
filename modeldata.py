from utilities import get_variable_name,convert_time_to_datetime, convert_lon_360_to_180
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker

class Dimension:
    def __init__(self, name: str, values: np.ndarray, units: str):
        self.name = name
        self.values = values
        self.units = units
        self.length = len(values)
        if name == 'time':
            self.datetime = convert_time_to_datetime(self.values,self.units)

class Quantity:
    def __init__(self, name: str, values: np.ndarray, units: str):
        self.values = values
        self.units = units
        self.name = name

class Quantity3D(Quantity):
    def __init__(self, name: str, values: np.ndarray, units: str):
        super().__init__(name,values,units)
        self.dimensions = ('time','lat','lon')

class Quantity4D(Quantity):
    def __init__(self, name: str, values: np.ndarray, units: str):
        super().__init__(name,values,units)
        self.dimensions = ('time','depth','lat','lon')

def create_depth_dimension(value=0,units='m'):
    return Dimension('depth',np.array([value]),units)

def netcdf_to_dimension(netcdf: Dataset, variable_name: str, new_name=None, i_use=None) -> Dimension:
    i_use = _all_slice_if_none(i_use)
    values = netcdf[variable_name][i_use].filled(fill_value=np.nan)
    if len(values.shape) == 2:
        unique_di = np.unique(np.round(np.diff(values[:,0]),3))
        unique_dj = np.unique(np.round(np.diff(values[0,:]),3))
        if len(dj) == 1 and unique_dj[0] == 0:
            values = values[:,0]
        elif len(di) == 1 and unique_di[0] == 0:
            values = values[0,:]
    try:
        units = netcdf[variable_name].units
    except:
        units = ''
    new_name = _new_name_is_variable_name_if_none(new_name,variable_name)
    return Dimension(new_name,values,units)

def netcdf_to_quantity(netcdf: Dataset, variable_name: str, new_name=None,
                       i_times=None, i_depths=None, i_lats=None, i_lons=None):
    if len(netcdf[variable_name].shape) == 3:
        return netcfd_to_quantity3D(netcdf,variable_name,new_name=new_name,
                                    i_times=i_times,i_lats=i_lats,i_lons=i_lons)
    if len(netcdf[variable_name].shape) == 4:
        return netcdf_to_quantity4D(netcdf,variable_name,new_name=new_name,
                                    i_times=i_times,i_depths=i_depths,
                                    i_lats=i_lats,i_lons=i_lons)

def netcfd_to_quantity3D(netcfd: Dataset, variable_name: str, new_name=None,
                         i_times=None, i_lats=None, i_lons=None) -> Quantity3D:
    i_times = _all_slice_if_none(i_times)
    i_lats = _all_slice_if_none(i_lats)
    i_lons = _all_slice_if_none(i_lons)
    values = netcfd[variable_name][i_times,i_lats,i_lons].filled(fill_value=np.nan)
    units = netcfd[variable_name].units
    new_name = _new_name_is_variable_name_if_none(new_name,variable_name)
    return Quantity3D(new_name,values,units)

def netcdf_to_quantity4D(netcdf : Dataset, variable_name : str, new_name=None,
                         i_times=None, i_depths=None,
                         i_lats=None, i_lons=None) -> Quantity4D:
    i_times = _all_slice_if_none(i_times)
    i_depths = _all_slice_if_none(i_depths)
    i_lats = _all_slice_if_none(i_lats)
    i_lons = _all_slice_if_none(i_lons)
    values = netcdf[variable_name][i_times,i_depths,i_lats,i_lons].filled(fill_value=np.nan)
    units = netcdf[variable_name].units
    new_name = _new_name_is_variable_name_if_none(new_name,variable_name)
    return Quantity4D(new_name,values,units)

def _new_name_is_variable_name_if_none(new_name,variable_name):
    if new_name is None:
        new_name = variable_name
    return new_name

def _all_slice_if_none(i):
    if i is not None:
        return i
    return slice(None,None,None)

class ModelData:
    def __init__(self,
        time: Dimension,
        depth: Dimension,
        lat: Dimension,
        lon: Dimension,
        u=None,
        v=None,
        w=None,
        temp=None,
        salinity=None,
        ssh=None,
        mld=None,
        sea_ice_cover=None,
        o2=None,
        u_stokes=None,
        v_stokes=None,
        Hs=None,
        Hs_dir=None,
        Tp=None,
        Tm=None,
        Hs_sea=None,
        Hs_sea_dir=None,
        Tm_sea=None,
        Hs_swell=None,
        Hs_swell_dir=None,
        Tm_swell=None,
        u10=None,
        v10=None,
        temp2m=None):
        self.time = time
        self.depth = depth
        self.lat = lat
        self.lon = lon
        self.u = u
        self.v = v
        self.w = w
        self.temp = temp
        self.salinity = salinity
        self.ssh = ssh
        self.mld = mld
        self.sea_ice_cover = sea_ice_cover
        self.o2 = o2        
        self.u_stokes = u_stokes
        self.v_stokes = v_stokes
        self.Hs = Hs
        self.Hs_dir = Hs_dir
        self.Tp = Tp
        self.Tm = Tm
        self.Hs_sea = Hs_sea
        self.Hs_sea_dir = Hs_sea_dir
        self.Tm_sea = Tm_sea
        self.Hs_swell = Hs_swell
        self.Hs_swell_dir = Hs_swell_dir
        self.Tm_swell = Tm_swell
        self.u10 = u10
        self.v10 = v10
        self.temp2m = temp2m

    def fill_variable(self,variable_name,variable):
        if hasattr(self,variable_name):
            setattr(self,variable_name,variable)
        else:
            raise ValueError(f'ModelData does not have {variable_name} attribute, failed to fill variable.')

    def convert_lon360_to_lon180(self):
        if self.lon.values.min() >= 0:
            self.lon.values,i_lon = convert_lon_360_to_180(self.lon.values)
            variable_names = self.get_variable_names()
            for variable_name in variable_names:
                variable = getattr(self,variable_name)
                if variable is None:
                    continue
                if len(variable.dimensions) == 3:
                    variable.values = variable.values[:,:,i_lon]
                if len(variable.dimensions) == 4:
                    variable.values = variable.values[:,:,:,i_lon]
                setattr(self,variable_name,variable)

    def sort_lat_ascending(self):
        if self.lat.values[0] > self.lat.values[-1]:
            i_lat = np.argsort(self.lat.values)
            self.lat.values = self.lat.values[i_lat]
            variable_names = self.get_variable_names()
            for variable_name in variable_names:
                variable = getattr(self,variable_name)
                if variable is None:
                    continue
                if len(variable.dimensions) == 3:
                    variable.values = variable.values[:,i_lat,:]
                if len(variable.dimensions) == 4:
                    variable.values = variable.values[:,:,i_lat,:]
                setattr(self,variable_name,variable)

    def plot_variable(self,variable_name,i_time=0,i_depth=0,output_path=None):
        variable = getattr(self,variable_name)
        if variable is None:
            raise ValueError(f'Requested variable {variable_name} to plot is empty in ModelData.')
        if len(variable.dimensions) == 3:
            plot_values = variable.values[i_time,:,:]
        if len(variable.dimensions) == 4:
            plot_values = variable.values[i_time,i_depth,:,:]
        mpl.rc('font',family='arial')
        fig = plt.figure(figsize=(12,8))
        ax = plt.gca(projection=ccrs.PlateCarree())
        c = ax.contourf(self.lon.values,self.lat.values,plot_values,transform=ccrs.PlateCarree())
        ax.coastlines()
        # colorbar
        pos = ax.get_position()
        cbax = fig.add_axes([pos.x0+pos.width+0.01,pos.y0,0.02,pos.height])
        cbar = fig.colorbar(c,label=f'{variable.name} [{variable.units}]',cax=cbax)
        # grid
        parallels = np.arange(-80,90,20)
        meridians = np.arange(-180,190,60)
        lon_formatter = cticker.LongitudeFormatter()
        lat_formatter = cticker.LatitudeFormatter()
        ax.set_xticks(meridians,crs=ccrs.PlateCarree())
        ax.set_xticklabels(meridians)
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.set_yticks(parallels,crs=ccrs.PlateCarree())
        ax.set_yticklabels(parallels)
        ax.yaxis.set_major_formatter(lat_formatter)        
        ax.set_frame_on(True)
        if output_path is not None:
            plt.savefig(output_path,bbox_inches='tight',dpi=300)
        plt.show()

    def get_variable_names(self):
        variable_names = list(set(self.__dict__.keys())-set(['time','depth','lat','lon']))
        return variable_names

    def get_output_path(self,output_dir,filename_format='%Y%m%d'):
        time_string = self.time.datetime[0].strftime(filename_format)
        output_path = output_dir+time_string+'.nc'
        return output_path

    def write_to_netcdf(self,output_dir,filename_format='%Y%m%d'):
        output_path = self.get_output_path(output_dir,filename_format=filename_format)
        nc = Dataset(output_path,'w',format='NETCDF4')
        # --- define dimensions ---
        nc.createDimension(self.time.name,self.time.length)
        nc.createDimension(self.depth.name,self.depth.length)
        nc.createDimension(self.lat.name,self.lat.length)
        nc.createDimension(self.lon.name,self.lon.length)
        # --- write dimensions ---
        nc_time = nc.createVariable(self.time.name,float,self.time.name,zlib=True)
        nc_time[:] = self.time.values
        nc_time.units = self.time.units
        nc_depth = nc.createVariable(self.depth.name,float,self.depth.name,zlib=True)
        nc_depth[:] = self.depth.values
        nc_depth.units = self.depth.units
        nc_lat = nc.createVariable(self.lat.name,float,self.lat.name,zlib=True)
        nc_lat[:] = self.lat.values
        nc_lat.units = self.lat.units
        nc_lon = nc.createVariable(self.lon.name,float,self.lon.name,zlib=True)
        nc_lon[:] = self.lon.values
        nc_lon.units = self.lon.units
        # --- define and write variables ---
        variable_names = self.get_variable_names()
        for variable_name in variable_names:
            variable = getattr(self,variable_name)
            if variable is None:
                continue
            nc_var = nc.createVariable(variable.name,float,variable.dimensions,zlib=True)
            nc_var[:] = variable.values
            nc_var.units = variable.units
            nc_var = None
        nc.close()
        return output_path
        
def from_downloaded(netcdf : Dataset, variables : list, model_name : str,
                    i_times=None, i_depths=None, i_lats=None, i_lons=None, depth_value=None) -> ModelData:
    time_name = get_variable_name(model_name,'time')
    time = netcdf_to_dimension(netcdf,time_name,new_name='time',i_use=i_times)
    depth_name = get_variable_name(model_name,'depth')
    if depth_name in netcdf.dimensions.keys():
        depth = netcdf_to_dimension(netcdf,depth_name,new_name='depth',i_use=i_depths)
    else:
        depth = create_depth_dimension(value=depth_value)
    lat_name = get_variable_name(model_name,'lat')
    lat = netcdf_to_dimension(netcdf,lat_name,new_name='lat',i_use=None)
    lon_name = get_variable_name(model_name,'lon')
    lon = netcdf_to_dimension(netcdf,lon_name,new_name='lon',i_use=i_lons)
    modeldata = ModelData(time,depth,lat,lon)
    for variable_name in variables:
        variable_name_model = get_variable_name(model_name,variable_name)
        variable = netcdf_to_quantity(netcdf,variable_name_model,new_name=variable_name,
                                        i_times=i_times,i_depths=i_depths,i_lats=i_lats,i_lons=i_lons)
        modeldata.fill_variable(variable_name,variable)
    return modeldata

def from_local_file(input_path : str, variables=None,
                    i_times=None, i_depths=None, i_lats=None, i_lons=None, depth_value=None) -> ModelData:
    netcdf = Dataset(input_path)
    time = netcdf_to_dimension(netcdf,'time',i_times)
    if 'depth' in netcdf.dimensions.keys():
        depth = netcdf_to_dimension(netcdf,'depth',i_depths)
    else:
        depth = create_depth_dimension(value=depth_value)
    lat = netcdf_to_dimension(netcdf,'lat',i_lats)
    lon = netcdf_to_dimension(netcdf,'lon',i_lons)
    modeldata = ModelData(time,depth,lat,lon)
    if variables is None:
        variables = list(set(netcdf.variables.keys())-set(['time','depth','lat','lon']))
    for variable_name in variables:
        variable = netcdf_to_quantity(netcdf,variable_name,
                                        i_times=i_times,i_depths=i_depths,i_lats=i_lats,i_lons=i_lons)
        modeldata.fill_variable(variable_name,variable)
    netcdf.close()
    return modeldata
