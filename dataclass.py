from utilities import get_variable_name,convert_time_to_datetime
from netCDF4 import Dataset
import numpy as np

class PhysData:
    def __init__(self,data,i_times,i_depths,variables,model_name):
        self.time = []
        self.time_units = []
        self.lon = []
        self.lon_units = []
        self.lat = []
        self.lat_units = []
        self.depth = []
        self.depth_units = []
        self.u = []
        self.u_units = []
        self.v = []
        self.v_units = []
        self.w = []
        self.w_units = []
        self.temp = []
        self.temp_units = []
        self.salinity = []
        self.salinity_units = []
        self.ssh = []
        self.ssh_units = []
        self.mld = []
        self.mld_units = []
        self.sea_ice_cover = []
        self.sea_ice_cover_units = []
        self.fill_dimensions(data,i_times,i_depths,model_name)
        self.fill_variables(data,i_times,i_depths,variables,model_name)

    def fill_dimensions(self,data,i_times,i_depths,model_name):
        model_time = get_variable_name(model_name,'time')
        self.time = data[model_time][i_times].filled(fill_value=np.nan)
        self.time_units = data[model_time].units
        model_lon = get_variable_name(model_name,'lon')
        self.lon = data[model_lon][:].filled(fill_value=np.nan)
        self.lon_units = data[model_lon].units
        model_lat = get_variable_name(model_name,'lat')
        self.lat = data[model_lat][:].filled(fill_value=np.nan)
        self.lat_units = data[model_lat].units
        model_depth = get_variable_name(model_name,'depth')
        self.depth = data[model_depth][i_depths].filled(fill_value=np.nan)
        if type(i_depths) == int:
            self.depth = [self.depth]
        self.depth_units = data[model_depth].units

    def fill_variables(self,data,i_times,i_depths,variables,model_name):
        for variable in variables:
            model_variable = get_variable_name(model_name,variable)
            if len(data[model_variable].shape) == 3:
                values = data[model_variable][i_times,:,:].filled(fill_value=np.nan)
            elif len(data[model_variable].shape) == 4:
                values = data[model_variable][i_times,i_depths,:,:].filled(fill_value=np.nan)
            units = data[model_variable].units
            setattr(self,variable,values)
            setattr(self,variable+'_units',units)

    def write_to_netcdf(self,output_dir):
        time = convert_time_to_datetime(self.time,self.time_units)
        time_string = time[0].strftime('%Y%m%d')
        output_path = output_dir+time_string+'.nc'
        nc = Dataset(output_path,'w',format='NETCDF4')
        # --- define dimensions ---
        nc.createDimension('time',len(self.time))
        nc.createDimension('lat',len(self.lat))
        nc.createDimension('lon',len(self.lon))
        nc.createDimension('depth',len(self.depth))
        variables_size = ('time','depth','lat','lon')
        # --- define and write variables ---
        # time
        nc_time = nc.createVariable('time',float,'time',zlib=True)
        nc_time[:] = self.time
        nc_time.units = self.time_units
        # lon
        nc_lon = nc.createVariable('lon',float,'lon',zlib=True)
        nc_lon[:] = self.lon
        nc_lon.units = self.lon_units
        # lat
        nc_lat = nc.createVariable('lat',float,'lat',zlib=True)
        nc_lat[:] = self.lat
        nc_lat.units = self.lat_units
        # depth
        nc_depth = nc.createVariable('depth',float,'depth',zlib=True)
        nc_depth[:] = self.depth
        nc_depth.units = self.depth_units
        # u-velocity
        if len(self.u) != 0:
            nc_u = nc.createVariable('u',float,variables_size,zlib=True)
            nc_u[:] = self.u
            nc_u.units = self.u_units
        # v-velocity
        if len(self.v) != 0:
            nc_v = nc.createVariable('v',float,variables_size,zlib=True)
            nc_v[:] = self.v
            nc_v.units = self.v_units
        # w-velocity
        if len(self.w) != 0:
            nc_w = nc.createVariable('w',float,variables_size,zlib=True)
            nc_w[:] = self.w
            nc_w.units = self.w_units
        # temperature
        if len(self.temp) != 0:
            nc_temp = nc.createVariable('temp',float,variables_size,zlib=True)
            nc_temp[:] = self.temp
            nc_temp.units = self.temp_units
        # salinity
        if len(self.salinity) != 0:
            nc_salinity = nc.createVariable('salinity',float,variables_size,zlib=True)
            nc_salinity[:] = self.salinity
            nc_salinity.units = self.salinity_units
        # ssh
        if len(self.ssh) != 0:
            nc_ssh = nc.createVariable('ssh',float,('time','lat','lon'),zlib=True)
            nc_ssh[:] = self.ssh
            nc_ssh.units = self.ssh_units
        # mld
        if len(self.mld) != 0:
            nc_mld = nc.createVariable('mld',float,variables_size,zlib=True)
            nc_mld[:] = self.mld
            nc_mld.units = self.mld_units
        nc.close()
        return output_path