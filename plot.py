from utilities import get_dir
from modeldata import from_local_file
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cftr
import cartopy.mpl.ticker as cticker
import numpy as np
import cmocean.cm as cmo

def _cartopy_map_and_grid(ax,lon_range=[-180,180],lat_range=[-80,90],dx=60,dy=20,resolution=None):
    if resolution is None:
        ax.add_feature(cftr.LAND,facecolor='#A3A3A3',edgecolor='k',zorder=2)
    else:
        ax.coastlines(resolution=resolution, color='k', linewidth=1, zorder=2)
    ax.set_extent([lon_range[0],lon_range[1],lat_range[0],lat_range[1]],ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    lat_formatter = cticker.LatitudeFormatter()
    meridians = np.arange(lon_range[0],lon_range[1],dx)
    parallels = np.arange(lat_range[0],lat_range[1],dy)
    ax.set_xticks(meridians,crs=ccrs.PlateCarree())
    ax.set_xticklabels(meridians)
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.set_yticks(parallels,crs=ccrs.PlateCarree())
    ax.set_yticklabels(parallels)        
    ax.yaxis.set_major_formatter(lat_formatter)    
    ax.grid(b=True,linewidth=2,color='k',linestyle='-')

def _pcolormesh_and_colorbar(fig,ax,lon,lat,variable_values,cbar_label,vmin=None,vmax=None,cmap='viridis'):
    c = ax.pcolormesh(lon,lat,variable_values,cmap=cmap,transform=ccrs.PlateCarree(),zorder=1,vmin=vmin,vmax=vmax)
    pos = ax.get_position()
    cbax = fig.add_axes([pos.x0+pos.width+0.01,pos.y0,0.02,pos.height])
    cbar = fig.colorbar(c,label=cbar_label,cax=cbax)

def _quiverplot(ax,lon,lat,u,v,thin=5,scale=25):
    if thin is not None:
        i = np.arange(0,u.shape[0],thin)
        j = np.arange(0,u.shape[1],thin)
        u = u[i][:,j]
        v = v[i][:,j]
        lon = lon[j]
        lat = lat[i]
    ax.quiver(lon,lat,u,v,scale=scale,transform=ccrs.PlateCarree(),zorder=3)

def _add_text_watermark(ax,wm_text='MegaMove',color='#EEEEEE',fontsize=130,alpha=0.5):
    ax.text(0.5,0.5,wm_text,fontsize=fontsize,color=color,
             ha='center',va='center',alpha=alpha)

def map_velocities(modeldata,cbar_label='Surface currents [m/s]',thin=5,scale=25,
                   lon_range=[112,118], lat_range=[-36,-25],
                   cmap=cmo.speed,output_path=None):
    lon = modeldata.lon.values
    lat = modeldata.lat.values
    u = np.copy(modeldata.u.values[0,:,:])
    v = np.copy(modeldata.v.values[0,:,:])
    vel = np.copy(modeldata.vel.values[0,:,:])
    u[u==0] = np.nan
    v[v==0] = np.nan
    vel[vel==0] = np.nan
    plt.style.use('input/plots.mplstyle')
    fig = plt.figure()
    ax = plt.gca(projection=ccrs.PlateCarree())
    _cartopy_map_and_grid(ax,lon_range=lon_range,lat_range=lat_range,dx=2,dy=2,resolution='50m')
    _pcolormesh_and_colorbar(fig,ax,lon,lat,vel,cbar_label,cmap=cmap,vmin=0,vmax=1.5)
    _quiverplot(ax,lon,lat,u,v,thin=thin,scale=scale)
    if output_path is not None:
        plt.savefig(output_path,bbox_inches='tight',dpi=300)
    plt.show()

def global_map_velocities(lon,lat,u,v,vel,cbar_label,thin=5,scale=25,
                          cmap='viridis',output_path=None):
    plt.style.use('input/plots.mplstyle')
    fig = plt.figure()
    ax = plt.gca(projection=ccrs.PlateCarree())
    _cartopy_map_and_grid(ax)
    _pcolormesh_and_colorbar(fig,ax,lon,lat,vel,cbar_label,cmap=cmap)
    _quiverplot(ax,lon,lat,u,v,thin=thin,scale=scale)
    _add_text_watermark(ax)
    if output_path is not None:
        plt.savefig(output_path,bbox_inches='tight',dpi=300)
    plt.show()

def global_map(lon,lat,variable_values,cbar_label,cmap='viridis',output_path=None):
    plt.style.use('input/plots.mplstyle')
    fig = plt.figure()
    ax = plt.gca(projection=ccrs.PlateCarree())
    _cartopy_map_and_grid(ax)
    _pcolormesh_and_colorbar(fig,ax,lon,lat,variable_values,cbar_label,cmap=cmap)
    _add_text_watermark(ax)
    if output_path is not None:
        plt.savefig(output_path,bbox_inches='tight',dpi=300)
    plt.show()

if __name__ == '__main__':
    cmap_style = 'normal'
    output_dir = get_dir('plot_output')
    modeldata = from_local_file(get_dir('plot_input'))
    lon = modeldata.lon.values
    lat = modeldata.lat.values
    # 2D variables
    variables = ['sst','temp2m','mld','o2','salinity','sea_ice_cover']
    cbar_labels = ['SST [$\degree$C]','Air temperature [$\degree$C]','Mixed layer depth [m]',
                   'Oxygen concentration [mmol/m$^3$]','Salinity [10$^{-3}^]','Sea-ice cover [fraction]']
    cmaps = [cmo.thermal,cmo.thermal,cmo.deep,cmo.tempo,cmo.haline,cmo.ice_r]
    for i in range(len(variables)):
        if len(getattr(modeldata,variables[i]).dimensions) == 3:
            values = getattr(modeldata,variables[i]).values[0,:,:]
        elif len(getattr(modeldata,variables[i]).dimensions) == 4:
            values = getattr(modeldata,variables[i]).values[0,0,:,:]
        if cmap_style == 'fancy':
            global_map(lon,lat,values,cbar_labels[i],cmap=cmaps[i],output_path=f'{output_dir}{variables[i]}.jpg')
        else:
            global_map(lon,lat,values,cbar_labels[i],output_path=f'{output_dir}{variables[i]}.jpg')
    # ocean currents
    u = modeldata.u.values[0,0,:,:]
    v = modeldata.v.values[0,0,:,:]
    vel = modeldata.vel.values[0,0,:,:]
    if cmap_style == 'fancy':
        global_map_velocities(lon,lat,u,v,vel,'Surface current velocity [m/s]',
                              cmap=cmo.speed,output_path=f'{output_dir}surface_currents.jpg')
    else:
        global_map_velocities(lon,lat,u,v,vel,'Surface current velocity [m/s]',
                              output_path=f'{output_dir}surface_currents.jpg')
    # wind
    u10 = modeldata.u10.values[0,:,:]
    v10 = modeldata.v10.values[0,:,:]
    vel10 = modeldata.vel10.values[0,:,:]
    if cmap_style == 'fancy':
        global_map_velocities(lon,lat,u10,v10,vel10,'Wind velocity [m/s]',thin=10,scale=350,
                              cmap=cmo.speed,output_path=f'{output_dir}wind_velocity.jpg')
    else:
        global_map_velocities(lon,lat,u10,v10,vel10,'Wind velocity [m/s]',thin=10,scale=350,
                              output_path=f'{output_dir}wind_velocity.jpg')
    