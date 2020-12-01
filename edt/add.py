from modeldata import from_local_file
from utilities import get_ncfiles_in_dir
import log

def add_eke_to_ncfiles(input_dir,output_dir,filename_format='%Y%m',log_file='edt/add_eke.log'):
    ncfiles = get_ncfiles_in_dir(input_dir)
    for ncfile in ncfiles:
        input_path = input_dir+ncfile
        log.info(log_file,f'Loading netcdf file: {input_path}')
        modeldata = from_local_file(input_path)
        log.info(log_file,'Adding EKE to modeldata...')
        modeldata.add_eke()
        output_path = modeldata.get_output_path(output_dir,filename_format=filename_format)
        log.info(log_file,f'Writing modeldata to netcdf file: {output_path}')
        _ = modeldata.write_to_netcdf(output_dir,filename_format=filename_format)
        modeldata = None