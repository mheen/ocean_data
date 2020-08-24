from datetime import datetime,timedelta
import numpy as np
import json
import os

# -----------------------------------------------
# General
# -----------------------------------------------
def get_dir(dirname,json_file='input/dirs.json'):
    with open(json_file,'r') as f:
        all_dirs = json.load(f)
    return all_dirs[dirname]

def get_variable_name(model,variable,json_file='input/variables.json'):
    with open(json_file,'r') as f:
        all_models = json.load(f)
        model = all_models[model]
    return model[variable]

def get_ecmwf_variable_code(variable,json_file='input/ecmwf_codes.json'):
    with open(json_file,'r') as f:
        all_codes = json.load(f)
    return all_codes[variable]

def get_variable_name_reverse(model,variable,json_file='input/variables.json'):
    with open(json_file,'r') as f:
        all_models = json.load(f)
        model = all_models[model]
    model_variable_names = list(model.values())
    variable_names = list(model.keys())
    i_variable = model_variable_names.index(variable)
    return variable_names[i_variable]

def get_urls(model,json_file='input/urls.json'):
    with open(json_file,'r') as f:
        all_urls = json.load(f)
    return all_urls[model]

def get_logins(model,json_file='input/logins.json'):
    with open(json_file,'r') as f:
        all_logins = json.load(f)
    return all_logins[model]

def get_ncfiles_in_dir(input_dir):
    ncfiles = []
    for filename in os.listdir(input_dir):
        if filename.endswith('.nc'):
            ncfiles.append(filename)
    return ncfiles

def get_ncfiles_in_time_range(input_dir,start_date,end_date,timeformat='%Y%m%d'):
        all_ncfiles = get_ncfiles_in_dir(input_dir)
        ndays = (end_date-start_date).days+1
        ncfiles = []
        for n in range(ndays):
            date = start_date+timedelta(days=n)
            for ncfile in all_ncfiles:
                if ncfile.startswith(date.strftime(timeformat)):
                    ncfiles.append(ncfile)
        return ncfiles

def get_closest_index(A,target):
    # A must be sorted!
    idx = A.searchsorted(target)
    idx = np.clip(idx,1,len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target-left < right-target
    return idx

def rename_ncfiles_in_dir(input_dir: str, filename_indices: list, filename_format: str, new_filename_format: str):
    ncfiles = get_ncfiles_in_dir(input_dir)
    for ncfile in ncfiles:
        input_path = input_dir+ncfile
        ncfile_date = datetime.strptime(ncfile[filename_indices[0]:filename_indices[1]],filename_format)
        output_path = input_dir+ncfile_date.strftime(new_filename_format)+'.nc'
        os.rename(input_path,output_path)        

# -----------------------------------------------
# Timeseries
# -----------------------------------------------
def get_time_index(time_array,time):
    '''Returns exact index of a requested time, raises
    error if this does not exist.'''
    t = np.where(time_array==time)[0]
    if len(t) > 1:
        raise ValueError('Multiple times found in time array that equal requested time.')
    elif len(t) == 0:
        raise ValueError('Requested time not found in time array.')
    else:
        return t[0]

def get_time_indices(timeseries,time):
    i_times = []
    for i,t in enumerate(timeseries):
        if t.date() == time.date():
            i_times.append(i)
    return i_times

def get_closest_time_index(time_array,time):
    '''Returns exact index of a requested time if is exists,
    otherwise returns the index of the closest time.'''
    dt = abs(time_array-time)
    i_closest = np.where(dt == dt.min())[0][0]
    return i_closest

def get_l_time_range(time,start_time,end_time):
    if type(start_time) is datetime.date:
        start_time = datetime.datetime(start_time.year,start_time.month,start_time.day)
    if type(end_time) is datetime.date:
        end_time = datetime.datetime(end_time.year,end_time.month,end_time.day)
    l_start = time >= start_time
    l_end = time <= end_time
    l_time = l_start & l_end
    return l_time

def get_n_months(start_date,end_date):
    n_months = end_date.month-start_date.month
    n_years = end_date.year-start_date.year    
    if not n_years == 0:
       n_months = n_months+12*n_years
    return n_months

def add_month_to_timestamp(timestamp,n_month):
    month = timestamp.month - 1 + n_month
    year = timestamp.year + month // 12
    month = month % 12 + 1
    return datetime(year,month,timestamp.day)

def convert_time_to_datetime(time_org,time_units):
    time = []
    i_start_time = time_units.index('since')+len('since')+1
    if 'T' in time_units: # YYYY-mm-ddTHH:MM format used by Parcels
        i_end_time = i_start_time+len('YYYY-mm-ddTHH:MM')
        base_time = datetime.strptime(time_units[i_start_time:i_end_time],'%Y-%m-%dT%H:%M')
    else: # YYYY-mm-dd format used by multiple numerical models
        i_end_time = i_start_time+len('YYYY-mm-dd')
        base_time = datetime.strptime(time_units[i_start_time:i_end_time],'%Y-%m-%d')
    if time_units.startswith('seconds'):
        if time_org.shape == ():
            time = base_time+timedelta(seconds=float(time_org))
            return time
        for t in time_org:
            if not np.isnan(t):
                time.append(base_time+timedelta(seconds=float(t)))
            else:
                time.append(np.nan)
        return np.array(time)
    elif time_units.startswith('hours'):
        if time_org.shape == ():
            time = base_time+timedelta(hours=float(time_org))
            return time
        for t in time_org:
            if not np.isnan(t):
                time.append(base_time+timedelta(hours=float(t)))
            else:
                time.append(np.nan)
        return np.array(time)
    elif time_units.startswith('days'):
        if time_org.shape == ():
            time = base_time+timedelta(days=float(time_org))
            return time
        for t in time_org:
            if not np.isnan(t):
                time.append(base_time+timedelta(days=float(t)))
            else:
                time.append(np.nan)
        return np.array(time)
    else:
        raise ValueError('Unknown time units for time conversion to datetime.')

def convert_datetime_to_time(time_org,time_units='seconds',time_origin=datetime(1995,1,1,12,0)):
    time = []
    if time_units == 'seconds':
        conversion = 1
    elif time_units == 'hours':
        conversion = 60*60
    elif time_units == 'days':
        conversion = 24*60*60
    else:
        raise ValueError('Unknown time units requested fro time conversion from datetime.')
    for t in time_org:        
        time.append((t-time_origin).total_seconds()/conversion)
    return np.array(time)

# -----------------------------------------------
# Coordinates
# -----------------------------------------------
def convert_lon_360_to_180(lon):
    lon180 = np.copy(lon)
    lon180[lon180>180] = lon180[lon180>180]-360
    i_lon = np.argsort(lon180)
    lon180 = lon180[i_lon]
    return lon180,i_lon