import numpy as np
import pandas as pd
import xarray as xr

def get_total_dep(ds, dep_var='Total_Hg_Deposition', units='kg/mo'):
    """ Calculate total deposition for each grid cell and return a new dataset. """
    
    days_in_year  = 365.25

    if units in ['kg/mo', 'kg mo-1']:
    # convert from ug/m2/a to ug/m2
        days_in_month = ds.time.dt.days_in_month
        weights = days_in_month/days_in_year

    elif units in ['kg/yr', 'kg yr-1']:
        weights = 1
    elif units in ['kg/day', 'kg d-1']:
        weights = 1/days_in_year
    else:
        raise ValueError('units must be one of kg/mo, kg/yr, or kg/day')
    
    assert 'AREA' in ds, 'AREA must be a variable in the dataset'
    assert ds['AREA'].units == 'm2', 'AREA must be in units of m2'
    
    total_area = ds['AREA'].sum(dim=['lat','lon']).values.item()
    ds = ds[dep_var]*ds['AREA']*1e-9
    ds = ds*(weights)
    ds = ds.sum(dim=['lat','lon'])
    
    return {'time':ds.time, f'{dep_var} [{units}]': ds.values, 'area [m2]': total_area}

def scale_ds(f_TGM, E_total):
    """ Used to scale deposition data for uncertainty analysis."""
    # -- emissions used in model simulations [kg] --
    E_TGM_model = 565.9
    E_TPM_model = 21.4
    # -- scale factors --
    scale_TGM = (E_total*f_TGM)/E_TGM_model
    scale_TPM = (E_total*(1-f_TGM))/E_TPM_model

    return {'scale_TGM': scale_TGM, 'scale_TPM': scale_TPM}

def subset_table(table, match_dict:dict):
    """ Return a subset of the table that matches the keys and values in match_dict. """
    subset = table.copy()
    for key in match_dict:
        subset = subset[subset[key]==match_dict[key]]
    return subset