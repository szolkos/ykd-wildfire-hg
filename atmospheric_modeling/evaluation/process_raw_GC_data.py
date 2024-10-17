# --------------
# Name:    process_raw_GC_data.py
# Purpose: A script to process raw GEOS-Chem output into a single xarray Dataset 
#          containing aggregated deposition diagnostics
# Notes:   This script currently does not process concentrations.
#          Script is designed to be called from the command line
#          using the `call_process_raw_GC_data.py` script.
# --------------

from datetime import date
import os
import pandas as pd
import xarray as xr
import json

# --------------
# DEFINE GLOBAL VARIABLES
# --------------
today     = date.today()
today_str = today.strftime("%d-%m-%Y")

seconds_per_year = 3.154e+7
cm2_per_m2       = 1e2*1e2 
avo_num          = 6.022e23
MW_Hg            = 200.59 # molar mass of mercury (g mol-1)
g_per_molec      = MW_Hg*(1./avo_num) # molar mass of Hg * avogadro's number
ug_per_kg        = 1e9

M_air   = 28.9647 # molar mass of dry air (g mol-1)
rho_air = 1292.2  # dry air density at 1 atm, 0 degrees C (g m-3)

# -- read in dictionary of output variable attributes
with open('./output_variable_collections.json') as f:
    attr_dict = json.load(f)

# --------------
# DEFINE THE FUNCTIONS
# --------------

def load_df(path_root:str, coll_names:list, time_slice=None):
    '''
    Description: 
        A function to load a dataframe of GEOS-Chem diagnostics 
    
    Arguments: 
        path_root  (str)  : a string designating a path to the diagnostics
        coll_names (list) : a list of diagnostic collection names, grouped over time using the *
        time_slice (list, optional) : a list containing two strings representting the start and end slice dates
    
    Returns: 
        An xarray Dataset containing merged output from all collections in coll_names
    
    Example Usage: 
        load_df(path_root='Standard/OutputDir/', coll_names=['GEOSChem.WetLossConv*', 'GEOSChem.WetLossLS*'], time_slice=['20150105','20150105'])
    '''

    # ensure type of coll_names is list, even if given as string
    coll_names = list(coll_names)

    if len(coll_names)==1:
        df = xr.open_mfdataset(path_root+coll_names[0])
        
    elif len(coll_names)>1:
        print(coll_names)
        df = xr.open_mfdataset(path_root+coll_names[0])

        for v in coll_names[1:]:
            print(v)
            tmp_df = xr.open_mfdataset(path_root+v)
            df = xr.merge([df, tmp_df])
            print(f'completed merge of {v}')

    # now slice time if time_slice != None
    if time_slice != None:
        df = df.sel(time=slice(time_slice[0], time_slice[1]))
        print(f'completed timeslice from {time_slice[0]} to {time_slice[1]}')
    return df

def aggregate_species(ds, output_variable_name:str):
    '''
    Purpose: 
        aggregate species into collections of species based on data in `attr_dict`

    Arguments:
        ds: xarray Dataset containing GEOS-Chem output
        output_variable_name: string representing the output variable to aggregate
    
    Returns:
        ds: xarray Dataset containing aggregated species
    '''

    # get species characteristics from geoschem_species_metadata.yml
    # -- use Is_Advected, Is_DryDep, Is_WetDep tags to aggregate species
    assert output_variable_name in list(attr_dict.keys()), f'output_variable_name must be one of {list(attr_dict.keys())}'

    units = attr_dict[output_variable_name]['expected_units']

    # add empty variable to hold aggregated species
    ds[output_variable_name] = xr.zeros_like(ds[list(attr_dict[output_variable_name]['variables'])[0]])

    # loop over all variables in attr_dict[output_variable_name]['variables']
    for v in attr_dict[output_variable_name]['variables']:
        # check that units of v match expected units
        assert ds[output_variable_name].units == units, f'units of {output_variable_name} and {v} do not match'
        # add species to output_variable_name
        ds[output_variable_name] += ds[v]
        # now drop species which was just added

    return ds


def process_all(path_root, time_slice1, time_slice2, processed_file_folder_name, 
                subset_domain=True, subset_domain_lat=[30,80], subset_domain_lon=[-175,-105]):
    # ----------------------------------------------------------------------
    # define time_slice and time_save_str
    # --
    if time_slice1 != time_slice2:
        time_slice = [time_slice1, time_slice2]
        time_save_str = f'{time_slice[0][0:4]}-{time_slice[1][0:4]}'
    # otherwise, set time_slice to None
    else:
        time_slice = None
        # save full date in format of YYYYMMDD (typical of GEOS-Chem output)
        time_save_str = f'{time_slice1[0:8]}'
    
    # ----------------------------------------------------------------------
    # define collections and collection_paths
    # --
    collections = ['WetLossConv', 'WetLossLS', 'DryDep', 'MercuryOcean', 'SpeciesConc', 'StateMet']  # ['DryDep'] # 'MercuryEmis',

    if time_slice != None:
        collection_paths = [f'/OutputDir/GEOSChem.{coll}*' for coll in collections]
    else:
        collection_paths = [f'/OutputDir/GEOSChem.{coll}.{time_slice1[0:8]}*.nc4' for coll in collections]

    # ----------------------------------------------------------------------
    # load output and slice time
    base = load_df(path_root = path_root, coll_names=collection_paths, time_slice=time_slice)

    output_variable_list = list(attr_dict.keys())

    for var_name in output_variable_list:
        base = aggregate_species(base, var_name)

    output_vars = list(attr_dict.keys())
    # add AREA and air-sea exchange variables
    output_vars.append('AREA')
    output_vars.append('FluxHg0fromOceanToAir')
    output_vars.append('FluxHg0fromAirToOcean')
    # subset to only output variables
    base = base[output_vars]

    # convert to standard units
    for var_name in output_variable_list:
        # update wet deposition variables
        if 'Wet_Dep' in var_name:
            assert attr_dict[var_name]['output_attributes']['units'] == '\u03bcg m-2 yr-1', f'units of {var_name} are not as expected'
            base[var_name] = base[var_name].sum('lev')/base['AREA'] # kg/m2/s
            base[var_name] = base[var_name]*ug_per_kg*seconds_per_year # μg/m2/yr
            base[var_name] = base[var_name].assign_attrs(attr_dict[var_name]['output_attributes'])
        
        # update dry deposition variables
        if 'Dry_Dep' in var_name:
            assert attr_dict[var_name]['output_attributes']['units'] == '\u03bcg m-2 yr-1', f'units of {var_name} are not as expected'
            base[var_name] = base[var_name]*cm2_per_m2*seconds_per_year*g_per_molec*1e6 # μg m-2 yr-1
            base[var_name] = base[var_name].assign_attrs(attr_dict[var_name]['output_attributes'])

        # update SpeciesConc variables
        if 'SpeciesConc' in var_name:
            assert attr_dict[var_name]['output_attributes']['units'] == 'ng m-3', f'units of {var_name} are not as expected'
            base[var_name] = base[var_name]*((MW_Hg*1e9)/M_air)*rho_air # mol mol-1 dry air --> ng Hg m-3 dry air @ STP (1 atm, 0 C)
            base[var_name] = base[var_name].assign_attrs(attr_dict[var_name]['output_attributes'])

    # add air-sea exchange variables
    ocean_attr_dict = {'Net_Ocean_Hg0_Uptake':  {'species':'gaseous elemental Hg', 'deposition type':'air-sea exchange', 'units':'μg m-2 yr-1'},
                       'FluxHg0fromOceanToAir': {'species':'gaseous elemental Hg', 'deposition type':'air-sea exchange', 'units':'μg m-2 yr-1'},
                       'FluxHg0fromAirToOcean': {'species':'gaseous elemental Hg', 'deposition type':'air-sea exchange', 'units':'μg m-2 yr-1'}}

    base['Net_Ocean_Hg0_Uptake'] = base['FluxHg0fromOceanToAir']-base['FluxHg0fromAirToOcean'] # kg/s
    for v in ['FluxHg0fromOceanToAir', 'FluxHg0fromAirToOcean', 'Net_Ocean_Hg0_Uptake']:
        base[v] = (base[v]/base['AREA'])*seconds_per_year*ug_per_kg # convert kg s-1 to μg m-2 yr-1
        base[v] = base[v].assign_attrs(ocean_attr_dict[v])

    # sum wet and dry deposition
    total_Hg_Dep_attrs = {'species':'all Hg species', 'deposition type':'wet+dry (no air-sea exchange)', 'units':'μg m-2 yr-1'}
    base['Total_Hg_Dep'] = base['Total_Hg_Wet_Dep']+base['Total_Hg_Dry_Dep']
    base['Total_Hg_Dep'] = base['Total_Hg_Dep'].assign_attrs()
    base['Total_Hg_Dep'] = base['Total_Hg_Dep'].assign_attrs(total_Hg_Dep_attrs)

    # ----------------------------------------------------------------------
    # subset domain
    if subset_domain:
        base = base.sel(lat=slice(subset_domain_lat[0], subset_domain_lat[1]), lon=slice(subset_domain_lon[0], subset_domain_lon[1]))

    # ----------------------------------------------------------------------
    # save output
    encoding_dict = {var: {'zlib': True, 'complevel': 9} for var in base.data_vars}
    base.to_netcdf(f'{path_root}{processed_file_folder_name}GC_aggregated_{time_save_str}.nc', encoding=encoding_dict)

    return base

if __name__ == '__main__':
    import argparse

    """ Example Usage: 
            python call_process_raw_GC_data.py -p /Volumes/Expansion/YKD_Fire_v14_3/Standard/ -f Processed_Files/
    """
    # use argparse to get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path_root', default='/Volumes/Expansion/YKD_Fire_v14_3/Standard/', help='path to the root directory containing model output')
    parser.add_argument('-t1', '--time_slice1', default='20150601', help='timestamp of earliest file to process')
    parser.add_argument('-t2', '--time_slice2', default='20150601', help='timestamp of latest file to process')
    parser.add_argument('-f', '--processed_file_folder_name', default='Processed_Files/', help='name of the folder containing processed files')
    parser.add_argument('-sflag', '--subset_domain_flag', default=True, help='whether to subset the domain')
    parser.add_argument('-slat', '--subset_domain_lat', default=[52,74], help='latitude range to subset')
    parser.add_argument('-slon', '--subset_domain_lon', default=[-170,-138.5], help='longitude range to subset')
    args = parser.parse_args()

    # set variables from command line arguments
    path_root   = args.path_root
    time_slice1 = args.time_slice1
    time_slice2 = args.time_slice2
    processed_file_folder_name = args.processed_file_folder_name
    subset_domain_flag = args.subset_domain_flag
    subset_domain_lat = args.subset_domain_lat
    subset_domain_lon = args.subset_domain_lon

    # pass arguments to `process_all()`
    process_all(path_root, time_slice1, time_slice2, processed_file_folder_name,
                subset_domain_flag, subset_domain_lat, subset_domain_lon)