import numpy as np
import xarray as xr

def load_ds(ds_fn, blank_fn, remove_blank=True, directory_path='./'):
    """ Load a dataset from a netCDF file and subtract the blank if remove_blank is True. """
    ds = xr.open_mfdataset(f'{directory_path}{ds_fn}', compat='override', engine='netcdf4')

    # fill attribute dictionary using dictionary comprehension
    attr_dict = {v: ds[v].attrs for v in ds.data_vars}

    area = ds['AREA'].copy()
    # Remove blank from ds
    assert type(remove_blank) == bool, 'remove_blank must be a boolean'
    if remove_blank:
        print('Removing blank... in load_ds()')
        ds_blank = xr.open_mfdataset(f'{directory_path}{blank_fn}', compat='override', engine='netcdf4')
        ds = ds - ds_blank

    ds['AREA'] = area
    
    # Add attributes back in after removing blank
    for v in ds.data_vars:
        ds[v].attrs = attr_dict[v]
    return ds

# return slice of data for a given time period
def get_time_slice(data, start, end):
    data = data.sel(time=slice(start, end))
    return data

def subset_domain(data, lon_min=-170, lon_max=-138.5, lat_min=52, lat_max=74):
    data = data.sel(lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max))
    return data

def get_total_deposition(data):
    assert 'Total_Hg_Wet_Dep' in data.data_vars, 'Total_Hg_Wet_Dep not in data'
    assert 'Total_Hg_Dry_Dep' in data.data_vars, 'Total_Hg_Dry_Dep not in data'
    assert 'FluxHg0fromAirToOcean' in data.data_vars, 'FluxHg0fromAirToOcean not in data'
    assert data['Total_Hg_Wet_Dep'].units == data['Total_Hg_Dry_Dep'].units, 'wet and dry deposition units do not match'
    assert data['Total_Hg_Wet_Dep'].units == data['FluxHg0fromAirToOcean'].units, 'wet deposition and air-sea flux units do not match'
    
    data['Total_Hg_Deposition'] = data['Total_Hg_Wet_Dep'] + data['Total_Hg_Dry_Dep'] + data['FluxHg0fromAirToOcean']
    data['Total_Hg_Deposition'].attrs = {'units': data['Total_Hg_Wet_Dep'].units, 'description': 'Total deposition of Hg (wet + dry + ocean Hg0 uptake)'}
    return data

def get_total_concentration(data):
    assert 'SpeciesConc_Hg0' in data.data_vars, 'SpeciesConc_Hg0 not in data'
    assert 'SpeciesConc_Hg2' in data.data_vars, 'SpeciesConc_Hg2 not in data'
    assert 'SpeciesConc_HgP' in data.data_vars, 'SpeciesConc_HgP not in data'
    assert data['SpeciesConc_Hg0'].units == data['SpeciesConc_Hg2'].units == data['SpeciesConc_HgP'].units, 'species concentration units do not match'

    data['SpeciesConc_TGM'] = data['SpeciesConc_Hg0'] + data['SpeciesConc_Hg2'] + data['SpeciesConc_HgP']
    data['SpeciesConc_TGM'].attrs = {'units': data['SpeciesConc_Hg0'].units, 'description': 'Total gaseous mercury concentration (GEM + GOM + PBM)'}
    return data

def prepare_deposition(directory_path, ds_fn, blank_fn, remove_blank, start, end, lon_min, lon_max, lat_min, lat_max, save_fn=None, fill_negative=False, mask_tiny=True):
    
    print(f'Processing {directory_path}{ds_fn}... remove_blank={remove_blank}... mask_tiny={mask_tiny}... fill_negative={fill_negative}')
    # Load data
    data = load_ds(ds_fn, blank_fn, remove_blank, directory_path)
    # Slice data to only include the time period of interest
    data = get_time_slice(data, start, end)
    # Subset data to only include the domain of interest
    data = subset_domain(data, lon_min, lon_max, lat_min, lat_max)
    # add deposition
    data = get_total_deposition(data)
    # check that deposition not below tolerance
    minimum_allowed = -0.1 # ug/m2/a
    for v in data.data_vars:
        if ('Dep' in v) or ('Uptake' in v):
            if fill_negative:
                assert (data[v] > minimum_allowed).all(), f'{v} below tolerance'
                data[v] = data[v].where(data[v] > 0, 0)
            if mask_tiny:
                data[v] = data[v].where(np.abs(data[v]) > 1e-5, 0)

    output_vars = ['Total_Hg_Deposition', 'Total_Hg_Wet_Dep', 'Total_Hg_Dry_Dep', 'Net_Ocean_Hg0_Uptake', 'FluxHg0fromOceanToAir', 'FluxHg0fromAirToOcean', 'AREA']
    data = data[output_vars]

    area_attrs = data['AREA'].attrs
    data['AREA'] = data['AREA'].mean(dim='time')
    data['AREA'].attrs = area_attrs

    if save_fn is not None:
        print(f'Saving to {save_fn}...')
        data.to_netcdf(save_fn)
    return data

def prepare_species_conc(directory_path, ds_fn, blank_fn, remove_blank, start, end, lon_min, lon_max, lat_min, lat_max, save_fn=None, fill_negative=False, mask_tiny=True):
    
    print(f'Processing {directory_path}{ds_fn}... remove_blank={remove_blank}... mask_tiny={mask_tiny}... fill_negative={fill_negative}')
    # Load data
    data = load_ds(ds_fn, blank_fn, remove_blank, directory_path)
    # Slice data to only include the time period of interest
    data = get_time_slice(data, start, end)
    # Subset data to only include the domain of interest
    data = subset_domain(data, lon_min, lon_max, lat_min, lat_max)
    # add concentration
    data = get_total_concentration(data)
    # check that species concentration not below tolerance
    minimum_allowed = -0.1 # ng m-3
    for v in data.data_vars:
        if 'SpeciesConc' in v:
            if fill_negative:
                assert (data[v] > minimum_allowed).all(), f'{v} below tolerance'
                data[v] = data[v].where(data[v] > 0, 0)
            if mask_tiny:
                data[v] = data[v].where(np.abs(data[v]) > 1e-5, 0)

    output_vars = ['SpeciesConc_TGM', 'SpeciesConc_Hg0', 'SpeciesConc_Hg2', 'SpeciesConc_HgP']
    data = data[output_vars]

    if save_fn is not None:
        print(f'Saving to {save_fn}...')
        data.to_netcdf(save_fn)
    return data

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Prepare data for publication')
    parser.add_argument('--directory_path', type=str, default='./', help='Path to data files')
    parser.add_argument('--ds_fn', type=str, help='Data file to process')
    parser.add_argument('--blank_fn', type=str, help='Blank file to process')
    parser.add_argument('--start', type=str, default='20150501', help='Start date for data')
    parser.add_argument('--end', type=str, default='20160501', help='End date for data')
    parser.add_argument('--lon_min', type=float, default=-170, help='Minimum longitude in domain')
    parser.add_argument('--lon_max', type=float, default=-138.5, help='Maximum longitude in domain')
    parser.add_argument('--lat_min', type=float, default=52, help='Minimum latitude in domain')
    parser.add_argument('--lat_max', type=float, default=74, help='Maximum latitude in domain')
    parser.add_argument('--remove_blank', help='Remove blank from data')
    parser.add_argument('--save_fn_deposition', type=str, default=None, help='Filename to save processed deposition data')
    parser.add_argument('--save_fn_species_conc', type=str, default=None, help='Filename to save processed species concentration data')
    parser.add_argument('--fill_negative', default="False", help='Fill negative values with 0')
    parser.add_argument('--mask_tiny', default="True", help='Mask values below 1e-5')
    args = parser.parse_args()

    if args.remove_blank in ["True", "TRUE", "true", True]:
        remove_blank = True
    elif args.remove_blank in ["False", "FALSE", "false", False]:
        remove_blank = False
    else:
        raise ValueError("remove_blank must be 'True' or 'False' ")
        
    if args.fill_negative in ["True", "TRUE", "true", True]:
        fill_negative = True
    elif args.fill_negative in ["False", "FALSE", "false", False]:
        fill_negative = False
    else:
        raise ValueError("fill_negative must be 'True' or 'False' ")
    
    if args.mask_tiny in ["True", "TRUE", "true", True]:
        mask_tiny = True
    elif args.mask_tiny in ["False", "FALSE", "false", False]:
        mask_tiny = False
    else:
        raise ValueError("mask_tiny must be 'True' or 'False' ")
    prepare_deposition(args.directory_path, args.ds_fn, args.blank_fn, remove_blank, args.start, args.end, args.lon_min, args.lon_max, args.lat_min, args.lat_max, args.save_fn_deposition, fill_negative, mask_tiny)
    prepare_species_conc(args.directory_path, args.ds_fn, args.blank_fn, remove_blank, args.start, args.end, args.lon_min, args.lon_max, args.lat_min, args.lat_max, args.save_fn_species_conc, fill_negative, mask_tiny)