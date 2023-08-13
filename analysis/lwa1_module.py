import copy
import os
import re
from cmath import nan
from datetime import datetime
from pathlib import Path

import healpy as hp
import numpy as np

now = datetime.now()

data_path = Path("../data/")
LWA1_source_path = Path(data_path/"LWA1_source/")
LWA1_source_path.mkdir(parents=True, exist_ok=True)
commander_in_path = Path(data_path/"commander_in/")
commander_in_path.mkdir(parents=True, exist_ok=True)



def change_coord(m, coord):
    """Change coordinates of a HEALPIX map

    Parameters
    ----------
    m : map or array of maps
      map(s) to be rotated
    coord : sequence of two character
      First character is the coordinate system of m, second character
      is the coordinate system of the output map. As in HEALPIX, allowed
      coordinate systems are 'G' (galactic), 'E' (ecliptic) or 'C' (equatorial)

    Example
    -------
    The following rotate m from galactic to equatorial coordinates.
    Notice that m can contain both temperature and polarization.
    >>>> change_coord(m, ['G', 'C'])
    """
    # Basic HEALPix parameters
    npix = m.shape[-1]
    nside = hp.npix2nside(npix)
    ang = hp.pix2ang(nside, np.arange(npix))

    # Select the coordinate transformation
    rot = hp.Rotator(coord=reversed(coord))

    # Convert the coordinates
    new_ang = rot(*ang)
    new_pix = hp.ang2pix(nside, *new_ang)

    return m[..., new_pix]

freqs_lwa1 = [35, 38, 40, 45, 50, 60, 70, 74, 80]
def create_dict_fits(freqs=freqs_lwa1):
    """Create a dictionary of fits files for the LWA1 maps.

    Parameters
    ----------
    freqs : list
        List of frequencies to create maps for.
        
    Returns
    -------
    dict_maps : dict
        Dictionary of intensity maps for each frequency.
    dict_errs : dict
        Dictionary of error/noise maps for each frequency.
        
    Notes
    -----
    The maps are saved as fits files in the 'commander_inputs' directory.
    The maps are also returned as a dictionary.
    The maps are in galactic/equatorial coordinates.
    The maps are masked for bad pixels.

    """
    dict_maps = {}
    dict_errs = {}
    for f in freqs:
        freq = f"{f}"
        print(f"Creating maps for {freq}mhz...")
        # Import healpix maps in equatorial coords
        file_map = f"healpix-all-sky-rav-wsclean-map-{f}.fits"
        file_err = f"healpix-all-sky-rav-wsclean-err-{f}.fits"
        map_equat = hp.fitsfunc.read_map(LWA1_source_path/file_map, dtype=None)
        err_equat = hp.fitsfunc.read_map(LWA1_source_path/file_err, dtype=None)

        #
        # Change coordinates from equatorial to galactic/celestial
        #
        map_celest = change_coord(map_equat, ['C', 'G'])
        err_celest = change_coord(err_equat, ['C', 'G'])

        # print(f"{f}mhz -- Map size: ", maps_lwa1[f"{f}mhz"].size)
    
        # Mask the bad pixels
        dict_maps[freq] = mask_badsig(map_celest)
        dict_errs[freq] = mask_badsig(err_celest)

    # Create and save mask-maps for both the intensity and error/noise maps
    dict_maps, dict_errs, _, _ = create_masks(dict_maps, dict_errs)

    # #
    # # Save the maps as fits files with the right names etc.
    # #
    # map_fits_file = Path(f"LWA1-00{freq}/LWA1-00{freq}_map.fits")
    # map_fits_file.mkdir(parents=True, exist_ok=True)
    # err_fits_file = Path(f"LWA1-00{freq}/LWA1-00{freq}_err.fits")
    # err_fits_file.mkdir(parents=True, exist_ok=True)
    # hp.fitsfunc.write_map(commander_in_path/map_fits_file, maps[freq], overwrite=True, dtype=None)
    # hp.fitsfunc.write_map(commander_in_path/err_fits_file, errs[freq], overwrite=True, dtype=None)
    
    return dict_maps, dict_errs

def mask_badsig(map):
    """Mask bad pixels in the LWA1 maps.

    Parameters
    ----------
    map : array
        Map to be masked.

    Returns
    -------
    map : array
        Masked map.

    Notes
    -----
    The bad pixels are set to 'nan' values.
    Doing these checks for bad pixels, aka pixel values <= 0. 

    """
    # print(f"\n{freq} -- BEFORE NAN-ing bad-pixels:")
    #
    # Find the nan-valued pixels
    #
    idx_mask = np.argwhere(np.isnan(map))
    # print(f"{freq} -- No. of nan pixels (to be masked): ", idx_mask.size)
    # print("indices for nan pixels (to be masked):\n", idx_mask)

    #
    # Find the <=0 pixels:
    #
    idx_badpix = np.where(map <= 0) 
    #
    # The above is a tuple with at least 1 element, i.e. an empty array.
    # So to check the size of the array contained by the tuple, I need to extract the array (as below)
    #
    idx_badpix = idx_badpix[0]
    # print(f"(1) No. of pixels <= 0: ", idx_badpix.size)
    # print(f"(2) No. of pixels <= 0: ", (map <= 0).sum()) # This is another way to check for pixels valued <=0
    # print("indices of pixels <= 0:\n", idx_badpix)

    #
    # if there are pixels <=0 (there are quite a few for the noise maps, for some reason), we set them to 'nan' as: 
    #
    # print(map[idx_badpix])
    map[idx_badpix] = nan 

    #
    # Now check again to see if there are any bad pixels left...
    #
    idx_badpix = np.where(map <= 1)[0]
    # print("AFTER NAN-ing bad pixels:")
    # print("(1) No. of pixels <= 0:", idx_badpix.size)
    # print("(2) No. of pixels <= 0:", (map <= 0).sum()) # This is another way to check for pixels valued <=0
    # print("indices of pixels <= 0:\n", idx_badpix)

    # 
    # Now we update our array of indices for nan-pixels
    # 
    idx_mask = np.argwhere(np.isnan(map))
    # print("AFTER: No. of nan pixels (to be masked):", idx_mask.size)

    # print(map[idx_mask])
    
    return map

def create_masks(maps, errs):
    """Create masks for the maps and errors.
    
    Args:
        maps (dict): Dictionary of maps
        errs (dict): Dictionary of errors

    Returns:
        maps (dict): Dictionary of maps with masks
        errs (dict): Dictionary of errors with masks
        mask_maps (dict): Dictionary of masks for the maps
        mask_errs (dict): Dictionary of masks for the errors

    Notes:
        The masks are created by checking if the pixel values in the maps and errors are the same.
        If they are not, then the pixel is masked.
        The masks are saved as fits files.

    TODO:
        - Check if the masks are the same for the maps and errors

    """
    for freq,map,err in zip(maps.keys(), maps.values(), errs.values()):
        # This for loop will iterate through the maps and errors, and create masks for each map and error.
        print(f"\nCreating mask-maps for {freq}mhz...")
        
        # Find the nan-valued pixels in the maps
        idx_mask_map = np.argwhere(np.isnan(map))
        idx_mask_err = np.argwhere(np.isnan(err))

        # print("Intensity vs. noise mask-map sizes: ", idx_mask_map.size, idx_mask_err.size)

        # Check if there are any mask-pix in the intensity map that are not in the noise map
        mapdifferr = np.setdiff1d(idx_mask_map, idx_mask_err)
        
        # print("Mask-pix in intensity NOT in noise: ", mapdifferr.size)

        # Check if there are any mask-pix in the noise map that are not in the intensity map
        errdiffmap = np.setdiff1d(idx_mask_err, idx_mask_map)
        
        # print("Mask-pix in noise NOT in intensity: ", errdiffmap.size)

        idx_mask_combined = np.union1d(idx_mask_map, idx_mask_err)
        
        # print("Size of combined pix-mask-map: ", idx_mask_combined.size)
        
        map[idx_mask_combined] = nan
        err[idx_mask_combined] = nan
        maps[freq] = map
        errs[freq] = err

        #
        # Save the maps as fits files with the right names etc.
        #
        map_fits_dir = Path(commander_in_path/f"LWA1-00{freq}/")
        err_fits_dir = Path(commander_in_path/f"LWA1-00{freq}/")
        map_fits_dir.mkdir(parents=True, exist_ok=True)
        err_fits_dir.mkdir(parents=True, exist_ok=True)
        map_fits_file = f"LWA1-00{freq}_map.fits"
        err_fits_file = f"LWA1-00{freq}_err.fits"
        hp.fitsfunc.write_map(map_fits_dir/map_fits_file, maps[freq], overwrite=True, dtype=None)
        hp.fitsfunc.write_map(err_fits_dir/err_fits_file, errs[freq], overwrite=True, dtype=None)

        # if errdiffmap.size > mapdifferr.size:
        #     #
        #     # The mask we'll use is the noise map mask:
        #     #
        #     maps[freq][errdiffmap] = nan
        #     idx_mask = idx_mask_err
        #     #
        #     # Save the map as a fits file with the right name etc.
        #     #
        #     map_fits_file = f"LWA1-00{freq}/LWA1-00{freq}_map.fits"
        #     map_fits_file.mkdir(parents=True, exist_ok=True)
        #     hp.fitsfunc.write_map(commander_in_path/map_fits_file, maps[freq], overwrite=True, dtype=None)
        # else:
        #     #
        #     # The mask we'll use is the intensity map mask:
        #     #
        #     errs[freq][mapdifferr] = nan
        #     idx_mask = idx_mask_map
        #     #
        #     # Save the map as a fits file with the right name etc.
        #     #
        #     err_fits_file = f"LWA1-00{freq}/LWA1-00{freq}_err.fits"
        #     err_fits_file.mkdir(parents=True, exist_ok=True)
        #     hp.fitsfunc.write_map(commander_in_path/err_fits_file, errs[freq], overwrite=True, dtype=None)

        goodmask_map = hp.mask_good(map)
        goodmask_err = hp.mask_good(err)
        idx_mask_map = np.where(goodmask_err==False)[0]
        idx_mask_err = np.where(goodmask_map==False)[0]
        
        # print("Intensity vs. noise mask-map sizes: ", idx_mask_map.size, idx_mask_err.size)

        mask_map = goodmask_map.astype(int)
        mask_err = goodmask_err.astype(int)

        mask_fits_dir = Path(commander_in_path/f"LWA1-00{freq}/")
        mask_fits_dir.mkdir(parents=True, exist_ok=True)
        mask_fits_file = f"LWA1-00{freq}_maskmap.fits"
        hp.fitsfunc.write_map(mask_fits_dir/mask_fits_file, mask_map, overwrite=True, dtype=None)

    return maps, errs, mask_map, mask_err


vals_kK = [[0,40], [0,30], [0,25], [0,20], [0,15], [0,10], [0,7], [0,6], [0,5]]
def init_dict_cbrange(dict_fits, vals_kK=vals_kK):
    """Initialize the dictionaries of colorbar ranges for the LWA1 maps.

    Args:
        dict_fits (dict): Dictionary of fits files for the LWA1 maps.
        vals_kK (list): List of lists of colorbar ranges for the LWA1 maps.

    Returns:
        cbrange_K (dict): Dictionary of colorbar ranges for the LWA1 maps in K.
        cbrange_kK (dict): Dictionary of colorbar ranges for the LWA1 maps in kK.

    Notes:
        The colorbar ranges are in Kelvin (K) and kilo-Kelvin (kK).

    """
    ### Create dictionaries of temperature ranges for healpy colorbars for each freq map ###
    # cb_ranges = dict.fromkeys(maps_lwa1.keys())
    # temp_ranges = {x:maps_lwa1[freq] for freq in maps_lwa1.keys()}
    vals_kK = [[float(y) for y in x] for x in vals_kK]
    vals_K = [[y*1000.0 for y in x] for x in vals_kK]

    cbrange_kK=dict(zip(dict_fits.keys(), vals_kK))
    cbrange_K=dict(zip(dict_fits.keys(), vals_K))

    # print(cbrange_K)
    
    return cbrange_K, cbrange_kK


def pct_errs(maps, errs):
    """Calculate the percentage errors in the maps.

    Args:
        maps (dict): Dictionary of maps
        errs (dict): Dictionary of errors

    Returns:
        pct_errs (dict): Dictionary of percentage errors

    Notes:
        The percentage errors are calculated as the errors divided by the maps, multiplied by 100.

    """
    pct_errs = dict.fromkeys(maps.keys())
    for freq,map,err in zip(maps.keys(), maps.values(), errs.values()):
        pct_errs[freq] = err/map * 100.0

    return pct_errs



#
# To find FWHM, we're given the beam size/FWHM in degrees, so just convert to radians (beamsize*pi/180)
# To calculate l_{max}: nside*3 - 1 
#   (but add a small value to the l_max to avoid aliasing??)
# Compute b_{lm} using healpy.spherfunc.gauss_beam with FWHM and l_{max} as input and beam window function [0, l_{max}].
# Use Healpy to write the calculated C(l) to a fits file for each frequency.
#
# Beam sizes:
# 35 MHz - Stokes I - 4.7° beam
# 38 MHz - Stokes I - 4.3° beam
# 40 MHz - Stokes I - 4.1° beam
# 45 MHz - Stokes I - 3.6° beam
# 50 MHz - Stokes I - 3.3° beam
# 60 MHz - Stokes I - 2.7° beam
# 70 MHz - Stokes I - 2.3° beam
# 74 MHz - Stokes I - 2.2° beam
# 80 MHz - Stokes I - 2.0° beam
#
# Create a dictionary of beam sizes for each frequency (using the data in the above comment, with the frequencies as string type keys and the beam sizes as values):
beam_sizes = {'35': 4.7, '38': 4.3, '40': 4.1, '45': 3.6, '50': 3.3, '60': 2.7, '70': 2.3, '74': 2.2, '80': 2.0}
# Create a function to calculate the beam window function for each frequency (without using a for loop, if possible), write the beam window function to a fits file, and return the beamcls as a dictionary (with the frequencies as keys and the beam window functions as values) and the lmax as a separate dictionary (with the frequencies as keys and the lmax as values):
def calc_beam_window_funcs(beam_sizes=beam_sizes):
    """Calculate the beam window function for each frequency and write the beam window function to a fits file.
    
    Args:
        beam_sizes (dict): Dictionary of beam sizes for each frequency. Defaults to beam_sizes in the global namespace.
        
    Returns:
        beamcls (dict): Dictionary of beam window functions for each frequency.
        lmax (dict): Dictionary of lmax for each frequency.

    Notes:
        The beam window function is calculated using the Healpy function hp.spherfunc.gauss_beam() and written to a fits file using the Healpy function hp.fitsfunc.write_cl().
        The beam window function is returned as a dictionary with the frequencies as keys and the beam window functions as values.
        The lmax is returned as a separate dictionary with the frequencies as keys and the lmax as values.

    """
    beamcls = {}
    lmax_dict = {}
    for freq,beamsize in beam_sizes.items():
        fwhm = beamsize*np.pi/180
        nside = 256
        lmax = (nside * 3 - 1) + 10 # Add a small value to the l_max to avoid aliasing?? Or is there another reason?

        beamcls[freq] = hp.sphtfunc.gauss_beam(fwhm=fwhm, lmax=lmax)
        lmax_dict[freq] = lmax

        beam_fits_dir = Path(commander_in_path/f"LWA1-00{freq}/")
        beam_fits_dir.mkdir(parents=True, exist_ok=True)
        beam_fits_file = f"LWA1-00{freq}_beam.fits"
        hp.fitsfunc.write_cl(filename=beam_fits_dir/beam_fits_file, cl=beamcls[freq], overwrite=True, dtype=None)
    
    return beamcls, lmax_dict


# Create a function to calculate the monopole (mp) and dipole (dp) using the maps_lwa1 dictionary at every frequency. Use the Healpy function hp.fit_dipole() to calculate the monopole and dipole values (collectively abbreviated as md). The function should return the md values as a dictionary with the frequencies as keys and the md values as values. Also, print the results, with Monopole and Dipole given separately.
def calc_md(maps):
    """Calculate the monopole and dipole for each frequency.
    
    Args:
        maps_lwa1 (dict): Dictionary of maps for each frequency.
        
    Returns:
        md (dict): Dictionary of monopole and dipole values for each frequency.

    Notes:
        The monopole and dipole values are calculated using the Healpy function hp.fit_dipole() and returned as a dictionary with the frequencies as keys and the monopole and dipole values as values.
        The monopole and dipole values are printed separately.

    """
    md = {}
    for freq,map in maps.items():
        md[freq] = hp.fit_dipole(map)
        print(f"{freq} MHz Monopole: {md[freq][0]}")
        print(f"{freq} MHz Dipole: {md[freq][1]}")
    
    return md



# Extract monopole/dipole from the map
# gcut = 0.0
# gcut = 15.0
# gcut = 30.0
# map_80 = copy.deepcopy(maps_lwa1['80'])
# md_80mlwa1 = hp.fit_dipole(map_80,  gal_cut=gcut)
# mp_80mlwa1 = hp.fit_monopole(map_80, gal_cut=gcut)

# print(md_80mlwa1)
# print(mp_80mlwa1)

# Write a function to perform gal_cut on maps for all nine LWA1 frequency bands (with the map data input as a dictionary), apply the gal_cut, find the monopole and dipole using hp.fit_dipole, and return the monopole and dipole values as a dictionary with the frequencies as keys and the monopole and dipole values as values. Also, print the results, with Monopole and Dipole given separately.
def gal_cut_md(maps, gcut):
    """Perform gal_cut on maps for all nine LWA1 frequency bands and calculate the monopole and dipole for each frequency.
    
    Args:
        maps (dict): Dictionary of maps for each frequency.
        gcut (float): The gal_cut value to use.
        
    Returns:
        md (dict): Dictionary of monopole and dipole values for each frequency.

    Notes:
        The monopole and dipole values are calculated using the Healpy function hp.fit_dipole() and returned as a dictionary with the frequencies as keys and the monopole and dipole values as values.
        The monopole and dipole values are printed separately.

    """
    md = {} # Initialize the dictionary to hold the monopole and dipole values.
    for freq,map in maps.items():
        map_gcut = copy.deepcopy(map) # Create a copy of the map to apply the gal_cut to.
        thetadiff = np.radians(gcut) # Convert the gal_cut value to radians.
        theta1 = np.pi/2 - thetadiff 
        theta2 = np.pi/2 + thetadiff
        ipix_strip = hp.query_strip(nside=256, theta1=theta1, theta2=theta2, inclusive=True) # Find the pixels in the strip.
        map_gcut[ipix_strip] = nan # Set the pixels in the strip to nan.
        md[freq] = hp.fit_dipole(map_gcut) # Calculate the monopole and dipole for the map with the gal_cut applied.
        print(f"{freq} MHz Monopole: {md[freq][0]}")
        print(f"{freq} MHz Dipole: {md[freq][1]}")

    return md

# # Calculate the monopole and dipole for each frequency using the gal_cut_md() function.
# gcut = 0.0
# gcut = 30.0
# gcut = 15.0
# md_gcut = gal_cut_md(maps_lwa1, gcut)





# Create a dictionary using the LWA1 survey parameters above with the frequency as the string type key and the centre frequency as the value, e.g. '35': 34.979
lwa1_centre_freqs = {'35': 34.979, '38': 38.042, '40': 40.052, '45': 44.933, '50': 50.005, '60': 59.985, '70': 70.007, '74': 73.931, '80': 79.960}
lwa1_defaults_dir = '/mn/stornext/d16/cmbco/AST9240/2022/jibran/commander_inputs/defaults/'
lwa1_0080_defaults_filename = 'LWA1-0080.defaults'
lwa1_0080_defaults_file = os.path.join(lwa1_defaults_dir, lwa1_0080_defaults_filename)      # join the path to the LWA1-0080.defaults file with the filename to get the full path to the LWA1-0080.defaults file
beamcls, lmax_dict = calc_beam_window_funcs(beam_sizes) # Calculate the beam window functions for the LWA1 freq bands
def create_lwa1_defaults_files(lwa1_centre_freqs=lwa1_centre_freqs, lwa1_defaults_template_file=lwa1_0080_defaults_file, template_freq=80, lwa1_defaults_dir=lwa1_defaults_dir):
    """Create .defaults files for the LWA1 freq bands.

    This function creates .defaults files for the LWA1 freq bands by modifying a template 
    .defaults file for a given freq band and returns a list of the .defaults files for the
    LWA1 freq bands.

    Parameters
    ----------
    lwa1_centre_freqs : dict
        A dictionary containing the centre frequencies of the LWA1 freq bands.
    lwa1_defaults_template : str
        The .defaults template file for an LWA1 freq band.
    template_freq : int
        The freq band of the .defaults template file.

    Returns
    -------
    lwa1_defaults_files : list
        A list of the .defaults files for the LWA1 freq bands.
    """
    with open(lwa1_defaults_template_file, 'r') as f: # open the LWA1-0080.defaults file in read mode
        lwa1_defaults = f.readlines() # read the lines in the LWA1-0080.defaults file and store them in a list called lwa1_80mhz_defaults
        # Ignore (but do not remove!) commented lines/newline characters/whitespaces in the defaults file and read the rest of the lines in the .defaults file for the purpose of reading the data in the .defaults file.
        # The reason for not removing the commented lines/newline characters/whitespaces is to remember the structure in order to put it back in the output .defaults file such that the output .defaults file for each freq band is identical to the input .defaults file for LWA1-0080 except for the parameters that need to be changed.
        # Ignore commented lines and read the rest of the lines in the .defaults file
        # lwa1_80mhz_defaults = [line for line in lwa1_80mhz_defaults if not line.startswith('#')] # ignore commented lines in the .defaults file
        # lwa1_80mhz_defaults = [line for line in lwa1_80mhz_defaults if line.strip()] # ignore newline characters in the .defaults file
        # lwa1_80mhz_defaults = [line for line in lwa1_80mhz_defaults if not line.isspace()] # ignore whitespaces in the .defaults file
        # print(lwa1_80mhz_defaults)

        for freq in lwa1_centre_freqs: # loop over the frequencies in the lwa1_centre_freqs dictionary (keys)
            # Create a new .defaults file for each freq band and save it in the same directory as the LWA1-0080.defaults file with the filename LWA1-00<freq>.defaults, i.e. prepend '00' to the freq band and append '.defaults' to the filename, e.g. LWA1-0035.defaults
            print("\n\nFrequency band: " + freq + " MHz")

            new_filename = 'LWA1-00' + freq + '.defaults' # create the new filename for each freq band
            new_file = os.path.join(lwa1_defaults_dir, new_filename) # join the path to the LWA1-0080.defaults file with the new filename to get the full path to the new .defaults file

            # If the new .defaults file already exists, delete it and create a new one
            if os.path.exists(new_file): # if the new .defaults file already exists
                os.remove(new_file) # delete the new .defaults file
                print("Deleted existing .defaults file for LWA1-00" + freq + " MHz...")
                print("Creating new .defaults file for LWA1-00" + freq + " MHz...")
            else: # if the new .defaults file does not exist
                print("Creating new .defaults file for LWA1-00" + freq + " MHz...")
            open(new_file, 'w+').close() # create a new .defaults file for each freq band. The 'w+' argument in the open() function creates the file if it does not exist and opens it in write mode. The close() function closes the file.

            with open(new_file, 'r+') as f: # open the new .defaults file in write mode and create the file if it does not exist
                # # It seems that the file is not being created. Check if the file exists and if it does not exist, print 'File does not exist' and if it does exist, print 'File exists'. Use os module to check if the file itself exists, not just the path to the file.
                # if not os.path.exists(new_file): # if the file does not exist
                #     print('File does not exist') # print 'File does not exist'
                # else: # if the file exists
                #     print('File exists') # print 'File exists'
                # # Another test to check if the file exists:
                # try:
                #     open(new_file, 'r')
                #     print('File exists')
                # except FileNotFoundError:
                #     print('File does not exist')

                # Write the modified lines to the new .defaults file
                print("Writing modified lines to the new .defaults file...")
                for line in lwa1_defaults: # loop over the lines in the lwa1_80mhz_defaults list
                    # The first line is a comment with the freq band, e.g. # 80 MHz. Replace the freq band in the first line with the freq band in the lwa1_centre_freqs dictionary.
                    if lwa1_defaults.index(line) == 0: # if the line is the first line in the .defaults file
                        print("First line found...")
                        print("Line: " + line)
                        
                        new_line = line.replace(f'{template_freq}', freq) # replace the string f'{template_freq}' with the freq band in the line
                        f.write(new_line) # replace the string f'{template_freq}' with the freq band in the line and write the modified line to the new .defaults file.
                        print("New line: ", new_line)

                        # Read and print the modified line from the new .defaults file to check if the line was modified correctly
                        f.seek(0) # move the file pointer to the beginning of the file
                        modified_line = f.readlines()[0] # read the first line in the new .defaults file
                        print("Modified line: ", modified_line)
                        
                    if 'BAND_LABEL&&&' in line: # if the line contains the string 'BAND_LABEL&&&'
                        print("BAND_LABEL&&& found in the line...")
                        print("Line: " + line)
                        line_number = lwa1_defaults.index(line) # save the line number of the line that contains the string 'BAND_LABEL&&&' in the .defaults file

                        new_line = line.replace(f'{template_freq}', freq) # replace the string f'{template_freq}' with the freq band in the line
                        f.write(new_line) # replace the string f'{template_freq}' with the freq band in the line and write the modified line to the new .defaults file. 
                        print("New line: ", new_line)
                        
                        # Read and print the modified line from the new .defaults file to check if the line was modified correctly
                        f.seek(0) # move the file pointer to the beginning of the file
                        modified_line = f.readlines()[line_number] # read the modified line from the new .defaults file
                        print("Modified line: " + modified_line)

                    elif 'BAND_LMAX&&&' in line: # if the line contains the string 'BAND_LMAX&&&'
                        print("BAND_LMAX&&& found in the line...")
                        print("Line: " + line)
                        line_number = lwa1_defaults.index(line) # save the line number of the line that contains the string 'BAND_LABEL&&&' in the .defaults file

                        # Replace the last string (that comes after the equal sign '=' and before the inline-comment sentence of arbitrary length; the comment is marked with a # sign) in the line with the lmax value for the freq band in the lmax_dict dictionary (converted to a string) and add a newline character to the end of the line. Make sure to leave the inline-comment sentence at the end of the line. Replace the whitespace characters before and after the lmax value with the same number of whitespace characters.
                        new_line = re.sub(r'= (\S*)(\s*)(#.*)', r'= ' + str(lmax_dict[freq]) + r'\2\3', line) 
                        
                        # Explanation of the regex: r'(\s*)(\S*)(\s*)(#.*)' - match any number of whitespace characters (\s*), any number of non-whitespace characters (\S*), any number of whitespace characters (\s*), and any number of characters that are not a newline character (#.*) - and replace them with any number of whitespace characters (\1), the lmax value for the freq band (str(lmax_dict[freq])), any number of whitespace characters (\3), and the inline-comment sentence at the end of the line (\4). The \<number> in the regex is a backreference to the group with the same number. The regex is case-sensitive.
                        
                        f.write(new_line) # replace the string f'{template_freq}' with the freq band in the line and write the modified line to the new .defaults file.
                        print("New line: ", new_line)

                        # Read and print the modified line from the new .defaults file to check if the line was modified correctly
                        f.seek(0) # move the file pointer to the beginning of the file
                        modified_line = f.readlines()[line_number] # read the modified line from the new .defaults file
                        print("Modified line: " + modified_line)

                        

                    elif 'BAND_MAPFILE&&&' in line: # if the line contains the string 'BAND_MAPFILE&&&'
                        print("BAND_MAPFILE&&& found in the line...")
                        print("Line: " + line)
                        line_number = lwa1_defaults.index(line) # save the line number of the line that contains the string 'BAND_LABEL&&&' in the .defaults file
                        
                        new_line = line.replace(f'{template_freq}', freq) # replace the string f'{template_freq}' with the freq band in the line
                        f.write(new_line) # replace the string f'{template_freq}' with the freq band in the line and write the modified line to the new .defaults file
                        print("New line: ", new_line)

                        # Read and print the modified line from the new .defaults file to check if the line was modified correctly
                        f.seek(0) # move the file pointer to the beginning of the file
                        modified_line = f.readlines()[line_number] # read the modified line from the new .defaults file
                        print("Modified line: " + modified_line)

                    elif 'BAND_NOISEFILE&&&' in line: # if the line contains the string 'BAND_NOISEFILE&&&'
                        print("BAND_NOISEFILE&&& found in the line...")
                        print("Line: " + line)
                        line_number = lwa1_defaults.index(line) # save the line number of the line that contains the string 'BAND_LABEL&&&' in the .defaults file

                        new_line = line.replace(f'{template_freq}', freq) # replace the string f'{template_freq}' with the freq band in the line
                        f.write(new_line) # replace the string f'{template_freq}' with the freq band in the line and write the modified line to the new .defaults file
                        print("New line: ", new_line)

                        # Read and print the modified line from the new .defaults file to check if the line was modified correctly
                        f.seek(0) # move the file pointer to the beginning of the file
                        modified_line = f.readlines()[line_number] # read the modified line from the new .defaults file
                        print("Modified line: " + modified_line)

                    elif 'BAND_MASKFILE&&&' in line: # if the line contains the string 'BAND_MASKFILE&&&'
                        print("BAND_MASKFILE&&& found in the line...")
                        print("Line: " + line)
                        line_number = lwa1_defaults.index(line) # save the line number of the line that contains the string 'BAND_LABEL&&&' in the .defaults file

                        new_line = line.replace(f'{template_freq}', freq) # replace the string f'{template_freq}' with the freq band in the line
                        f.write(new_line) # replace the string f'{template_freq}' with the freq band in the line and write the modified line to the new .defaults file
                        print("New line: ", new_line)

                        # Read and print the modified line from the new .defaults file to check if the line was modified correctly
                        f.seek(0) # move the file pointer to the beginning of the file
                        modified_line = f.readlines()[line_number] # read the modified line from the new .defaults file
                        print("Modified line: " + modified_line)

                    elif 'BAND_BEAM_B_L_FILE&&&' in line: # if the line contains the string 'BAND_BEAM_B_L_FILE&&&'
                        print("BAND_BEAM_B_L_FILE&&& found in the line...")
                        print("Line: " + line)
                        line_number = lwa1_defaults.index(line) # save the line number of the line that contains the string 'BAND_LABEL&&&' in the .defaults file

                        new_line = line.replace(f'{template_freq}', freq) # replace the string f'{template_freq}' with the freq band in the line
                        f.write(new_line) # replace the string f'{template_freq}' with the freq band in the line and write the modified line to the new .defaults file
                        print("New line: ", new_line)

                        # Read and print the modified line from the new .defaults file to check if the line was modified correctly
                        f.seek(0) # move the file pointer to the beginning of the file
                        modified_line = f.readlines()[line_number] # read the modified line from the new .defaults file
                        print("Modified line: " + modified_line)

                    elif 'BAND_NOMINAL_FREQ&&&' in line: # if the line contains the string 'BAND_NOMINAL_FREQ&&&'
                        print("BAND_NOMINAL_FREQ&&& found in the line...")
                        print("Line: " + line)
                        line_number = lwa1_defaults.index(line) # save the line number of the line that contains the string 'BAND_LABEL&&&' in the .defaults file
                        
                        # Replace the last string (that comes after the equal sign '=' and before the inline-comment sentence of arbitrary length; the comment is marked with a # sign) in the line with the lwa_centre_freqs value for the freq band in the lwa_centre_freqs dictionary (converted to a string) and add a newline character to the end of the line. Make sure to leave the inline-comment sentence at the end of the line. Replace the whitespace characters before and after the lmax value with the same number of whitespace characters.
                        new_line = re.sub(r'= (\S*)(\s*)(#.*)', r'= ' + str(round(lwa1_centre_freqs[freq]/1000, 6)) + r'\2\3', line)
                        # new_line = line.replace('0.079960', str(lwa1_centre_freqs[freq]/1000)) # replace the string f'{template_freq}' with the freq band in the line

                        f.write(new_line) # replace the string '0.079960' with the centre frequency of the freq band in the line and write the modified line to the new .defaults file
                        print("New line: ", new_line)

                        # Read and print the modified line from the new .defaults file to check if the line was modified correctly
                        f.seek(0) # move the file pointer to the beginning of the file
                        modified_line = f.readlines()[line_number] # read the modified line from the new .defaults file
                        print("Modified line: " + modified_line)

                    else: # if the line does not contain any of the strings above
                        f.write(line) # write the line to the new .defaults file as it is

    # Create a list of the .defaults files for the LWA1 freq bands
    lwa1_defaults_files = [lwa1_defaults_template_file.replace(f'{template_freq}', freq) for freq in lwa1_centre_freqs] # replace the string f'{template_freq}' in the filename with the freq band and append the modified filename to the lwa1_defaults_files list
    # print(lwa1_defaults_files)

    return lwa1_defaults_files