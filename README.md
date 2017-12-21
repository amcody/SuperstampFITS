# SuperstampFITS
These scripts can be used to create individual FITS files of K2 superstamp regions from the original
target pixel files. To use this version of the astrometry.net API and determine a WCS solution, you must have 
simplejson in your python path, as well as version 2 of the Kepler PyKEtools (https://keplergo.arc.nasa.gov/PyKE.shtml).

Follow these steps to reproduce the superstamp fits files that we have released at the MAST archive at
archive.stsci.edu/prepds/k2superstamp/

0.) Pick a superstamp for which you want fits files.

1.) Gather all target pixel files associated with the superstamp of interest and place them in a directory.
    Go into the appropriate assemble_fits_C*.py file and change the file path in line 13 ('filenames = ...')
    to reflect that directory path. Run assemble_fits_C*.py.

2.) Make a list of the resulting fits files and put them in a text file called filenames_XX.txt, where XX is your cluster (e.g., M67).

3.) To determine an accurate WCS solution for your files you need an account with astrometry net, and an
    associated API key. (Alternatively you can follow their instructions to do a local installation and
    skip this step). Go into the file fixwcs_api_C*.py and insert your API key in the cli.login() line. Do
    the same in the file myclient.py for the variable apiKey. Run fixwcs_api_C*.py.

4.) Make a list of the results files *wcs.fits in the text file filenames_XX_wcsfinal.txt, where XX is again your cluster.

4a.) If you are interested in a very accurate WCS solution for the Lagoon region (Campaign 9), run the intermediate
    script adjust_wcs_C9-Lagoon.py. You need to put the input fits files in a list filenames_M8_wcs.txt. Run the script and 
    then save a list of the resulting *wcs2.fits files into filenames_M8_wcsfinal.txt.

5.) To do the final header clean-up that was necessary for us to meet the criteria for a High Level Science Product, run
    the cleanheader_c*.py script.

