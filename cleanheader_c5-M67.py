import numpy as np
from astropy.io import fits
from astropy.wcs.utils import pixel_to_skycoord
from astropy.wcs import WCS
import datetime
from astropy.time import Time
import glob
import os

def main():

    prefix='hlsp_k2superstamp_k2_photometer_c0-m67'
    suffix='_kepler_v1_image'
    oldprefix='M67'

    filenames = np.loadtxt('filenames_M67_wcsfinal.txt',dtype='str')

    tpfs = glob.glob('/Volumes/Work/Field_5/M67/TPFs/ktwo200*targ.fits')

    for oldoutfile in filenames:

     time = oldoutfile[7:19]
     outfile = prefix + '-bjd' + time + suffix + '.fits'

     file = fits.open(oldoutfile,mode='readonly',memmap=True)
     cards0 = file[0].header.cards
     bigarr = file[0].data

     hdu0 = fits.PrimaryHDU(bigarr)

     for n in range(len(cards0)):
            try:
                if cards0[n].keyword not in hdu0.header.keys():
                    hdu0.header[cards0[n].keyword] = (cards0[n].value,
                                                      cards0[n].comment)
                else:
                    hdu0.header.cards[cards0[n].keyword].comment = cards0[n].comment
            except:
                pass

# Remove old WCS keywords

     hdu0.header.remove('PC1_1')
     hdu0.header.remove('PC1_2')
     hdu0.header.remove('PC2_1')
     hdu0.header.remove('PC2_2')

     for k in range(len(cards0)):
      if (cards0[k].keyword[0] == '_'): hdu0.header.remove(cards0[k].keyword)

     hdu0.header.remove('CHECKSUM')
     hdu0.header.remove('DATASUM')
     hdu0.header.remove('LATPOLE')
     hdu0.header.remove('LONPOLE')
     hdu0.header.remove('IMAGEW')
     hdu0.header.remove('IMAGEH')
     hdu0.header.remove('COMMENT')
     hdu0.header.remove('HISTORY')

# Get RA, Dec of image center
     wcs = WCS(oldoutfile)

     xmid = 200.5
     ymid = 200.5
     skypos = pixel_to_skycoord(xmid,ymid,wcs=wcs)
     ramid = skypos.ra.deg
     decmid = skypos.dec.deg

     hdu0.header['RA_OBJ'] = ramid
     hdu0.header['DEC_OBJ'] = decmid

     hdu0.header.cards['RA_OBJ'].comment = '[deg] right ascension of field center'
     hdu0.header.cards['DEC_OBJ'].comment = '[deg] declination of field center'

# Add some history of file generation

     hdu0.header.set('HISTORY','This FITS file is a mosaic of single cadence K2 data')
     hdu0.header.set('HISTORY','  extracted from a collection of target pixel files (TPFs).')
     hdu0.header.set('HISTORY','It was generated from the following TPFs:')
     for file in tpfs:
       hdu0.header.set('HISTORY',' '+file[31:])
     hdu0.header.set('HISTORY','Target pixel file stamps stitched together according to the')
     hdu0.header.set('HISTORY','  1CRV4P and 2CRV4P keywords')
     hdu0.header.set('HISTORY','WCS solution provided by the astrometry.net suite.')
     hdu0.header.set('HISTORY','For more details, see http://astrometry.net.')
     hdu0.header.set('HISTORY','The scripts used to execute these steps are available at')
     hdu0.header.set('HISTORY','  https://github.com/amcody/SuperstampFITS')
     hdu0.header.set('HISTORY','This file originated from https://archive.stsci.edu/prepds/k2superstamp/')

     hdu0.writeto(outfile,checksum=True)


if __name__ == '__main__':
  main()

