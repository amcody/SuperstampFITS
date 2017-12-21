from pyraf import iraf
from iraf import kepler
import kepprf_AMC
import numpy as np
import os
import glob
from shutil import copyfile
from astropy.wcs.utils import skycoord_to_pixel
from astropy import coordinates
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits

def main():

 files = np.loadtxt("filenames_M8_wcs.txt",dtype='str') 

 for infile in files:
   print infile
   ccmapfile = open('Lagoon_2MASS_xy_radec.dat', 'w')
   newfile = infile.replace('wcs','wcs2')
   copyfile(infile,newfile)

   coords = np.loadtxt('Lagoon_2MASS_radec_amended.dat',dtype=str)
   ra = coords[:,0]
   ra = np.array([float(i) for i in ra])
   dec = coords[:,1]
   dec = np.array([float(i) for i in dec])

#   open image:
   instr = fits.open(infile,mode='readonly',memmap=True)
   pixels = instr[0].data[:]
   crval1p = instr[0].header['CRVAL1P']
   crval2p = instr[0].header['CRVAL2P']   
   instr.close()

   for i in np.arange(len(ra)):
     skyposition = SkyCoord(ra[i], dec[i], unit=('deg','deg'), frame='icrs')
     wcs = WCS(infile)
     pixelpos = skycoord_to_pixel(skyposition, wcs=wcs)
     columns = str(int(round(pixelpos[0])))
     rows = str(int(round(pixelpos[1])))

## Guess at flux:
     fluxguess = str(pixels[int(rows),int(columns)])

     result = kepprf_AMC.kepprf_AMC(infile,'1',columns,rows,fluxguess,border=1,background='yes',focus='no', \
       prfdir='/Users/acody/Data/Kepler',xtol=0.0001,ftol=0.01,verbose=False,logfile='kepprf.log')

     newx = result[1]-crval1p+1
     newy = result[2]-crval2p+1

     print>>ccmapfile, newx, newy, ra[i], dec[i]

   ccmapfile.close()


# Call ccmap to compute new wcs
   iraf.ccmap("Lagoon_2MASS_xy_radec.dat","Lagoon_coordfit.db",images=newfile,lngunit="degrees",latunit="degrees",update="yes",verbose="no",interactive="no")

   os.remove('Lagoon_2MASS_xy_radec.dat')
   os.remove('Lagoon_coordfit.db')


# Standard boilerplate to call the main() function.
if __name__ == '__main__':
  main()
