import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io.fits.card import Undefined, UNDEFINED
from astropy.wcs.utils import pixel_to_skycoord
from astropy.wcs import WCS
import datetime 
from astropy.time import Time
import glob

def main():

    prefix='hlsp_k2superstamp_k2_photometer_c9-lagoon'
    suffix='_kepler_v1_image'
    oldprefix='M8'

    filenames = glob.glob('/Volumes/Work/Field_9/TPFs_part1_superstamp/ktwo200*targ.fits')
#    filenames = glob.glob('/Volumes/Work/Field_9/TPFs_part2_superstamp/ktwo200*targ.fits')

    midfile = fits.open(filenames[0],mode='readonly',memmap=True)
    cards0 = midfile[0].header.cards
    cards1 = midfile[1].header.cards
    cards2 = midfile[2].header.cards

    time = midfile[1].data.field('TIME')[:] + midfile[1].header['BJDREFI'] + midfile[1].header['BJDREFF']
    timecorr = midfile[1].data.field('TIMECORR')[:]
    cadenceno = midfile[1].data.field('CADENCENO')[:]
    quality = midfile[1].data.field('QUALITY')[:]
    cosmic_rays = midfile[1].data.field('COSMIC_RAYS')[:]
    pos_corr1 = midfile[1].data.field('POS_CORR1')[:]
    pos_corr2 = midfile[1].data.field('POS_CORR2')[:]

    bigarr = np.empty([135,220])
    bigarr[:] = np.nan
    bigarr = bigarr.astype('float32')

# For part1, images 311-313 have a header problem...
    for i in range(0,1):
#    for i in range(0,310):
#    for i in range(314,1290):

# For part2, images 503-505 have a header problem...
#    for i in range(0,502):
#    for i in range(502,503):
#    for i in range(506,2022):

# construct output primary extension
       outfile = prefix + '-bjd%.4f' % time[i] + suffix + '.fits'
       oldoutfile = 'Images_wcs2/'+oldprefix + '_BJD%.4f' % time[i] + '_wcs2.fits'
       print outfile
       print oldoutfile

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

# Fix a few things
       hdu0.header['DATE'] = str(datetime.date.today())
       hdu0.header['CREATOR'] = 'Ann Marie Cody'
       hdu0.header.set('OBJECT','M8','target')
#       hdu0.header['KEPLERID'] = UNDEFINED
#       hdu0.header['RA_OBJ'] = UNDEFINED
#       hdu0.header['DEC_OBJ'] = UNDEFINED
       hdu0.header['CHANNEL'] = UNDEFINED
       hdu0.header.cards['CHANNEL'].comment = 'double-valued for this superstamp'
       hdu0.header['OUTPUT'] = UNDEFINED
       hdu0.header.cards['OUTPUT'].comment = 'double-valued for this superstamp'

       hdu0.header.remove('PROCVER')
       hdu0.header.remove('FILEVER')
       hdu0.header.remove('TIMVERSN')

       hdu0.header.remove('KEPLERID')
       hdu0.header.remove('TTABLEID')
       hdu0.header.remove('DATA_REL')
       hdu0.header.remove('PMRA')
       hdu0.header.remove('PMDEC')
       hdu0.header.remove('PMTOTAL')
       hdu0.header.remove('PARALLAX')
       hdu0.header.remove('GLON')
       hdu0.header.remove('GLAT')
       hdu0.header.remove('GMAG')
       hdu0.header.remove('RMAG')
       hdu0.header.remove('IMAG')
       hdu0.header.remove('ZMAG')
       hdu0.header.remove('JMAG')
       hdu0.header.remove('HMAG')
       hdu0.header.remove('KMAG')
       hdu0.header.remove('KEPMAG')
       hdu0.header.remove('GRCOLOR')
       hdu0.header.remove('JKCOLOR')
       hdu0.header.remove('GKCOLOR')
       hdu0.header.remove('TEFF')
       hdu0.header.remove('LOGG')
       hdu0.header.remove('FEH')
       hdu0.header.remove('EBMINUSV')
       hdu0.header.remove('AV')
       hdu0.header.remove('RADIUS')
       hdu0.header.remove('TMINDEX')

# construct output image extension
       for fn in filenames:
           with fits.open(fn) as f:
               x = f[1].header['1CRV4P']
               y = f[1].header['2CRV4P']
               ch = f[0].header['CHANNEL']
               mod = f[0].header['MODULE']
               out = f[0].header['OUTPUT']
               alldata = f[1].data['FLUX']

               dim = alldata.shape
               if (out == 1):
                 bigarr[y-868+10:y-868+10+dim[1],x-1002:x-1002+dim[2]] = alldata[i]
               if (out == 2):
                 bigarr[y-858:y-858+dim[1],1112-(x+dim[2])+110:1112-x+110] = np.fliplr(alldata[i])

       for k in range(len(cards1)):
            if (cards1[k].keyword not in hdu0.header.keys() and
                cards1[k].keyword[:4] not in ['TTYP','TFOR','TUNI','TDIS','TDIM','WCAX','1CTY',
                                          '2CTY','1CRP','2CRP','1CRV','2CRV','1CUN','2CUN',
                                          '1CDE','2CDE','1CTY','2CTY','1CDL','2CDL','11PC',
                                          '12PC','21PC','22PC','WCSN','TFIE']):
                hdu0.header.set(cards1[k].keyword, cards1[k].value, cards1[k].comment)

       try:
           int_time = cards1['INT_TIME'].value
       except:
           print 'WARNING -- KEPIMAGES: cannot find INT_TIME keyword'
       try:
           frametim = cards1['FRAMETIM'].value
       except:
           print 'WARNING -- KEPIMAGES: cannot find FRAMETIM keyword'
       try:
           num_frm = cards1['NUM_FRM'].value
       except:
           print 'WARNING -- KEPIMAGES: cannot find NUM_FRM keyword'
       try:
           hdu0.header.set('TELAPSE',frametim * num_frm,'[s] elapsed time for exposure')
       except:
           hdu0.header.set('TELAPSE',-999,'[s] elapsed time for exposure')
       try:
           hdu0.header.set('LIVETIME',int_time * num_frm,'[s] TELAPSE multiplied by DEADC')
       except:
           hdu0.header.set('LIVETIME',-999,'[s] TELAPSE multiplied by DEADC')
       try:
           hdu0.header.set('EXPOSURE',int_time * num_frm,'[s] time on source')
       except:
           hdu0.header.set('EXPOSURE',-999,'[s] time on source')
       try:
           hdu0.header.set('MIDTIME',time[i],'[BJD] mid-time of exposure')
       except:
           hdu0.header.set('MIDTIME',-999,'[BJD] mid-time of exposure')
       try:
           hdu0.header.set('TIMECORR',timecorr[i],'[d] barycenter - timeslice correction')
       except:
           hdu0.header.set('TIMECORR',-999,'[d] barycenter - timeslice correction')
       try:
           hdu0.header.set('CADENCEN',cadenceno[i],'unique cadence number')
       except:
           hdu0.header.set('CADENCEN',-999,'unique cadence number')
       try:
           hdu0.header.set('QUALITY',quality[i],'pixel quality flag')
       except:
           hdu0.header.set('QUALITY',-999,'pixel quality flag')
       try:
               pc1 = str(pos_corr1[i])
               pc2 = str(pos_corr2[i])
               hdu0.header.set('POSCORR1',pc1,'[pix] column position correction')
               hdu0.header.set('POSCORR2',pc2,'[pix] row position correction')
       except:
               hdu0.header.set('POSCORR1','','[pix] column position correction')
               hdu0.header.set('POSCORR2','','[pix] row position correction')

       for j in range(len(cards2)):
           try:
               if cards2[j].keyword not in hdu0.header.keys():
                   hdu0.header.set(cards2[j].keyword, cards2[j].value, cards2[j].comment)
           except:
               pass

# Fix a few things
       hdu0.header.remove('PCOUNT')
       hdu0.header.remove('GCOUNT')
       hdu0.header.remove('TNULL4')
       hdu0.header.remove('NPIXSAP')
       hdu0.header.remove('NPIXMISS')
       hdu0.header.remove('TIERABSO')
       hdu0.header.remove('LC_START')
       hdu0.header.remove('LC_END')
       hdu0.header.remove('CDPP3_0')
       hdu0.header.remove('CDPP6_0')
       hdu0.header.remove('CDPP12_0')
       hdu0.header.remove('CROWDSAP')
       hdu0.header.remove('FLFRCSAP')
       hdu0.header.remove('DBCOLCO')
       hdu0.header.remove('DBTHRES')

       hdu0.header.set('OBJECT','M8','target')
#       hdu0.header['KEPLERID'] = UNDEFINED
#       hdu0.header['RA_OBJ'] = UNDEFINED
#       hdu0.header['DEC_OBJ'] = UNDEFINED
#       hdu0.header['LC_START'] = UNDEFINED
#       hdu0.header['LC_END'] = UNDEFINED
       hdu0.header['GAIN'] = UNDEFINED
       hdu0.header.cards['GAIN'].comment = 'double-valued for this superstamp'
       hdu0.header['READNOIS'] = UNDEFINED
       hdu0.header.cards['READNOIS'].comment = 'double-valued for this superstamp'
#       hdu0.header['NPIXSAP'] = UNDEFINED
#       hdu0.header['NPIXMISS'] = UNDEFINED
#       hdu0.header['TNULL'] = UNDEFINED
       hdu0.header['TSTART'] = time[i] - frametim/3600./24./2. * num_frm
       hdu0.header['TSTOP'] = time[i] + frametim/3600./24./2. * num_frm
# Do beginning of obs here:
       temptime = Time(time[i] - frametim/3600./24./2. * num_frm - 2400000.5-timecorr[i]+(0.25 + 0.62*(5-cards1['TIMSLICE'].value))/86400., format='mjd')
       temptime = str(temptime.datetime)
       hdu0.header['DATE-OBS'] = temptime.replace(' ','T')+'Z'

# Do end of obs here:
       temptime = Time(time[i] + frametim/3600./24./2. * num_frm - 2400000.5-timecorr[i]+(0.25 + 0.62*(5-cards1['TIMSLICE'].value))/86400., format='mjd')
       temptime = str(temptime.datetime)
       hdu0.header['DATE-END'] = temptime.replace(' ','T')+'Z'    

# Remove current WCS keywords

       hdu0.header.remove('CDELT1')
       hdu0.header.remove('CDELT2')
       hdu0.header.remove('PC1_1')
       hdu0.header.remove('PC1_2')
       hdu0.header.remove('PC2_1')
       hdu0.header.remove('PC2_2')

# Open oldoutfile here and get wcs info out of it...
       wcsfile = fits.open(oldoutfile,mode='readonly',memmap=True)
       cards1a = wcsfile[0].header.cards

       hdu0.header.set(cards1a['CRVAL1'][0],cards1a['CRVAL1'][1],cards1a['CRVAL1'][2])
       hdu0.header.set(cards1a['CRVAL2'][0],cards1a['CRVAL2'][1],cards1a['CRVAL2'][2])
       hdu0.header.set(cards1a['CRPIX1'][0],cards1a['CRPIX1'][1],cards1a['CRPIX1'][2])
       hdu0.header.set(cards1a['CRPIX2'][0],cards1a['CRPIX2'][1],cards1a['CRPIX2'][2])
       hdu0.header.set(cards1a['CD1_1'][0],cards1a['CD1_1'][1],cards1a['CD1_1'][2])
       hdu0.header.set(cards1a['CD1_2'][0],cards1a['CD1_2'][1],cards1a['CD1_2'][2])
       hdu0.header.set(cards1a['CD2_1'][0],cards1a['CD2_1'][1],cards1a['CD2_1'][2])
       hdu0.header.set(cards1a['CD2_2'][0],cards1a['CD2_2'][1],cards1a['CD2_2'][2])
       hdu0.header.set(cards1a['A_ORDER'][0],cards1a['A_ORDER'][1],cards1a['A_ORDER'][2])
       hdu0.header.set(cards1a['B_ORDER'][0],cards1a['B_ORDER'][1],cards1a['B_ORDER'][2])
       hdu0.header.set(cards1a['AP_ORDER'][0],cards1a['AP_ORDER'][1],cards1a['AP_ORDER'][2])
       hdu0.header.set(cards1a['BP_ORDER'][0],cards1a['BP_ORDER'][1],cards1a['BP_ORDER'][2])

# Get RA, Dec of image center
       wcs = WCS(oldoutfile)

       xmid = 109.5
       ymid = 67
       skypos = pixel_to_skycoord(xmid,ymid,wcs=wcs)
       ramid = skypos.ra.deg
       decmid = skypos.dec.deg

       hdu0.header['RA_OBJ'] = ramid
       hdu0.header['DEC_OBJ'] = decmid

       hdu0.header.cards['RA_OBJ'].comment = '[deg] right ascension of field center'
       hdu0.header.cards['DEC_OBJ'].comment = '[deg] declination of field center'

       hdu0.header.remove('XTENSION')
       hdu0.header.remove('EXTEND')
       hdu0.header.remove('NEXTEND')
       hdu0.header.remove('EXTVER')

# Add some history of file generation

       hdu0.header.set('HISTORY','This superstamp was generated from the following target pixel files:')
       for file in filenames:
         hdu0.header.set('HISTORY',' '+file[44:])
       hdu0.header.set('HISTORY','Target pixel file stamps stitched together according to 1CRV4P and 2CRV4P keywords')
       hdu0.header.set('HISTORY','Initial WCS solution provided by the astrometry.net suite.')
       hdu0.header.set('HISTORY','For more details, see http://astrometry.net.')
       hdu0.header.set('HISTORY','Fine tuning of WCS solution subsequently carried out using 2MASS star centroids and the IRAF CCMAP task.')

# write output file
#       outstr = fits.HDUList(hdu0)
#       outstr.writeto(outfile, checksum=True)

       hdu0.writeto(outfile,checksum=True)


if __name__ == '__main__':
  main()
