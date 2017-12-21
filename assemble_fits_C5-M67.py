#from __future__ import division, print_function

import numpy as np
from astropy.io import fits
from astropy.io.fits.card import Undefined, UNDEFINED
import datetime
from astropy.time import Time
import glob

def main():

    prefix='M67'

    filenames = glob.glob('/Volumes/Work/Field_5/M67/TPFs/ktwo200*targ.fits')

    midfile = fits.open(filenames[0],mode='readonly',memmap=True)
    cards0 = midfile[0].header.cards
    cards1 = midfile[1].header.cards
    cards2 = midfile[2].header.cards

    time = midfile[1].data.field('TIME')[:] + 2454833.0
    timecorr = midfile[1].data.field('TIMECORR')[:]
    cadenceno = midfile[1].data.field('CADENCENO')[:]
    quality = midfile[1].data.field('QUALITY')[:]
    cosmic_rays = midfile[1].data.field('COSMIC_RAYS')[:]
    pos_corr1 = midfile[1].data.field('POS_CORR1')[:]
    pos_corr2 = midfile[1].data.field('POS_CORR2')[:]

    bigarr = np.zeros([400,400])
    bigarr[:] = np.nan
    bigarr = bigarr.astype('float32')

    for i in range(0,3663):

# only continue if there is valid data at this timestamp
     if (~np.isnan(time[i]) and len(np.where((midfile[1].data['FLUX'])[i] != 0.0)[0]) != 0 and len(np.where((~np.isnan(midfile[1].data['FLUX'])[i]))[0]) != 0):

# construct output primary extension
       outfile = prefix + '_BJD%.4f' % time[i] + '.fits'
       print outfile

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
                 bigarr[y-187:y-137,x-801:x-801+dim[2]] = alldata[i]
               if (out == 2):
                 bigarr[y-187:y-137,1112-(x+dim[2])+1112-801:1112-x+1112-801] = np.fliplr(alldata[i])
     
       hdu0 = fits.PrimaryHDU(bigarr)

# add in primary keywords
       for n in range(len(cards0)):
            try:
                if cards0[n].keyword not in hdu0.header.keys():
                    hdu0.header[cards0[n].keyword] = (cards0[n].value,
                                                      cards0[n].comment)
                else:
                    hdu0.header.cards[cards0[n].keyword].comment = cards0[n].comment
            except:
                pass

# add additional keywords
       for k in range(len(cards1)):
            if (cards1[k].keyword not in hdu0.header.keys() and
                cards1[k].keyword[:4] not in ['TTYP','TFOR','TUNI','TDIS','TDIM','WCAX','1CTY',
                                          '2CTY','1CRP','2CRP','1CRV','2CRV','1CUN','2CUN',
                                          '1CDE','2CDE','1CTY','2CTY','1CDL','2CDL','11PC',
                                          '12PC','21PC','22PC','WCSN','TFIE','XTEN','EXTN',
                                          'PCOU','GCOU','TNUL','INHE']):
                hdu0.header.set(cards1[k].keyword, cards1[k].value, cards1[k].comment)

# and a few more keywords:
       for j in range(len(cards2)):
           try:
               if cards2[j].keyword not in hdu0.header.keys():
                   hdu0.header.set(cards2[j].keyword, cards2[j].value, cards2[j].comment)
           except:
               pass


# pull some additional information out of the TPF headers

       try:
           int_time = cards1['INT_TIME'].value
       except:
           print 'WARNING -- cannot find INT_TIME keyword'
       try:
           frametim = cards1['FRAMETIM'].value
       except:
           print 'WARNING -- cannot find FRAMETIM keyword'
       try:
           num_frm = cards1['NUM_FRM'].value
       except:
           print 'WARNING -- cannot find NUM_FRM keyword'
       try:
           hdu0.header.set('TELAPSE',frametim * num_frm,'[s] elapsed time for exposure')
       except:
           hdu0.header.set('TELAPSE',-999,'[s] elapsed time for exposure')
       try:
           hdu0.header.set('LIVETIME',int_time * num_frm,'[s] TELASPE multiplied by DEADC')
       except:
           hdu0.header.set('LIVETIME',-999,'[s] TELASPE multiplied by DEADC')
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
           hd01.header.set('TIMECORR',-999,'[d] barycenter - timeslice correction')
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
               hdu0.header.set('POSCORR1',-999,'[pix] column position correction')
               hdu0.header.set('POSCORR2',-999,'[pix] row position correction')

# Edit and delete some keywords
       hdu0.header['DATE'] = str(datetime.date.today())
       hdu0.header['CREATOR'] = 'Ann Marie Cody'
       hdu0.header.cards['CREATOR'].comment = 'file creator'
       hdu0.header.set('OBJECT','M67','target')
       hdu0.header['CHANNEL'] = UNDEFINED
       hdu0.header.cards['CHANNEL'].comment = 'double-valued for this superstamp'
       hdu0.header['OUTPUT'] = UNDEFINED
       hdu0.header.cards['OUTPUT'].comment = 'double-valued for this superstamp'

       hdu0.header.remove('PCOUNT')
       hdu0.header.remove('GCOUNT')
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
       hdu0.header.remove('EXTEND')
       hdu0.header.remove('NEXTEND')
       hdu0.header.remove('EXTVER')
       hdu0.header.remove('XTENSION')
       hdu0.header.remove('INHERIT')
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
#       hdu0.header.set('NEXTEND', 1, 'number of extensions')

       hdu0.header.set('OBJECT','M67','target')
       hdu0.header['TSTART'] = time[i] - frametim/3600./24./2. * num_frm
       hdu0.header['TSTOP'] = time[i] + frametim/3600./24./2. * num_frm
       hdu0.header.cards['TSTART'].comment = 'observation start time in BJD'
       hdu0.header.cards['TSTOP'].comment = 'observation start time in BJD'

# Calculate TIME-OBS from BJD time:
       temptime = Time(time[i] - frametim/3600./24./2. * num_frm - 2400000.5-timecorr[i]+(0.25 + 0.62*(5-cards1['TIMSLICE'].value))/86400., format='mjd')
       temptime = str(temptime.datetime)
       hdu0.header['DATE-OBS'] = temptime.replace(' ','T')+'Z'

# Calculate DATE-END:
       temptime = Time(time[i] + frametim/3600./24./2. * num_frm - 2400000.5-timecorr[i]+(0.25 + 0.62*(5-cards1['TIMSLICE'].value))/86400., format='mjd')
       temptime = str(temptime.datetime)
       hdu0.header['DATE-END'] = temptime.replace(' ','T')+'Z'


# write output file
       hdu0.writeto(outfile,checksum=True)

if __name__ == '__main__':
  main()
