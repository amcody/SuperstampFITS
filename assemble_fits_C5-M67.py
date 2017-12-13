#from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import pyfits

from pyfits import *
import glob

def main():

    prefix='M67'

    filenames = glob.glob('/Volumes/Work/Field_5/M67/TPFs/ktwo200*targ.fits')

    midfile = pyfits.open(filenames[28],mode='readonly',memmap=True)
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

# Images 503-505 have a header problem... as do 2615-2617
    for i in range(2618,3663):
#    for i in range(506,2614):
#    for i in range(0,3663):
#    for i in range(0,502):

# construct output primary extension
       outfile = prefix + '_BJD%.4f' % time[i] + '.fits'
       print outfile

       for fn in filenames:
           with pyfits.open(fn) as f:
               x = f[1].header['1CRV4P']
               y = f[1].header['2CRV4P']
               ch = f[0].header['CHANNEL']
               mod = f[0].header['MODULE']
               out = f[0].header['OUTPUT']
               alldata = f[1].data['FLUX']
               dim = alldata.shape
#               print x,y
               if (out == 1):
                 bigarr[y-187:y-137,x-801:x-801+dim[2]] = alldata[i]
               if (out == 2):
                 bigarr[y-187:y-137,1112-(x+dim[2])+1112-801:1112-x+1112-801] = np.fliplr(alldata[i])
#          -->  x: {801,1101}
#          -->  y: {187,537}
#              (50x50 and some 11x50 stamps)

       hdu0 = pyfits.PrimaryHDU(bigarr)


       for j in range(len(cards2)):
           try:
               if cards2[j].key not in hdu0.header.keys():
                   hdu0.header.update(cards2[j].key, cards2[j].value, cards2[j].comment)
           except:
               pass

       hdu0.header.remove('XTENSION')

       for k in range(len(cards1)):
            if (cards1[k].key not in hdu0.header.keys() and
                cards1[k].key[:4] not in ['TTYP','TFOR','TUNI','TDIS','TDIM','WCAX','1CTY',
                                          '2CTY','1CRP','2CRP','1CRV','2CRV','1CUN','2CUN',
                                          '1CDE','2CDE','1CTY','2CTY','1CDL','2CDL','11PC',
                                          '12PC','21PC','22PC','WCSN','TFIE']):
                hdu0.header.update(cards1[k].key, cards1[k].value, cards1[k].comment)

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
       hdu0.header.update('EXTNAME','IMAGE','name of extension')
       try:
           hdu0.header.update('TELAPSE',frametim * num_frm,'[s] elapsed time for exposure')
       except:
           hdu0.header.update('TELAPSE',-999,'[s] elapsed time for exposure')
       try:
           hdu0.header.update('LIVETIME',int_time * num_frm,'[s] TELASPE multiplied by DEADC')
       except:
           hdu0.header.update('LIVETIME',-999,'[s] TELASPE multiplied by DEADC')
       try:
           hdu0.header.update('EXPOSURE',int_time * num_frm,'[s] time on source')
       except:
           hdu0.header.update('EXPOSURE',-999,'[s] time on source')
       try:
           hdu0.header.update('MIDTIME',time[i],'[BJD] mid-time of exposure')
       except:
           hdu0.header.update('MIDTIME',-999,'[BJD] mid-time of exposure')
       try:
           hdu0.header.update('TIMECORR',timecorr[i],'[d] barycenter - timeslice correction')
       except:
           hd01.header.update('TIMECORR',-999,'[d] barycenter - timeslice correction')
       try:
           hdu0.header.update('CADENCEN',cadenceno[i],'unique cadence number')
       except:
           hdu0.header.update('CADENCEN',-999,'unique cadence number')
       try:
           hdu0.header.update('QUALITY',quality[i],'pixel quality flag')
       except:
           hdu0.header.update('QUALITY',-999,'pixel quality flag')
       try:
               pc1 = str(pos_corr1[i])
               pc2 = str(pos_corr2[i])
               hdu0.header.update('POSCORR1',pc1,'[pix] column position correction')
               hdu0.header.update('POSCORR2',pc2,'[pix] row position correction')
       except:
               hdu0.header.update('POSCORR1',-999,'[pix] column position correction')
               hdu0.header.update('POSCORR2',-999,'[pix] row position correction')

       hdu0.header.remove('XTENSION')

       for k in range(len(cards0)):
            try:
              if cards0[k].key not in hdu0.header.keys():
                  hdu0.header.update(cards0[k].key, cards0[k].value, cards0[k].comment)
              else:
                  hdu0.header.cards[cards0[k].key].comment = cards0[k].comment
            except:
              pass

       hdu0.header.remove('EXTEND')
       hdu0.header.remove('NEXTEND')
       hdu0.header.remove('EXTVER')

       hdu0.header.update('OBJECT','M67','target')


# write output file
       hdu0.writeto(outfile,checksum=True)

      
if __name__ == '__main__':
  main()
