import pylab, numpy, pyfits, scipy
from pylab import *
from numpy import *
#from pyfits import *
from astropy.io import fits
from pyraf import iraf
from iraf import kepler
import kepio, kepmsg, kepkey, kepfit, keparray, kepfunc, kepstat
import sys, time, re, math, glob
from scipy import interpolate, optimize, ndimage, stats
from scipy.optimize import fmin_powell
from scipy.interpolate import RectBivariateSpline, interp2d
from scipy.ndimage import interpolation
from scipy.ndimage.interpolation import shift, rotate
from scipy.stats import nanmean
from astropy.nddata.utils import Cutout2D

# -----------------------------------------------------------
# core code

def kepprf_AMC(infile,rownum,columns,rows,fluxes,border,background,focus,prfdir,xtol,ftol,verbose,logfile): 

# input arguments

    status = 0
    seterr(all="ignore") 

# open FITS file and get header info + data

    instr = fits.open(infile,mode='readonly',memmap=True)

    crval1p = instr[0].header['CRVAL1P']
    crval2p = instr[0].header['CRVAL2P']

# construct inital guess vector for fit 

    if status == 0:
        guess = []
        try:
            f = fluxes.strip().split(',')
            x = columns.strip().split(',')
            y = rows.strip().split(',')
            for i in xrange(len(f)):
                f[i] = float(f[i])
        except:
            f = fluxes
            x = columns
            y = rows
        nsrc = len(f)
        for i in xrange(nsrc):
            try:
                guess.append(float(f[i]))
            except:
                message = 'ERROR -- KEPPRF: Fluxes must be floating point numbers'
                status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            if len(x) != nsrc or len(y) != nsrc:
                message = 'ERROR -- KEPFIT:FITMULTIPRF: Guesses for rows, columns and '
                message += 'fluxes must have the same number of sources'
                status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            for i in xrange(nsrc):
                try:
                    guess.append(float(x[i])+crval1p)
                except:
                    message = 'ERROR -- KEPPRF: Columns must be floating point numbers'
                    status = kepmsg.err(logfile,message,verbose)
        if status == 0:
            for i in xrange(nsrc):
                try:
                    guess.append(float(y[i])+crval2p)
                except:
                    message = 'ERROR -- KEPPRF: Rows must be floating point numbers'
                    status = kepmsg.err(logfile,message,verbose)
        if status == 0 and background:
            if border == 0:
                guess.append(0.0)
            else:
                for i in range((border+1)*2):
                    guess.append(0.0)
        if status == 0 and focus:
            guess.append(1.0); guess.append(1.0); guess.append(0.0)

# Get data from image; then close it

    fluxpixels = instr[0].data[:]
    errpixels = np.array([sqrt(val) for val in fluxpixels])

    module = instr[0].header['MODULE']    
    output = instr[0].header['OUTPUT']

    instr.close()

    x[0] = str(int(x[0])+crval1p)
    y[0] = str(int(y[0])+crval2p)
    column = int(x[0]) - 3
    row = int(y[0]) - 3

##    xdim = 11
##    ydim = 11
##    npix = 121

    xdim = 7
    ydim = 7
    npix = 49

# construct input pixel image

    if status == 0:
        flux = fluxpixels[row-crval2p:row-crval2p+ydim,column-crval1p:column-crval1p+xdim]
        ferr = errpixels[row-crval2p:row-crval2p+ydim,column-crval1p:column-crval1p+xdim]

        isize = numpy.shape(flux)[0]
        jsize = numpy.shape(flux)[1]
        flux = numpy.reshape(flux,(isize*jsize))
        ferr = numpy.reshape(ferr,(isize*jsize))

        DATx = arange(column,column+xdim)
        DATy = arange(row,row+ydim)
#        if numpy.nanmin > 420000.0: flux -= 420000.0

# image scale and intensity limits of pixel data

    if status == 0:
        n = 0
        DATimg = empty((ydim,xdim))
        ERRimg = empty((ydim,xdim))
        for i in range(ydim):
            for j in range(xdim):
                DATimg[i,j] = flux[n]
                ERRimg[i,j] = ferr[n]
                n += 1

# determine suitable PRF calibration file

    if status == 0:
        if int(module) < 10:
            prefix = 'kplr0'
        else:
            prefix = 'kplr'
        prfglob = prfdir + '/' + prefix + str(module) + '.' + str(output) + '*' + '_prf.fits'
        try:
            prffile = glob.glob(prfglob)[0]
        except:
            message = 'ERROR -- KEPPRF: No PRF file found in ' + prfdir
            status = kepmsg.err(logfile,message,verbose)

# read PRF images

    if status == 0:
        prfn = [0,0,0,0,0]
        crpix1p = numpy.zeros((5),dtype='float32')
        crpix2p = numpy.zeros((5),dtype='float32')
        crval1p = numpy.zeros((5),dtype='float32')
        crval2p = numpy.zeros((5),dtype='float32')
        cdelt1p = numpy.zeros((5),dtype='float32')
        cdelt2p = numpy.zeros((5),dtype='float32')
        for i in range(5):
            prfn[i], crpix1p[i], crpix2p[i], crval1p[i], crval2p[i], cdelt1p[i], cdelt2p[i], status \
                = kepio.readPRFimage(prffile,i+1,logfile,verbose) 
        prfn = array(prfn)
        PRFx = arange(0.5,shape(prfn[0])[1]+0.5)
        PRFy = arange(0.5,shape(prfn[0])[0]+0.5)
        PRFx = (PRFx - size(PRFx) / 2) * cdelt1p[0]
        PRFy = (PRFy - size(PRFy) / 2) * cdelt2p[0]

# interpolate the calibrated PRF shape to the target position

    if status == 0:
        prf = zeros(shape(prfn[0]),dtype='float32')
        prfWeight = zeros((5),dtype='float32')
        for i in xrange(5):
            prfWeight[i] = sqrt((column - crval1p[i])**2 + (row - crval2p[i])**2)
            if prfWeight[i] == 0.0:
                prfWeight[i] = 1.0e-6
            prf = prf + prfn[i] / prfWeight[i]
        prf = prf / nansum(prf) / cdelt1p[0] / cdelt2p[0]

# location of the data image centered on the PRF image (in PRF pixel units)

    if status == 0:
        prfDimY = int(ydim / cdelt1p[0])
        prfDimX = int(xdim / cdelt2p[0])
        PRFy0 = (shape(prf)[0] - prfDimY) / 2
        PRFx0 = (shape(prf)[1] - prfDimX) / 2

# interpolation function over the PRF

    if status == 0:
        splineInterpolation = scipy.interpolate.RectBivariateSpline(PRFx,PRFy,prf)

# construct mesh for background model

    if status == 0 and background:
        bx = numpy.arange(1.,float(xdim+1))
        by = numpy.arange(1.,float(ydim+1))
        xx, yy = numpy.meshgrid(numpy.linspace(bx.min(), bx.max(), xdim),
                                numpy.linspace(by.min(), by.max(), ydim))

# fit PRF model to pixel data

    if status == 0:
        if focus and background:
            args = (DATx,DATy,DATimg,ERRimg,nsrc,border,xx,yy,splineInterpolation,float(x[0]),float(y[0]))
            ans = fmin_powell(kepfunc.PRFwithFocusAndBackground,guess,args=args,xtol=xtol,
                              ftol=ftol,disp=False)
        elif focus and not background:
            args = (DATx,DATy,DATimg,ERRimg,nsrc,splineInterpolation,float(x[0]),float(y[0]))
            ans = fmin_powell(kepfunc.PRFwithFocus,guess,args=args,xtol=xtol,
                              ftol=ftol,disp=False)                    
        elif background and not focus:
            args = (DATx,DATy,DATimg,ERRimg,nsrc,border,xx,yy,splineInterpolation,float(x[0]),float(y[0]))
            ans = fmin_powell(kepfunc.PRFwithBackground,guess,args=args,xtol=xtol,
                              ftol=ftol,disp=False)
        else:
            args = (DATx,DATy,DATimg,ERRimg,nsrc,splineInterpolation,float(x[0]),float(y[0]))
            ans = fmin_powell(kepfunc.PRF,guess,args=args,xtol=xtol,
                              ftol=ftol,disp=False)

# pad the PRF data if the PRF array is smaller than the data array 

    if status == 0:
        flux = []; OBJx = []; OBJy = []
        PRFmod = numpy.zeros((prfDimY,prfDimX))
        if PRFy0 < 0 or PRFx0 < 0.0:
            PRFmod = numpy.zeros((prfDimY,prfDimX))
            superPRF = zeros((prfDimY+1,prfDimX+1))
            superPRF[abs(PRFy0):abs(PRFy0)+shape(prf)[0],abs(PRFx0):abs(PRFx0)+shape(prf)[1]] = prf
            prf = superPRF * 1.0
            PRFy0 = 0
            PRFx0 = 0

# rotate the PRF model around its center

        if focus:
            angle = ans[-1]
            prf = rotate(prf,-angle,reshape=False,mode='nearest')

# iterate through the sources in the best fit PSF model

        for i in range(nsrc):
            flux.append(ans[i])
            OBJx.append(ans[nsrc+i])
            OBJy.append(ans[nsrc*2+i]) 

# calculate best-fit model

            y = (OBJy[i]-mean(DATy)) / cdelt1p[0]
            x = (OBJx[i]-mean(DATx)) / cdelt2p[0]
            prfTmp = shift(prf,[y,x],order=3,mode='constant')
            prfTmp = prfTmp[PRFy0:PRFy0+prfDimY,PRFx0:PRFx0+prfDimX]
            PRFmod = PRFmod + prfTmp * flux[i]
            wx = 1.0
            wy = 1.0
            angle = 0
            b = 0.0

# write out best fit parameters

            if verbose:
                txt = 'Flux = %10.2f e-/s ' % flux[i]
                txt += 'X = %9.4f pix ' % OBJx[i]
                txt += 'Y = %9.4f pix ' % OBJy[i]
                kepmsg.log(logfile,txt,True)

        if background:
            bterms = border + 1
            if bterms == 1:
                b = ans[nsrc*3]
            else:
                bcoeff = array([ans[nsrc*3:nsrc*3+bterms],ans[nsrc*3+bterms:nsrc*3+bterms*2]]) 
                bkg = kepfunc.polyval2d(xx,yy,bcoeff)
                b = nanmean(bkg.reshape(bkg.size))
            txt = '\n   Mean background = %.2f e-/s' % b
        if verbose and background:
              kepmsg.log(logfile,txt,True)
        if focus:
            wx = ans[-3]
            wy = ans[-2]
            angle = ans[-1]
        if verbose and focus:
            if not background: kepmsg.log(logfile,'',True)
            kepmsg.log(logfile,' X/Y focus factors = %.3f/%.3f' % (wx,wy),True)
            kepmsg.log(logfile,'PRF rotation angle = %.2f deg' % angle,True)

# measure flux fraction and contamination

    if status == 0:
        PRFall = kepfunc.PRF2DET(flux,OBJx,OBJy,DATx,DATy,wx,wy,angle,splineInterpolation)
        PRFone = kepfunc.PRF2DET([flux[0]],[OBJx[0]],[OBJy[0]],DATx,DATy,wx,wy,angle,splineInterpolation)

# constuct model PRF in detector coordinates

    if status == 0:
        PRFfit = PRFall + 0.0
        if background and bterms == 1:
            PRFfit = PRFall + b
        if background and bterms > 1:
            PRFfit = PRFall + bkg

# calculate residual of DATA - FIT

    if status == 0:
        PRFres = DATimg - PRFfit
        FLUXres = numpy.nansum(PRFres) / npix
    
# calculate the sum squared difference between data and model

    if status == 0:
        Pearson = abs(numpy.nansum(numpy.square(DATimg - PRFfit) / PRFfit))
        Chi2 = numpy.nansum(numpy.square(DATimg - PRFfit) / numpy.square(ERRimg))
        DegOfFreedom = npix - len(guess) - 1
        if verbose:
         try:
            kepmsg.log(logfile,'\n       Residual flux = %.2f e-/s' % FLUXres,True)
            kepmsg.log(logfile,'Pearson\'s chi^2 test = %d for %d dof' % (Pearson,DegOfFreedom),True)
         except:
            pass
         kepmsg.log(logfile,'          Chi^2 test = %d for %d dof' % (Chi2,DegOfFreedom),True)

        
    result = np.array([flux[0], OBJx[0], OBJy[0]])

    return result

