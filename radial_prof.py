#! /usr/bin/env python
import argparse
import numpy as np
import pyfits
import matplotlib.pyplot as plt

def radialprofile(data, center=None, binsize=1, interpnan = True):
    # get indices of data
    if isinstance(data, str):
        data = pyfits.getdata(data)
    y, x = np.indices(data.shape)
    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    r = np.hypot(x - center[0], y - center[1])

    nbins = int(np.round(r.max() / binsize) + 1)
    maxbins = nbins * binsize
    bins = np.linspace(0,maxbins,nbins+1)
    bin_centers = (bins[1:] + bins[:-1])/2.0

    radial_prof = np.histogram(r,bins, weights = data)[0] / np.histogram(r,bins, weights=np.ones(data.shape))[0]

    if interpnan:
        radial_prof = np.interp(bin_centers,bin_centers[radial_prof==radial_prof],radial_prof[radial_prof==radial_prof])    

    return (bin_centers, radial_prof)

def interpol(x,y,edgeclip=False):
    yp = np.interp(x,x[y==y],y[y==y])
    if edgeclip:
        return yp[:-1]
    else:
        return yp

def radialForTerry(data,center=None,binsize=1,interpnan=True):
    # get indices of data
    if isinstance(data, str):
        data = pyfits.getdata(data)
    y, x = np.indices(data.shape)
    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    r = np.hypot(x - center[0], y - center[1])

    nbins = int(np.round(r.max() / binsize) + 1)
    maxbins = nbins * binsize
    bins = np.linspace(0,maxbins,nbins+1)
    bin_centers = (bins[1:] + bins[:-1])/2.0

    y = np.histogram(r,bin_centers, weights = data)[0] / np.histogram(r,bin_centers, weights=np.ones(data.shape))[0]

    n = np.histogram(r,bins=bin_centers)[0]
    sy = np.histogram(r,weights=data,bins=bin_centers)[0]
    sy2 = np.histogram(r,weights=data*data,bins=bin_centers)[0]
    
    if interpnan:
        y = interpol(bin_centers,y,edgeclip=True)

    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)
    #skip nans
    idx = np.where(y==y)
    bin_centers = np.array(bin_centers)[idx]
    y = y[idx]
    std = std[idx]
    std = interpol(bin_centers[idx],std)/2

    return bin_centers, y, std


    
def radialSeq(data,binseq,center=None):
    # get indices of data
    if isinstance(data, str):
        data = pyfits.getdata(data)
    y, x = np.indices(data.shape)
    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    r = np.hypot(x - center[0], y - center[1])

    n = np.histogram(r,bins=binseq)[0]
    sy = np.histogram(r,weights=data,bins=binseq)[0]
    
    sy2 = np.histogram(r,weights=data*data,bins=binseq)[0]
    #sy2 = interpol(binseq,sy2)
    
    y = sy/np.histogram(r,weights=np.ones(data.shape),bins=binseq)[0]
    y = interpol(binseq,y,edgeclip=True)
    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)

    #skip nans
    idx = np.where(y==y)
    x = np.array(binseq)[idx]
    y = y[idx]
    std = std[idx]
    #mean = mean[idx]  # same as y!

    # interp over 0
    '''
    idx = np.where(y==y)
    x = np.array(binseq)[idx]
    y = y[idx]
    std = np.interp(binseq[0:-1],x,std)
    #'''

    std = interpol(x[idx],std)/2
    #plt.errorbar(x,y,yerr=std,fmt='o')
    
    return (x,y,std)
    
def radialInt(data,center=None):
    # get indices of data
    if isinstance(data, str):
        data = pyfits.getdata(data)
    y, x = np.indices(data.shape)
    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    r = np.hypot(x - center[0], y - center[1])
    r = r.astype(np.int)
    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialp = tbin / nr
    print nr
    
    return radialp
    
    
def radialAllPix(data, center=None):
    # get indices of data
    if isinstance(data, str):
        data = pyfits.getdata(data)
    y, x = np.indices(data.shape)
    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    r = np.hypot(x - center[0], y - center[1])
    yarr = data.flatten()
    xarr = r.flatten()
    
    return xarr,yarr
    '''
    radDict = {}
    for x,y in zip(xarr,yarr):
        if x in radDict:
            radDict[x].append(y)
        else:
            radDict[x] = [y]

    xrad = []
    yrad = []
    for k,v in radDict.items():
        xrad.append(k)
        yrad.append(np.mean(v))
        
    return xrad,yrad
    #'''

def radialAverage(image, center=None, stddev=False, returnAz=False, return_naz=False, 
                  binsize=1.0, weights=None, steps=False, interpnan=False, left=None, right=None,
                  mask=None, symmetric=None ):
    """
    Calculate the radially averaged azimuthal profile.
    (this code has not been optimized; it could be speed boosted by ~20x)

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
    None, which then uses the center of the image (including 
    fractional pixels).
    stddev - if specified, return the radial standard deviation instead of the average
    returnAz - if specified, return (azimuthArray,azimuthal_profile)
    return_naz   - if specified, return number of pixels per azimuth *and* azimuth
    binsize - size of the averaging bin.  Can lead to strange results if
    non-binsize factors are used to specify the center and the binsize is
    too large
    weights - can do a weighted average instead of a simple average if this keyword parameter
    is set.  weights.shape must = image.shape.  weighted stddev is undefined, so don't
    set weights and stddev.
    steps - if specified, will return a double-length bin array and azimuthal
    profile so you can plot a step-form azimuthal profile (which more accurately
    represents what's going on)
    interpnan - Interpolate over NAN values, i.e. bins where there is no data?
    left,right - passed to interpnan; they set the extrapolated values
    mask - can supply a mask (boolean array same size as image with True for OK and False for not)
    to average over only select data.

    If a bin contains NO DATA, it will have a NAN value because of the
    divide-by-sum-of-weights component.  I think this is a useful way to denote
    lack of data, but users let me know if an alternative is prefered...
    
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if center is None:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])
    theta = np.arctan2(x - center[0], y - center[1])
    theta[theta < 0] += 2*np.pi
    theta_deg = theta*180.0/np.pi
    maxangle = 360

    if weights is None:
        weights = np.ones(image.shape)
    elif stddev:
        raise ValueError("Weighted standard deviation is not defined.")

    if mask is None:
        # mask is only used in a flat context
        mask = np.ones(image.shape,dtype='bool').ravel()
    elif len(mask.shape) > 1:
        mask = mask.ravel()

    # allow for symmetries
    if symmetric == 2:
        theta_deg = theta_deg % 90
        maxangle = 90
    elif symmetric == 1:
        theta_deg = theta_deg % 180
        maxangle = 180

    # the 'bins' as initially defined are lower/upper bounds for each bin
    # so that values will be in [lower,upper)  
    nbins = int(np.round(maxangle / binsize))
    maxbin = nbins * binsize
    bins = np.linspace(0,maxbin,nbins+1)
    # but we're probably more interested in the bin centers than their left or right sides...
    bin_centers = (bins[1:]+bins[:-1])/2.0

    # Find out which azimuthal bin each point in the map belongs to
    whichbin = np.digitize(theta_deg.flat,bins)

    # how many per bin (i.e., histogram)?
    # there are never any in bin 0, because the lowest index returned by digitize is 1
    nr = np.bincount(whichbin)[1:]

    # recall that bins are from 1 to nbins (which is expressed in array terms by arange(nbins)+1 or xrange(1,nbins+1) )
    # azimuthal_prof.shape = bin_centers.shape
    if stddev:
        azimuthal_prof = np.array([image.flat[mask*(whichbin==b)].std() for b in xrange(1,nbins+1)])
    else:
        azimuthal_prof = np.array([(image*weights).flat[mask*(whichbin==b)].sum() / weights.flat[mask*(whichbin==b)].sum() for b in xrange(1,nbins+1)])

    #import pdb; pdb.set_trace()

    if interpnan:
        azimuthal_prof = np.interp(bin_centers,
                                   bin_centers[azimuthal_prof==azimuthal_prof],
                                   azimuthal_prof[azimuthal_prof==azimuthal_prof],
                                   left=left,right=right)

    if steps:
        xarr = np.array(zip(bins[:-1],bins[1:])).ravel() 
        yarr = np.array(zip(azimuthal_prof,azimuthal_prof)).ravel() 
        return xarr,yarr
    elif returnAz: 
        return bin_centers,azimuthal_prof
    elif return_naz:
        return nr,bin_centers,azimuthal_prof
    else:
        return azimuthal_prof


def main():
    parser = argparse.ArgumentParser(description='Return azimuthally-averaged radially profile of target object')
    parser.add_argument('data',type=str,help='FITS file for profiling')
    parser.add_argument('-coords',type=float,nargs=2,default=None,help='x y coords of object. Default = center pixel of array')
    parser.add_argument('-bin',type=float,default=1,help='Binsize. Default=1')
    parser.add_argument('--interp',action='store_true',help='If specified, interpolate over NaN values')

    args = parser.parse_args

    profile = radialprofile(*args)

    print zip(profile[0],profile[1])

    return 0

if __name__ == '__main__':
    main()
