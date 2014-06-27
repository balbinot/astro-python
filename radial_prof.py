#! /usr/bin/env python
import argparse
import numpy as np
import pyfits

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
