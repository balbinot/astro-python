#! /usr/bin/env python
import argparse
import numpy as np
import pyfits
import matplotlib.pyplot as plt

def sector_mask(shape,center,radius,angle_range):
    """
    Return a boolean mask for a circular sector. The start/stop angles in  
    `angle_range` should be given in clockwise order.
    """

    x,y = np.ogrid[:shape[0],:shape[1]]
    cx,cy = center
    tmin,tmax = np.deg2rad(angle_range)

    # ensure stop angle > start angle
    if tmax < tmin:
            tmax += 2*np.pi

    # convert cartesian --> polar coordinates
    r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy)
    theta = np.arctan2(x-cx,y-cy) - tmin

    # wrap angles between 0 and 2*pi
    theta %= (2*np.pi)

    # circular mask
    circmask = r2 <= radius*radius

    # angular mask
    anglemask = theta <= (tmax-tmin)

    return circmask*anglemask

def sectorprofile(data, center=None, rad=15, binsize=30, interpnan = True):
    if isinstance(data, str):
        data = pyfits.getdata(data)

        
    y, x = np.indices(data.shape)
    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    #r = np.hypot(x - center[0], y - center[1])
    
    sector_prof = []
    bin_centers = []
    for start,stop in zip(np.arange(0,361,binsize),np.arange(binsize,361+binsize,binsize)):
        mask = sector_mask(data.shape,center,radius=rad,angle_range= (start,stop))
        sector = np.ma.MaskedArray(data,~mask)
        sector_prof.append(sector.sum()/sector.count())
        bin_centers.append(np.mean([start,stop]))

    
    return (bin_centers, sector_prof)

def main():
    parser = argparse.ArgumentParser(description='Return sector profile of target object')
    parser.add_argument('data',type=str,help='FITS file for profiling')
    parser.add_argument('-coords',type=float,nargs=2,default=None,help='x y coords of object. Default = center pixel of array')
    parser.add_argument('-bin',type=float,default=30,help='Bin angle (deg). Default=30')
    parser.add_argument('--interp',action='store_true',help='If specified, interpolate over NaN values')

    args = parser.parse_args

    profile = sectorprofile(*args)

    print zip(profile[0],profile[1])

    return 0

if __name__ == '__main__':
    main()
