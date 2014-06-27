#! /usr/bin/env python
import numpy as np
import argparse
from astropy.modeling import models, fitting
import pyfits

fitter = fitting.LevMarLSQFitter()

def centroid(data, coords, rad = 30, returnFit = False):
    if isinstance(data,str):
        data = pyfits.getdata(data)
        
    # Transpose x and y b/c reasons
    center_y,center_x = coords

    dslice = data[center_x-rad:center_x+rad,center_y-rad:center_y+rad]
    x,y = np.mgrid[0:dslice.shape[0],0:dslice.shape[1]]
    x -= dslice.shape[0]/2.
    y -= dslice.shape[1]/2.
                
    p_init = models.Gaussian2D(np.max(dslice),0,0,rad,rad)
    p = fitter(p_init,x,y,dslice)
    
    # Rescale coordinates to match data
    p.x_mean = center_y - p.x_mean
    p.y_mean = center_x - p.y_mean

    if returnFit:
        return p.x_mean.value, p.y_mean.value, p
    
    else:
        return p.x_mean.value, p.y_mean.value


def main():
    parser = argparse.ArgumentParser(description='Return central coordinates of object based on initial guess.')
    parser.add_argument('data',type=str,nargs=1,help='FITS file with object to be located')
    parser.add_argument('coords',type=int,nargs=2,help='x y coordinates of initial guess')
    parser.add_argument('-rad',type=int,default=30,help='Search radius [in pixels] from initial guess.')
    parser.add_argument('--fit',action='store_true',help='Return fitting params')

    args = parser.parse_args()

    center = centroid(*args)

    print '%s\t%f\t%f' % (args.data,center[0],center[1])
    if args.fit:
        print center[2]

    return 0

if __name__ == '__main__':
    main()
