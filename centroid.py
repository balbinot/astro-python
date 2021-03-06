#! /usr/bin/env python
import numpy as np
import argparse
from astropy.modeling import models, fitting
import pyfits
from scipy.optimize import curve_fit,leastsq
from scipy import stats


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


def centroid_airy(data, coords, rad = 30, returnFit = False):
    if isinstance(data,str):
        data = pyfits.getdata(data)

    # Transpose x and y b/c reasons
    center_y,center_x = coords
    dslice = data[center_x-rad:center_x+rad,center_y-rad:center_y+rad]

    # Construct a grid of coordinates
    x,y = np.mgrid[0:dslice.shape[0],0:dslice.shape[1]]
    x -= dslice.shape[0]/2.
    y -= dslice.shape[1]/2.
                
    p_init = models.AiryDisk2D(np.max(dslice),0,0,rad)
    p = fitter(p_init,x,y,dslice)

    # Rescale coordinates to match data
    px = center_y + p.y_0
    py = center_x + p.x_0

    if returnFit:
        return px, py, p
    
    else:
        return px, py


def r2(p, x, y):
    return p[4] * (x - p[2])**2 + p[5] * (y - p[3])**2

def moffat(p, rsq, satlevel):
    return np.clip(p[0] + p[1] / np.abs(1 + rsq)**p[6], 0, satlevel)

def errfunc(p, x, y, f, satlevel):
    rsq = r2(p, x, y)   
    return (moffat(p, rsq, satlevel) - f) 

#def moffat(xy,scale,a,b):
#    x,y = xy
#    return scale*(b-1.0)/(np.pi*a**2)*(1.0+(x**2+y**2)/a**2)**(-b)

def centroid_moffat(data,coords, rad=30):
    if isinstance(data,str):
        data = pyfits.getdata(data)

    # Transpose x and y b/c reasons
    center_y,center_x = coords
    dslice = data[center_x-rad:center_x+rad,center_y-rad:center_y+rad]
    dim = np.max(dslice.shape)

    # Generate x,y grids
    x = np.linspace(0, dim - 1., dim) - dim // 2
    y = np.linspace(0, dim - 1., dim) - dim // 2
    x, y = np.meshgrid(x, y)

    dslice = np.reshape(dslice, (-1))
    x = np.reshape(x, (-1))
    y = np.reshape(y, (-1))

    # leastsq requires a saturation level?
    satlevel = stats.scoreatpercentile(dslice, 99)
    
    #p_init = [vertshift,scaling,xcent,ycent,alpha,beta]
    p_init = [np.min(dslice),np.max(dslice),center_y,center_x,5e-4,1,8]
    p,success = leastsq(errfunc,p_init[:],args=(x,y,dslice,satlevel))
    print p
    return (p[2],p[3])


        
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
