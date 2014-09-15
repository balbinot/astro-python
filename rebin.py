#! /usr/bin/env python
import argparse
import pyfits
import numpy as np
import os

def rebin(filename,func,b):
    if isinstance(filename,str):
        data = pyfits.getdata(filename)
    else:
        data = filename

    # num of rows and cols
    M = data.shape[0]
    N = data.shape[1]

    m = M / b
    n = N / b

    bin_im = np.zeros((m,n))

    # loop over rows and columns of binned image, filling each one pixel
    #  with the sum of the corresponding bin in full-size images
    for y in np.arange(0, bin_im.shape[0]):
        for x in np.arange(0, bin_im.shape[1]):
            # each (y,x) location in the binned image comes from a b x b box
            # on the full-size image.  In that box some number of
            # pixels from 0 up to b**2 can be NaNs.  If at least half
            # of them are finite (not NaNs), use the sum of the finite
            # values.  If more than half are not finite, that means this
            # b x b box covers too much of a masked region (or region of
            # many bad pixels) and therefore will be written off as not
            # useable (make it all NaNs).
            by = b*y ; bx = b*x
            goodvals = [val for val in data[by:by+b, bx:bx+b].ravel()] # if np.isfinite(val)]
            
            if func == 'sum':
                bin_im[y,x] = np.nansum(goodvals)
            elif func == 'median':
                bin_im[y,x] = np.nanmedian(goodvals)
            else:
                bin_im[y,x] = np.nanmean(goodvals)
                
                #if len(goodvals) >= (b*b)/2+1:
                #    bin_im[y,x] = np.median(goodvals)
                #else:
                #    bin_im[y,x] = np.nan

    hdr = pyfits.getheader(filename)
    hdr['BINSIZE']  = (b,'num pixels binned')

    newname = os.path.splitext(filename)
    newname = ''.join([newname[0],'_bin',newname[1]])

    print 'Writing to %s' % newname
    pyfits.writeto(newname,bin_im,header=hdr,clobber=True)
    return newname
            
def main():
    parser = argparse.ArgumentParser(description='Rebin image(s) with input shape.')
    parser.add_argument('filelist',nargs='+',help='List of images to rebin.')

    parser.add_argument('-size',type=int,required=True,help='Size of box to bin')
    parser.add_argument('--func',choices=['sum','median','mean'],default='sum',help='Binning function (default = sum).')

    args = parser.parse_args()

    print 'Rebinning files to %i x %i bins' % (args.size,args.size)
    
    for filename in args.filelist:
        rebin(filename,args.func,args.size)

    return 0

if __name__ == '__main__':
    main()
