#!/usr/bin/env python
import numpy as np
import pyfits
#import glob
import matplotlib.pyplot as plt
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Median combine sky fibers and subtract from target images.')
    parser.add_argument('filelist',nargs='+',help='List of images to apply skysubtraction.')
    parser.add_argument('-sky',nargs='+',required=True,help='List of sky fibers to median-combine.')
    parser.add_argument('--patch',nargs='*',type=float,default=None,help='Wavelength(s) to patch over.')
    parser.add_argument('--width',default=15.0,type=float,dest='dlambda',help='Width to patch over.')

        
 #   skyfiles = glob.glob('*.sky.fits')
 #   image = '204.IRC10420.fits'

    args = parser.parse_args()

    skyspecs = []
    for filename in args.sky:
        f = pyfits.open(filename)
        spectrum = f[0].data[1]
        skyspecs.append(spectrum)

    medspec = np.median(skyspecs,axis=0)

    for image in args.filelist:
        if 'sky' in image:
            continue
        #print image
        f = pyfits.open(image)
        imspec = f[0].data[1]
        header = f[0].header
        #print header
        #exit()

        newspec = imspec - medspec

        newname = os.path.splitext(image)
        newname = ''.join([newname[0],'.Gskysub',newname[1]])


        #Patch over
        if args.patch is not None:
            
            # construct waves
            size = np.max(newspec.shape)
            lambda1 = header['crval1']
            lambda2 = np.round(size*header['cdelt1']+header['crval1'])

            waves = np.arange(lambda1, lambda2, header['cdelt1'])
            waves = waves[:size]

            newspec = newspec
            for lambda_cen in args.patch:
                rep_vals = np.median([y for x,y in zip(waves,newspec[0]) if x < (lambda_cen + args.dlambda) and x > (lambda_cen - args.dlambda)])

                newspec[0] = [y if (x > (lambda_cen + args.dlambda) or x < (lambda_cen - args.dlambda)) else rep_vals for x,y in zip(waves,newspec[0])]

        print 'Writing to %s...' % newname
        pyfits.writeto(newname,newspec,header,clobber=True)
        #pyfits.writeto('IRC10420_RED.fits',newspec,header,clobber=True)


    '''
    lambda1 = header['crval1']
    lambda2 = len(spectrum)*header['cdelt1']+header['crval1']
    waves = np.linspace(lambda1,lambda2,num=len(spectrum))

    plt.plot(waves,medspec,'r')
    plt.figure()
    plt.plot(waves,imspec,'b')

    plt.show()
    '''

if __name__ == '__main__':
    main()
