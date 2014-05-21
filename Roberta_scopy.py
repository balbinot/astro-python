#! /usr/bin/env python
import argparse
import pyfits
import numpy as np
import os

def main():
    parser = argparse.ArgumentParser(description='Copy first extension of input files to target directory.')
    
    parser.add_argument('filelist',metavar='files',nargs='+',help='List of files with 1D spectra.')

    parser.add_argument('-target',metavar='target',required=True,help='Target directory to copy files.')

    args = parser.parse_args()

    if args.target[-1] != '/':
        args.target = ''.join([args.target,'/'])

   
    fcounts = 0
    for filename in args.filelist:
        f = pyfits.open(filename)

        newname = os.path.basename(filename).split('.')
        newname = '.'.join(newname[1:])
        newfitsname = os.path.join(args.target,newname)

        newfitsname = newfitsname.replace('.ms.fits','.fits')
        newtxtname = os.path.splitext(newfitsname)[0]+'.txt'
        
        data = f[0].data[0]
        hdr = f[0].header
        f.close()

        # get waves
        size = np.max(data.shape)
        lambda1 = hdr['crval1']
        lambda2 = np.round(size*hdr['cdelt1']+hdr['crval1'])

        wave = np.arange(lambda1, lambda2, hdr['cdelt1'])
        wave = wave[:size]
        
        print 'Writing to %s...' % (newfitsname)
        print 'Writing to %s...\n' % (newtxtname)

        pyfits.writeto(newfitsname,data,header=hdr,clobber=True)

        f = open(newtxtname,'w')
        for x,y in zip(wave,data[0]):
            f.write('%f  %f\n' % (x,y))
        f.close()
        fcounts += 1

    if fcounts == 1:
        print '1 spectrum written'
    else:
        print '%i spectra written' % (fcounts)


if __name__ == '__main__':
    main()

    
