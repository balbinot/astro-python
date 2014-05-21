#! /usr/bin/env python
import argparse
import pyfits
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Copy first extension of input files to target directory.')
    
    parser.add_argument('filelist',metavar='files',nargs='*',help='List of files with 1D spectra.')

    parser.add_argument('-target',metavar='target',help='Target directory to copy files.')

    parser.add_argument('-prefix',metavar='prefix',help='Filename prefix for output files.')
    
    args = parser.parse_args()

    if args.target[-1] != '/':
        args.target = ''.join([args.target,'/'])

    fcounts = 0
    for filename in args.filelist:
        if 'opptar' in filename or 'unused' in filename:
            continue

        f = pyfits.open(filename)
        newname = f[0].header['CATOBJ'].replace('(','-').replace(')','')
        newname = '%s.%i' % (newname, f[0].header['FIBER'])
        newfitsname = ''.join([args.target,args.prefix,newname,'.fits'])
        starname = ''.join([args.target,'star/',args.prefix,newname,'.txt'])
        
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
        print 'Writing to %s...\n' % (starname)
        pyfits.writeto(newfitsname,data,header=hdr,clobber=True)

        f = open(starname,'w')
        for x,y in zip(wave,data):
            f.write('%f  %f\n' % (x,y))
        f.close()
        fcounts += 1

    if fcounts == 1:
        print '1 spectrum written'
    else:
        print '%i spectra written' % (fcounts)


if __name__ == '__main__':
    main()

    
