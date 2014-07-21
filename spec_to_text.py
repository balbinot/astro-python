#! /usr/bin/env python
import argparse
import pyfits
import numpy as np
import os

def readSpectrum(filename):
    '''
    Read in a spectrum.  Return wavelengths, intensities, and header.
    '''
    data, header = pyfits.getdata(filename, 0, header=True)
    if len(data) != int(header['naxis1']):
        data = data[0]

    waves = np.arange(0,len(data))*header['CDELT1'] + header['CRVAL1']
  
    return (waves, data, header)


def main():
    parser = argparse.ArgumentParser(description='Copy first extension of input files to target directory.')
    
    parser.add_argument('filelist',nargs='+',help='List of files with 1D spectra.')

    parser.add_argument('-outdir',default='.',help="Target directory to copy files (default = './')")


    args = parser.parse_args()

    fcounts = 0
    for filename in args.filelist:
        waves,data,hdr = readSpectrum(filename)
        
        #outname
        outname = os.path.splitext(filename)[0] + '.txt'
        outname = os.path.join(args.outdir,outname)

        print 'Writing to %s' % outname
        f = open(outname,'w')
        for x,y in zip(waves,data):
            f.write('%f  %f\n' % (x,y))
        f.close()

        fcounts += 1

    print '%i spectra written' % fcounts

if __name__ == '__main__':
    main()
