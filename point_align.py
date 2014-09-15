#! /usr/bin/env python
import pyfits
from scipy.ndimage.interpolation import shift
from centroid import centroid
import argparse


def align(data,coords,rad):

    


def main():
    parser = argparse.ArgumentParser(description='Align input images given an initial guess.')
    
    parser.add_argument('filelist',nargs='+',help='List of filenames for alignment.')
    parser.add_argument('--coords',nargs=2,type=list,help='Initial guess of coordinates: x y')
    parser.add_argument('-r',type=int,default=40,help='Search radius (default is 40 pixels)')
    
    args = parser.parse_args()

    for filename in args.filelist:
        align(filename,args.coords,args.r)

    return 0

if __name__ == '__main__':
    main()
