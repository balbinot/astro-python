#! /usr/bin/env python
import numpy as np
import PyGuide
import argparse
import sys
import pyfits

ccdinfo = PyGuide.CCDInfo(1,1,1,1)

def centroid(item,coords,rad=40):
    if isinstance(item,str):
        item = pyfits.getdata(item)

    satmask = np.zeros_like(item,dtype=bool)
    satmask[coords] = True
    cent = PyGuide.centroid(item,mask=None,satMask=satmask,xyGuess=coords,rad=rad,ccdInfo=ccdinfo)
    xyctr = np.array(cent.xyCtr)
    return xyctr

def main():

    parser = argparse.ArgumentParser(description='Find centroid of object given initial guess.')

    parser.add_argument('data',nargs='+',help='List of filenames or numpy arrays.')
    parser.add_argument('coords',help='Initial guess of coordinates [x, y]')
    parser.add_argument('-r',type=int,default=40,help='Search radius (default is 40 pixels)')

                                
    args = parser.parse_args()

    return centroid(args.data,args.coords)


if __name__ == '__main__':
    sys.exit(main())
    
