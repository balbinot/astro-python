#! /usr/bin/env python
import pyregion
import argparse
import ds9
from numpy import genfromtxt


def main():
    parser = argparse.ArgumentParser(description='Send regions to ds9')
    parser.add_argument('infile',type=str,help='[filename x y size]')
    

    args = parser.parse_args()
    
    d = ds9.ds9(target='pyds9')

    data = genfromtxt(args.infile,names=['fname','x','y','rad'],dtype=['a100','f8','f8','f8'],autostrip=True)

    print 'Printing regions to %i frames' % len(data['fname'])
    for idx,dat in enumerate(data):
        d.set('frame %i' % (idx+1))
        d.set('file %s' % dat['fname'])
        
        d.set('zscale')
        d.set('zoom to fit')

        d.set('regions command  "image; circle %f %f %f"' % (dat['x'],dat['y'],dat['rad']))
    d.set('match scale')
    d.set('match colorbar')
    d.set('frame 1')
        

    return 0


if __name__ == '__main__':
    main()
