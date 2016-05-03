#! /usr/bin/env python
import ds9
import numpy as np
import argparse
import sys
import pyfits
import glob



def load_filelist(d, filelist,regions=None,loadregions=None,scale='zscale'):
    for idx,f in enumerate(filelist):
        d.set("frame %i" % (idx+1))
        d.set("file %s" % f)

        if regions is not None:
            if loadregions == 0 or loadregions-1 == idx:
                d.set("regions load %s" % regions)
        d.set(scale)
        d.set("zoom to fit")

    d.set("match colorbar")
    d.set("match scale")
    d.set("frame 1")
    
    return d
        
##########
### pydisplay
##########
def pydisplay(datalist,regions=None,loadregions=0,target=None,scale='zscale'):
    if target == None:
        target = 'pyds9'
        d = ds9.ds9(target=target) #target='DS9:pyds9',verify=False)
        
        '''
        if ds9.ds9_targets() != None:
            target = ds9.ds9_targets()[0].split()[1]
            d = ds9.ds9(target=target)
        
        else:
            d = ds9.ds9(target='DS9:pyds9',verify=False)
        '''

    else:
        d = ds9.ds9(target=target, wait=15)



    # if pydisplay(), return instance
    if datalist==None:
        return d
    #d = ds9.ds9()  #target='pyds9',start=True)
    #try:
    #except:
    #    d = ds9.ds9('pyds9',wait=20)

    # If filelist is just a single file string
    if isinstance(datalist,str):
        
        # Check if regular expression
        if '*' in datalist:
            datalist = glob.glob(datalist)
            return load_filelist(d,datalist,regions=None,loadregions=0,scale='zscale')
        
        d.set("file %s" % datalist)
        d.set(scale)
        d.set("zoom to fit")
        if regions is not None:
            d.set("regions load %s" % regions)
        return target

    # If datalist is a single array
    if isinstance(datalist,np.ndarray):
        d.set("frame 1")
        d.set_np2arr(datalist,dtype=float)
        d.set(scale)
        d.set("zoom to fit")
        if regions is not None:
            d.set("regions load %s" % regions)
        return target

    # If datalist is a single HDU
    '''
    if isinstance(datalist,pyfits.HDUList):
        d.set("frame 1")
        d.set_pyfits(datalist)
        d.set("zscale")
        d.set("zoom to fit")
        if regions is not None:
            d.set("regions load %s" % regions)
        return target
    '''
    
    for idx,f in enumerate(datalist):

        # Check if regular expression
        if '*' in f:
            print f
            datalist = glob.glob(f)
            return load_filelist(d,datalist,regions=None,loadregions=0,scale='zscale')
        
        d.set("frame %i" % (idx+1))
        if isinstance(f,str):
            d.set("file %s" % f)
       # elif isinstance(f,pyfits.HDUList):
       #     d.set_pyfits(f)
        elif isinstance(f,np.ndarray):
            d.set_np2arr(f,dtype=float)
        else:
            try:
                d.set_np2arr(f.data,dtype=float)
            except AttributeError:
                d.set_np2arr(f[0].data,dtype=float)
        if regions is not None:
            if loadregions == 0 or loadregions-1 == idx:
                d.set("regions load %s" % regions)
        d.set(scale)
        d.set("zoom to fit")

    d.set("match colorbar")
    d.set("match scale")
    d.set("frame 1")

    return d


def main():

    parser = argparse.ArgumentParser(description='pyds9 wrapper to display list of files or numpy arrays.')

    parser.add_argument('datalist',nargs='*',help='List of filenames or numpy arrays.')

    parser.add_argument('-r','--regions',help='Region file to display in frames.')
    parser.add_argument('--loadregions',type=int, help='Specify in which frame to load regions. 0 (default) is all frames.')

    parser.add_argument('-t','--target',help='Specify in which ds9 window to display frames.  Default is most recent or, if none currently open, new.')

    args = parser.parse_args()

    pydisplay(args.datalist, regions=args.regions, loadregions=args.loadregions,target=args.target)

    return 0


if __name__ == '__main__':
    sys.exit(main())
    
