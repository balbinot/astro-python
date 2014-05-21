#! /usr/bin/env python
import sys
import argparse
import glob

# Try importing required modules. On fail, exit with install instructions.
failed_modules = []
try:
    import numpy as np
except:
    failed_modules.append(('numpy','http://www.scipy.org/install.html'))

try:
    from scipy.ndimage.interpolation import rotate
except:
    failed_modules.append(('scipy','http://www.scipy.org/install.html'))

try:
    import ds9
except:
    failed_modules.append(('pyds9','http://hea-www.harvard.edu/RD/ds9/pyds9/'))

try:
    import pyfits
except:
    failed_modules.append(('pyfits','http://www.stsci.edu/institute/software_hardware/pyfits/Download'))

if failed_modules:
    for module,site in failed_modules:
        print "'%s' module required. Visit %s for install instructions" % (module,site)
    sys.exit('Dependencies required')



def display(datalist,regions=None):
    d = ds9.ds9()

    for idx,f in enumerate(datalist):
        d.set("frame %i" % (idx+1))
        if isinstance(f,str):
            d.set("file %s" % f)
        elif isinstance(f,np.ndarray):
            d.set_np2arr(f,dtype=float)
            
        if regions is not None:
            for region in regions:
                d.set("regions load %s" % region)

        d.set("zscale")
        d.set("zoom to fit")

    d.set("match colorbar")
    d.set("match scale")
    d.set("frame 1")
    return


def get_filenum(filename,prefix,ext):
        return int(filename[len(prefix):-len(ext)])

def main():
    parser = argparse.ArgumentParser(description='Convenience function for sending files to ds9.')
    
    parser.add_argument('filelist',metavar='files',nargs='*',help='List of files to display. Ex: "a0001.fits a0002.fits"')
    
    parser.add_argument('-s',dest='sequence',nargs=4,metavar=('prefix','first','last','ext'),help='Sequence of files to display. Ex: "a000 4 6 .fits" will display files a0004.fits, a0005.fits, and a0006.fits.')

    parser.add_argument('-r',dest='regions',nargs='+',metavar='regions',help='Region file(s). The region(s) will be loaded into all frames.')

    parser.add_argument('--angle',nargs=1,type=float,help='Specify angle in degrees to rotate image(s) before display. Note: rotated images are reshaped, and outside boundaries are filled with zeros.')

    parser.add_argument('--transpose',action='store_true',help='Specify if image(s) should be transposed upon display. Note: if --angle is specified, image will be rotated before transposition.')
    
    args = parser.parse_args()

    datalist = []
    if args.sequence is None:
        datalist = args.filelist

    else:
        # Find all files in sequence
        prefix = args.sequence[0]
        first_frame = int(args.sequence[1])
        last_frame = int(args.sequence[2])
        ext = args.sequence[3]
        
        for f in glob.glob(''.join([prefix,'*',ext])):
            try:
                filenum = get_filenum(f,prefix,ext)
            except:
                continue

            if (filenum <= last_frame) and (filenum >= first_frame):
                datalist.append(f)

        datalist = sorted(datalist, key=lambda x: get_filenum(x,prefix,ext))
        datalist = list(args.filelist) + datalist
    print datalist
    
    if (args.angle is not None) or (args.transpose) :
        # Open each file and rotate/transpose
        for data in datalist:
            data = pyfits.open(data)
            if args.angle is not None:
                data[0].data = rotate(data[0].data,np.float(args.angle[0]))
                data[0].header['ANGLE'] = ('%f' % args.angle[0],'Angle of rotation applied')
            if args.transpose:
                data[0].data = data[0].data.transpose()
        
    display(datalist,regions=args.regions)

    return
    
    
if __name__ == '__main__':
    sys.exit(main())
    
    
    
