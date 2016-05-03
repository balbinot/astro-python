#! /usr/bin/env python
from pyds9 import DS9
import argparse
import numpy as np
from photutils.morphology import centroid_2dg,centroid_com,cutout_footprint,data_properties
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.wcs.utils import proj_plane_pixel_scales
from pyraf import iraf
from astropy.io import fits
import tempfile
import os

class iDS9(DS9):
    data = None
    header = None
    wcs = None
    radius = None
    tempdatafile = None
    tempphotfile = None
    tempcoofile = None

    def __init__(self,filename=None,target='DS9:*',start=True,wait=10,verify=True,radius=15):
        super(self.__class__,self).__init__(target,start,wait,verify)
        self._setup_parser()
        self.radius = radius
        
        if filename:
            self.set(' '.join(['file',filename]))
            
        self._instantiate()

    def _instantiate(self,hdul=None):
        if hdul is None:
            try:
                hdul = self.get_pyfits()
            except (AttributeError,IOError):
                return

        self.data = hdul[0].data
        self.header = hdul[0].header
        self.wcs = WCS(self.header)

        self.tempdatafile = self.__generate_temp_file(data=hdul,suffix='.fits',ref=self.tempdatafile)
        self.tempphotfile = self.__generate_temp_file(suffix='.dat',ref=self.tempphotfile)
        self.tempcoofile = self.__generate_temp_file(suffix='.coo',ref=self.tempcoofile)
        
        
    def _setup_parser(self):
        # set up parser
        self.parser=argparse.ArgumentParser(description='Type key commands in DS9 window',prog='DS9 imexam',add_help=False,prefix_chars=' ',usage=argparse.SUPPRESS)
        keygroup = self.parser.add_argument_group('key commands', '')
        keygroup.add_argument(' h',action='store_true',required=False,help='show help')
        keygroup.add_argument(' c',action='store_true',required=False,help='Centroid')
        keygroup.add_argument(' r',action='store_true',required=False,help='Radial profile')
        keygroup.add_argument(' q',action='store_true',required=False,help='Quit')

    def __generate_temp_file(self,data=None,suffix='',ref=None):
        if ref is not None:
            os.remove(ref)
            
        with tempfile.NamedTemporaryFile(prefix='ipyds9_tmp',suffix=suffix,delete=False) as f:
            if data:
                data.writeto(f)
        return f

    def __clear_temp_files(self):
        if self.tempdatafile is not None:
            print self.tempdatafile.name, 'removed'
            os.remove(self.tempdatafile.name)

        if self.tempphotfile is not None:
            print self.tempphotfile.name, 'removed'
            os.remove(self.tempphotfile.name)

        if self.tempcoofile is not None:
            print self.tempcoofile.name, 'removed'
            os.remove(self.tempcoofile.name)



    def __del__(self):
        self.__clear_temp_files()
    
    def interact(self):
        if self.tempdatafile is None:
            self._instantiate()
            if self.tempdatafile is None:
                print 'No data in ds9'
                return
            
        self.parser.print_help()
        loop = True
                
        while loop:
            query = self.get('imexam key coordinate')
            key,x,y = query.split()

            if key not in ['h','r','c','q']:
                print query
                continue

            key = [' ' + key]


            args = self.parser.parse_args(key)

            if args.h:
                self.parser.print_help()

            if args.c:
                with open(self.tempcoofile.name,'w') as fcoo:
                    fcoo.write(' '.join([x,y]))
                #iraf.apphot.center(self.tempdatafile.name,coords=self.tempcoofile.name,output=self.tempphotfile.name,interactive=False,mode='h')
                centroid_iraf(self.tempdatafile.name,coords=self.tempcoofile.name,output=self.tempphotfile.name,radius=self.radius)
                tbl = Table.read(self.tempphotfile.name,format='ascii.daophot')
                tbl.pprint()


                #tbl = centroid(self.data,(x,y),box_size=self.radius,wcs=self.wcs)
                #tbl.pprint()

            if args.r:
                tbl = centroid(self.data,(x,y),box_size=self.radius,wcs=self.wcs)
                coord = (tbl[0]['xcen'],tbl[0]['ycen'])
                cutout,_,_,_ = cutout_footprint(self.data,coord,box_size=self.radius)
                xr,yr = radial_profile_allpix(cutout)

                '''
                if self.wcs:
                    cd1,cd2 = proj_plane_pixel_scales(self.wcs)
                    cu1,cu2 = self.wcs.wcs.cunit
                    cd1 = u.Quantity(cd1,unit=cu1)
                    cd2 = u.Quantity(cd2,unit=cu2)
                    rscale = np.sqrt(cd1**2+cd2**2).to(u.arcsec)
                    xr *= rscale
                #'''

                plt.scatter(xr,yr)
                plt.gca().yaxis.get_major_formatter().set_powerlimits((0, 1))
                #plt.gca().relim()
                #y_lims = plt.gca().get_ylim()
                #if self.wcs:
                #    plt.xlabel('Radius [arcsec]')
                #else:
                #    plt.xlabel('Radius [pix]')

                #plt.xlim([0,xr.max()])
                #plt.gca().set_ylim(y_lims)
                plt.show()




            if args.q:
                plt.close()
                loop = False

        return


def centroid_iraf(image,coords,output,radius):
    iraf.apphot.center.setParam('image', image)
    iraf.apphot.center.setParam('coords', coords)
    iraf.apphot.center.setParam('output', output)
    iraf.apphot.center.setParam('interactive', 'no')
    iraf.apphot.centerpars.setParam('calgorithm','centroid')
    iraf.apphot.centerpars.setParam('cbox', radius)
    iraf.apphot.centerpars.saveParList()
    iraf.apphot.center(image,coords=coords,output=output,mode='h')
    return output


def pix2wcs(coord,wcs):
    ncoord = wcs.all_pix2world([coord],1,ra_dec_order=True)
    ncoord = SkyCoord(ncoord,unit=(u.deg,u.deg),frame='icrs') #ra,dec
    ra,dec = ncoord.to_string('hmsdms',sep=':')[0].split()

    return ra,dec


def get_mask(data,coord,box_size=10):
    _,_,_,slices = cutout_footprint(data,coord,box_size=box_size)
    mask = np.ones_like(data,dtype=bool)
    mask[slices] = False
    return mask

def centroid(data,coord,box_size=10,wcs=None):
    coord = [np.float(coord[0]),np.float(coord[1])]
    ## make cutout mask
    mask = get_mask(data,coord,box_size)
    
    # first, try gaussian
    xc,yc = centroid_2dg(data,mask=mask)
    if np.isnan((xc,yc)).any():
        print 'Failed to centroid'
        xc,yc = coord
        
    if wcs:
        ra,dec = pix2wcs([xc,yc],wcs)
        tbl = Table([[xc],[yc],[ra],[dec]],names=['xcen','ycen','ra','dec'])
    else:
        tbl = Table([[xc],[yc]],names=['xcen','ycen'])

    return tbl



def radial_profile(data,center=None,binsize=1,interpnan=True):
    y, x = np.indices(data.shape)
    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    r = np.hypot(x - center[0], y - center[1])
    nbins = int(np.round(r.max() / binsize) + 1)
    maxbins = nbins * binsize
    bins = np.linspace(0,maxbins,nbins+1)
    bin_centers = (bins[1:] + bins[:-1])/2.0

    radial_prof = np.histogram(r,bins, weights = data)[0] / np.histogram(r,bins, weights=np.ones(data.shape))[0]

    if interpnan:
        radial_prof = np.interp(bin_centers,bin_centers[radial_prof==radial_prof],radial_prof[radial_prof==radial_prof])    

    return (bin_centers, radial_prof)


def radial_profile_allpix(data, center=None):
    y, x = np.indices(data.shape)
    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    r = np.hypot(x - center[0], y - center[1])
    yarr = data.flatten()
    xarr = r.flatten()
    
    return xarr,yarr


#for testing
def main():
    parser = argparse.ArgumentParser(description='Open ds9 window and interact with data')
    parser.add_argument('file',nargs='?',help='Open image in ds9')
    args = parser.parse_args()
    
    d = iDS9(args.file)
    d.interact()

    
if __name__ == '__main__':
    main()
