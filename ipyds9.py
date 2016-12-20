#! /usr/bin/env python
###
## launched from 'source activate irafenv'

from pyds9 import DS9
import argparse
import numpy as np
from photutils.morphology import centroid_2dg,centroid_com,cutout_footprint,data_properties
from astropy.table import Table,Column,join
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.wcs.utils import proj_plane_pixel_scales
from pyraf import iraf
from astropy.io import fits
import tempfile
from collections import OrderedDict
import os
from astropy.modeling import models,fitting
import warnings
import re
from utils import RadProf

C_COLMAP = OrderedDict([('XINIT','xi'),('YINIT','yi'),('XCENTER','xc'),('YCENTER','yc'),('XSHIFT','xs'),('YSHIFT','ys'),('CERROR','err')])
P_COLMAP = OrderedDict([('Column','xi'),('Line','yi'),('PA','pa'),('Ellip','ellip'),('FWHM','hwhm'),('GFWHM','ghwhm'),('MFWHM','mhwhm')])
A_COLMAP = OrderedDict([('XCENTER','xc'),('YCENTER','yc'),('AREA','area'),('MSKY','msky'),('FLUX','flux'),('PERROR','err')])
R_COLMAP = ['prad', 'frad', 'flux', 'bgd', 'mpeak', 'ellip', 'pa', 'mbeta', '2dfwhm', 'mfwhm', 'fwhm']

#GKIDECODE expressions
fitRe = re.compile('color=1\s*\n\s*(?:polyline.*\s*([\d\.\d*\s]*))*flush')
pairRe = re.compile('(\d\.\d*\s\d\.\d*)')
dataRe = re.compile('color=1\n(?:(?:polyline.*\n)\s(.*\n))+flush')
extractRe = re.compile('xavg=(\d\.\d\d).*yavg=(\d\.\d\d)')
ndcRe = re.compile('set_wcs.*\n(.*)')




class iDS9(DS9):
    data = None
    header = None
    wcs = None
    radius = None
    tempdatafile = None
    tempphotfile = None
    tempcoofile = None
    tempradfile = None

    def __init__(self,filename=None,target='DS9:*',start=True,wait=10,verify=True,radius=15,arad=3,irad=5,orad=7):
        super(self.__class__,self).__init__(target,start,wait,verify)
        self._setup_parser()
        self.radius = radius
        self.arad = 3
        self.irad = 5
        self.orad = 7
        
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
        if self.wcs:
            cd1,cd2 = proj_plane_pixel_scales(self.wcs)
            cu1,cu2 = self.wcs.wcs.cunit
            cd1 = u.Quantity(cd1,unit=cu1)
            cd2 = u.Quantity(cd2,unit=cu2)
            try:
                self.pxscale = np.sqrt(cd1*cd2).to(u.arcsec)/u.pix
            except Exception as e:
                print 'Error parsing wcs'
                print e,'\n'
                self.wcs = None
                self.pxscale = None
        else:
            self.pxscale = None
            


        self.tempdatafile = self.__generate_temp_file(data=hdul,suffix='.fits',ref=self.tempdatafile)
        self.tempphotfile = self.__generate_temp_file(suffix='.dat',ref=self.tempphotfile)
        self.tempcoofile = self.__generate_temp_file(suffix='.coo',ref=self.tempcoofile)
        self.tempradfile = self.__generate_temp_file(suffix='.prf',ref=self.tempradfile)
        
        
    def _setup_parser(self):
        # set up parser
        self.parser=argparse.ArgumentParser(description='Type key commands in DS9 window',prog='DS9 imexam',add_help=False,prefix_chars=' ',usage=argparse.SUPPRESS)
        keygroup = self.parser.add_argument_group('key commands', '')
        keygroup.add_argument(' h',action='store_true',required=False,help='show help')
        keygroup.add_argument(' c',action='store_true',required=False,help='Centroid')
        keygroup.add_argument(' a',action='store_true',required=False,help='Quick photometry')
        keygroup.add_argument(' r',action='store_true',required=False,help='Radial profile')
        keygroup.add_argument(' p',action='store_true',required=False,help='PSF measure')
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

        if self.tempradfile is not None:
            print self.tempradfile.name, 'removed'
            os.remove(self.tempradfile.name)


    def _write_coord_file(self,x,y):
        with open(self.tempcoofile.name,'w') as fcoo:
            fcoo.write(' '.join([x,y]))
        return self.tempcoofile.name

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
        last_command = None

        #radprofiler
        r = RadProf(self.tempdatafile.name)
        
        while loop:
            query = self.get('imexam key coordinate')
            key,x,y = query.split()

            key = [' ' + key]


            args = self.parser.parse_args(key)

            if args.h:
                self.parser.print_help()

            if args.c:
                #with open(self.tempcoofile.name,'w') as fcoo:
                #    fcoo.write(' '.join([x,y]))
                self._write_coord_file(x,y)
                tbl = centroid_iraf(self.tempdatafile.name,coords=self.tempcoofile.name,radius=self.radius,wcs=self.wcs)

                # only print col names if fresh
                if last_command == 'c':
                    tbl.pprint(show_name=False,show_unit=False)
                else:
                    tbl.pprint(show_unit=False)


            if args.a:
                self._write_coord_file(x,y)
                tbl = qphot_iraf(self.tempdatafile.name,coords=self.tempcoofile.name,radius=self.radius,arad=self.arad,irad=self.irad,orad=self.orad)
                
            if args.r:
                prof,tbl = r.profile((x,y))
                if last_command == 'r':
                    tbl.pprint(show_name=False,show_unit=False)
                else:
                    tbl.pprint()

                prof[0].plot(self.pxscale)

                
                
                ###centroid first
                ###self._write_coord_file(x,y)
                ###ctbl = centroid_iraf(self.tempdatafile.name,coords=self.tempcoofile.name,radius=self.radius)
                ###xc,yc = ctbl['xc'][0], ctbl['yc'][0]
                #tbl = radialForTerry(self.data,(xc,yc),self.radius)
                #tbl = radprof(self.data,(xc,yc),radius=self.radius,pxscale=self.pxscale)
                #####tbl = radial_profile_allpix(self.data,(xc,yc),self.radius,pxscale=self.pxscale)
                #get psfmeasure for moffat fit
                ###ptbl = psfmeasure(self.tempdatafile.name,coords=self.tempcoofile.name,radius=self.radius,size='all',pxscale=self.pxscale)
                ###ptbl.pprint()

                #ptbl = radprof_iraf(self.tempdatafile.name,coords=self.tempcoofile.name,radius=self.radius,verbose=True,test=True,plotfile=self.tempradfile.name,pxscale=self.pxscale)
                ####tbl = radial_profile_allpix2(self.data,(xc,yc),self.radius,pxscale=self.pxscale)
                #data = '\n'.join(['%f %f' % (x,y) for x,y in zip(tbl['radius'],tbl['counts'])])
                #self.set('plot new "Radial Profile" scatter')

                #moffat
                ###mof = fit_moffat_1d(ptbl['radius'],ptbl['intensity'])
                #mof = fit_moffat_1d(ptbl['radius'],ptbl['intensity'])
                #mhwhm = mof.gamma*np.sqrt(2.**(1./mof.alpha)-1)
                #print mhwhm
                #print mof.x_0

                ###UNCOMMENT THIS NEXT BIT TO PLOT
                #plot_rad(tbl,ptbl,mof,pxscale=self.pxscale)


                #tbl = radprof_iraf(self.tempdatafile.name,coords=self.tempcoofile.name,radius=self.radius,plotfile=self.tempradfile.name,verbose=True)
                
                #tbl = Table.read(tbl,format='ascii.daophot')

                #tbl = Table.read(self.tempphotfile.name,format='ascii.daophot')
                #tbl.pprint()
                '''
                tbl = centroid(self.data,(x,y),box_size=self.radius,wcs=self.wcs)
                coord = (tbl[0]['xcen'],tbl[0]['ycen'])
                cutout,_,_,_ = cutout_footprint(self.data,coord,box_size=self.radius)
                xr,yr = radial_profile_allpix(cutout)
                #'''

                '''
                if self.wcs:
                    cd1,cd2 = proj_plane_pixel_scales(self.wcs)
                    cu1,cu2 = self.wcs.wcs.cunit
                    cd1 = u.Quantity(cd1,unit=cu1)
                    cd2 = u.Quantity(cd2,unit=cu2)
                    rscale = np.sqrt(cd1**2+cd2**2).to(u.arcsec)
                    xr *= rscale
                #'''

                #plt.scatter(xr,yr)
                #plt.gca().yaxis.get_major_formatter().set_powerlimits((0, 1))
                #plt.gca().relim()
                #y_lims = plt.gca().get_ylim()
                #if self.wcs:
                #    plt.xlabel('Radius [arcsec]')
                #else:
                #    plt.xlabel('Radius [pix]')

                #plt.xlim([0,xr.max()])
                #plt.gca().set_ylim(y_lims)
                #plt.show()
                

            if args.p:
                self._write_coord_file(x,y)
                tbl = psfmeasure(self.tempdatafile.name,coords=self.tempcoofile.name,radius=self.radius,pxscale=self.pxscale)
                
                # only print col names if fresh
                if last_command == 'p':
                    tbl.pprint(show_name=False,show_unit=False)
                else:
                    tbl.pprint()

            if args.q:
                plt.close()
                loop = False

            last_command = key[0].strip()
        return


def centroid_iraf(image,coords,radius,output='STDOUT',verbose=False,wcs=None):
    iraf.apphot.centerpars.setParam('maxshift',radius)
    results = iraf.apphot.center(image,coords=coords,output='STDOUT',verbose=verbose,verify=False,interactive=False,Stdout=1)
    tbl = Table.read(results,format='ascii.daophot')
    tbl.keep_columns(C_COLMAP.keys())
    for k,v in C_COLMAP.iteritems():
        tbl.rename_column(k,v)

    if wcs:
        pixcoord = (tbl['xc'][0],tbl['yc'][0])
        wcscoord = wcs.all_pix2world([pixcoord],1,ra_dec_order=True)
        ra,dec = SkyCoord(wcscoord,unit=(u.deg,u.deg)).to_string('hmsdms',sep=':',precision=2)[0].split()
        tbl.replace_column('xs',Column([ra],name='ra'))
        tbl.replace_column('ys',Column([dec],name='dec'))
        tbl.rename_column('xs','ra')
        tbl.rename_column('ys','dec')
    return tbl


def qphot_iraf(image,coords,radius=15,arad=3,irad=5,orad=7,verbose=False):
    results = iraf.apphot.qphot(image,coords=coords,cbox=radius,annulus=irad,dannulus=orad,apertures=arad,interactive=False,output='STDOUT',Stdout=1)
    tbl = Table.read(results,format='ascii.daophot')
    print tbl.colnames
    tbl.keep_columns(A_COLMAP.keys())
    for k,v in A_COLMAP.iteritems():
        tbl.rename_column(k,v)
    tbl.pprint()

def radprof(data,coords,radius,output='STDOUT',pxscale=None):
    xc,yc = int(coords[0]), int(coords[1])
    # just grab the data box we want from the image
    data_chunk = data[yc-radius:yc+radius,
                      xc-radius:xc+radius]

    y, x = np.indices((data_chunk.shape))
    r = np.sqrt((x - radius)**2 + (y - radius)**2)
    r = r.astype(np.int)
    # add up the flux in integer bins
    tbin = np.bincount(r.ravel(), data_chunk.ravel())
    nr = np.arange(len(tbin))
    print nr
    print tbin
    plt.scatter(nr,tbin)
    plt.show()


def psfmeasure(image,coords,radius,output='STDOUT',size='FWHM',pxscale=None):
    iraf.noao(_doprint=0)
    iraf.obsutil(_doprint=0)
    if size == 'all':
        sname = ['FWHM','GFWHM','MFWHM']
        tbls = []
        for s in sname:
            stbl = psfmeasure(image,coords,radius,output,s,pxscale)
            tbls.append(stbl)
        tbl = join(tbls[0],tbls[1])
        tbl = join(tbl,tbls[2])
        return tbl
        
        
    results = iraf.obsutil.psfmeasure(image,coords='mark1',size=size,imagecur=coords,display='no',radius=radius,logfile='',Stdout=1,StdoutG="dev$null",scale=1)
    #print results
    results = results[2:4]
    results[0] = results[0].strip()
    results = {k:v for k,v in zip(results[0].split(),results[1].split())}
    for k in results.keys():
        try:
            results[k] = np.float(results[k])
        except ValueError:
            continue
    tbl = Table([results])
    tbl = tbl[[k for k in P_COLMAP.keys() if k in tbl.colnames]]
    for k,v in P_COLMAP.iteritems():
        if k in tbl.colnames:
            tbl.rename_column(k,v)

    fwhmcol = tbl.colnames[-1]
    tbl[fwhmcol] /= 2.0  # make hwhm
    if pxscale:
        tbl.add_column(Column(tbl[fwhmcol]*u.pix*pxscale,name='%s_s'%fwhmcol,unit=u.arcsec,format='%.3f'))
        tbl[fwhmcol].unit=u.pix
    
    return tbl


def extract_gkidata(text):
    fitdata = fitRe.findall(text)
    fitdata = pairRe.findall(fitdata[1]) # get pair vals
    fitdata = [x.split() for x in fitdata]
    fitdata = [(np.float(x[0]),np.float(x[1])) for x in fitdata]
    xf,yf = zip(*fitdata)

    data = dataRe.search(text).group(0)
    data = extractRe.findall(data)
    data = [(np.float(x[0]),np.float(x[1])) for x in data]

    radius,pixelval = zip(*data)

    #get units
    unitline = ndcRe.findall(text)[0]
    unitline = unitline.split()
    #id x1 x2 y1 y2 x1 x2 y1 y2
    x1r,x2r = np.float(unitline[1]),np.float(unitline[2])
    y1r,y2r = np.float(unitline[3]),np.float(unitline[4])

    x1n,x2n = np.float(unitline[5]),np.float(unitline[6])
    y1n,y2n = np.float(unitline[7]),np.float(unitline[8])

    # linear fit in x
    xfit = np.polyfit((x1n,x2n),(x1r,x2r),deg=1)
    xfit = np.poly1d(xfit)
    radius = xfit(radius)
    xf = xfit(xf)
    return radius,pixelval,xf,yf

    

def radprof_iraf(image,coords,radius,step=0.1,output='STDOUT',plotfile=None,verbose=False,test=False,pxscale=None):
    iraf.apphot.centerpars.setParam('maxshift',radius)  # search radius
    if test:
        results = iraf.imexamine(image,image=image,use_display='no',logfile='',defkey='r',imagecur=coords,frame=1,mode='h',StdoutG=plotfile,Stdout=1)
        results = [np.float(x) for x in results[0].split()]
        # prad frad flux bgd mpeak ellip pa mbeta 2dfwhm mfwhm fwhm
        results = OrderedDict(zip(R_COLMAP,results))
        tbl = Table([results])
        tbl = tbl[R_COLMAP]
        tbl['flux'].unit = u.ct
        tbl['bgd'].unit = u.ct
        tbl['mpeak'].unit = u.ct
        tbl['pa'].unit = u.deg
        #tbl['2dfwhm'].unit = u.pix
        tbl['mfwhm'].unit = u.pix
        tbl['fwhm'].unit = u.pix

        #remove 2dfwhm cuz wut?
        tbl.remove_column('2dfwhm')

        if pxscale:
            #add arcsec columns
            for fwhmcol in tbl.colnames[-2:]:
                tbl.add_column(Column(tbl[fwhmcol]*pxscale,name='%s_s'%fwhmcol,unit=u.arcsec,format='%.3f'),tbl.index_column(fwhmcol)+1)

        gki = iraf.gkidecode(plotfile,Stdout=1,verbose=True)
        gki = '\n'.join(gki)
        xr,yr,xf,yf = extract_gkidata(gki)
        #plt.scatter(xr,yr)
        #plt.plot(xf,yf)
        #plt.show()

        return tbl,xr,yr,xf,yf

    results = iraf.apphot.radprof(image,radius=radius,step=step,coords=coords,output=output,verbose=verbose,verify=False,interactive=False,Stdout=1)

    #remove datalines
    data = []
    pradline = None
    for idx,line in enumerate(results):
        if '*\\' in line:
            line = line.split()[0:3]
            line = [np.float(x) for x in line]
            data.append((idx,line))
        #remove pradius line descriptors
        elif 'PRADIUS' in line:
            pradline = idx
    idx,data = zip(*data)
    #print data
    #remove from min idx onward
    ridx = np.min(idx)
    del results[ridx:]
    del results[pradline:pradline+4]
    # remove '\' from last line
    results[-1] = results[-1][0:-1]
    results = '\n'.join(results[1:])
    
    rtbl = Table.read(results,format='ascii.daophot')
    #print rtbl['PFWHM']/2.0
    ptbl = Table(rows=data,names=('radius','intensity','tintensity'))
    ptbl['intensity'] *= rtbl['INORM']
    #tbl.pprint()
    #plt.plot(tbl['radius'],tbl['intensity']*ptbl['INORM'])
    return ptbl


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



def interpol(x,y,edgeclip=False):
    yp = np.interp(x,x[y==y],y[y==y])
    if edgeclip:
        return yp[:-1]
    else:
        return yp

def radialForTerry(data,center,radius,binsize=1,interpnan=True):
    data_chunk,_,_,_ = cutout_footprint(data,center,box_size=radius)
    # get indices of data
    y, x = np.indices(data_chunk.shape)
    #new center is center of chunk
    cen = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    r = np.hypot(x - cen[0], y - cen[1])
    nbins = int(np.round(r.max() / binsize) + 1)
    maxbins = nbins * binsize
    bins = np.linspace(0,maxbins,nbins+1)
    bin_centers = (bins[1:] + bins[:-1])/2.0
    print bin_centers
    
    y = np.histogram(r,bin_centers, weights = data_chunk)[0] / np.histogram(r,bin_centers, weights=np.ones(data_chunk.shape))[0]
    
    if interpnan:
        y = interpol(bin_centers,y,edgeclip=True)
    #skip nans
    idx = np.where(y==y)
    bin_centers = np.array(bin_centers)[idx]
    y = y[idx]

    plt.scatter(bin_centers,y)
    plt.show()
    return bin_centers, y


def radial_profile_allpix2(data,center,radius,pxscale=None):
    y,x = np.indices(data.shape)
    r = np.hypot(x - center[0],y-center[1])
    yarr = data.flatten()
    xarr = r.flatten()
    idx = np.where(xarr <= radius)
    xarr = xarr[idx]
    yarr = yarr[idx]
    tbl = Table([xarr,yarr],names=('radius','counts'))
    if pxscale:
        tbl.add_column(Column(tbl['radius']*pxscale,name='radius_s',unit=u.arcsec))
        tbl['radius'].unit = u.pix
    
    return tbl
    
    

def radial_profile_allpix(data, center, radius,pxscale=None):
    data_chunk,data_mask,_,data_slices = cutout_footprint(data,center,box_size=radius)
    y, x = np.indices(data_chunk.shape)
    ### new center is center of array
    cen = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    r = np.hypot(x - cen[0], y - cen[1])
    ###r = np.hypot(x - center[0], y - center[1])
    yarr = data_chunk.flatten()
    xarr = r.flatten()
    tbl = Table([xarr,yarr],names=('radius','counts'))
    if pxscale:
        tbl.add_column(Column(tbl['radius']*pxscale,name='radius_s',unit=u.arcsec))
        tbl['radius'].unit = u.pix
    
    return tbl


def plot_rad(tbl,ptbl,mof,pxscale=None):
    plt.scatter(tbl['radius'],tbl['counts'])
    plt.xlabel('Radius [pix]')
    plt.ylabel('Counts')

    plt.plot(ptbl['radius'],mof(ptbl['radius']),'r')
    plt.plot(ptbl['radius'],ptbl['intensity'],'k')
    
    if pxscale:
        ax1 = plt.gca()
        ax2 = ax1.twiny()
        ax2Ticks = [tick*pxscale for tick in ax1.get_xticks()]
        ax2Labels = ['%.2f'%tick.value for tick in ax2Ticks]
        ax2.set_xticks(ax1.get_xticks())
        ax2.set_xbound(ax1.get_xbound())
        ax2.set_xticklabels(ax2Labels)
        ax2.set_xlabel('Radius [arcsec]')
    plt.show()

def fit_moffat_1d(x, y, gamma=2., alpha=1.):
    """Fit a 1D moffat profile to the data and return the fit."""
    # Fit model to data
    fit = fitting.LevMarLSQFitter()

    x_0 = x[0]
    
    # Moffat1D
    model = models.Moffat1D(amplitude=max(y), x_0=x_0, gamma=gamma, alpha=alpha)
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        results = fit(model, x, y)

    # previous yield amp, ycenter, xcenter, sigma, offset
    return results


def generate_temp_file(self,data=None,suffix='',ref=None):
    if ref is not None:
        os.remove(ref)
            
    with tempfile.NamedTemporaryFile(prefix='ipyds9_tmp',suffix=suffix,delete=False) as f:
        if data:
            data.writeto(f)
    return f


#for testing
def main():
    parser = argparse.ArgumentParser(description='Open ds9 window and interact with data')
    parser.add_argument('file',nargs='?',help='Open image in ds9')
    parser.add_argument('-srad',metavar='search_radius',type=float,default=15,help='Search radius for centroiding in pixels (default=15)')
    parser.add_argument('-arad',metavar='ap_radius',type=float,default=3,help='Aperture radius in pixels (default=3)')
    parser.add_argument('-irad',metavar='inner_radius',type=float,default=5,help='Inner radius of sky annulus in pixels (default=5)')
    parser.add_argument('-orad',metavar='outer_radius',type=float,default=7,help='Inner radius of sky annulus in pixels (default=7)')
                                 
    args = parser.parse_args()
    
    d = iDS9(args.file,radius=args.srad,arad=args.arad,irad=args.irad,orad=args.orad)
    d.interact()

    
if __name__ == '__main__':
    main()
