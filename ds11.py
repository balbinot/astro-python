#! /usr/bin/env python
#version 12/19/2016
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table,hstack,Column
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.visualization.mpl_normalize import ImageNormalize
from ds9norm import DS9Normalize
from astropy.visualization import SqrtStretch, LinearStretch
#from photutils.morphology import centroid_com,centroid_1dg,centroid_2dg,cutout_footprint,data_properties
from photutils import properties_table, CircularAperture, CircularAnnulus, aperture_photometry,SkyCircularAperture, SkyAperture
from sys import exit
import numpy as np
from zscale import zscale
from pyds9 import DS9
#from matplotlib.backend_bases import NavigationToolbar2, Event

#def new_home(self,*args,**kwargs):
#    s = 'home_event'
#    event = Event(s,self)
#    self.canvas.callbacks.process(s,event)
#    home(self,*args,**kwargs)
#NavigationToolbar2.home = new_home

class ds11(object):
    filename = None
    data = None
    header = None
    wcs = None
    fig = None
    ax = None
    norm = None
    normed=None
    apertures = []
    skyapertures = []


    def __init__(self, data, header=None, wcs=None):
        if isinstance(data,str):
            # data is filename
            self.filename = data
            self.data, self.header = fits.getdata(data,header=True)

            if header is not None:
                self.header = header            
        else:
            # data is array
            self.data = data
            self.header = header

        if wcs:
            self.wcs = wcs
        else:
            self.wcs = self._parse_wcs(self.header)

        if wcs == -1:
            self.wcs = None
        if self.wcs:
            self.fig,self.ax = self._setup_figure()

    def _parse_wcs(self, header):
        try:
            wcs = WCS(header,naxis=2)
        except:
            print 'Failed to parse wcs'
            return None
        return wcs

    def convert_pix(self,coord):
        ncoord = self.wcs.all_pix2world([coord],1,ra_dec_order=True)
        ncoord = SkyCoord(ncoord,unit=(u.deg,u.deg),frame='icrs') #ra,dec
        ra,dec = ncoord.to_string('hmsdms',sep=':')[0].split()

        return ra,dec
    

    def add_aperture(self,aperture,name):
        aperture.name = name
        self.apertures.append(aperture)
        return

    def add_skyaperture(self,aperture,name):
        aperture.name = name + '_sky'
        self.skyapertures.append(aperture)
        return

    @staticmethod
    def make_sky_circular_aperture(coord,rad):
        return SkyCircularAperture(coord,rad)

    def make_regions(self):
        regions = []
        for ap in self.apertures:
            if isinstance(ap.positions,SkyCoord):
                x,y = ap.positions.to_pixel(self.wcs)
                r = ap.r.value
            else:
                x,y = ap.positions[0]
                r = ap.r
            name = ap.name
            circle = 'image;circle(%f,%f,%f) # color=green text={%s} source' % (x,y,r,name)
            regions.append(circle)
        for ap in self.skyapertures:
            x,y = ap.positions[0]
            r_in,r_out = ap.r_in,ap.r_out
            name = ap.name
            annulus = 'image;annulus(%f,%f,%f,%f) # color=white text={%s} dash=1 background' % (x,y,r_in,r_out,name)
            regions.append(annulus)

        return regions

    
    def save_regions(self,regions,filename):
        with open(filename,'w') as f:
            for region in regions:
                f.write(region+'\n')
        
        
        
    def _setup_figure(self):
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=self.wcs)
        lon = ax.coords['ra']
        lat = ax.coords['dec']
        lon.set_axislabel('RA')
        lat.set_axislabel('DEC')
        lon.set_major_formatter('hh:mm:ss.s')
        lat.set_major_formatter('dd:mm:ss.s')
        ax.coords.grid(color='white')
        plt.style.context('seaborn-talk')

        # disable default key handler
        fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
        # add key handler
        keycid = fig.canvas.mpl_connect('key_press_event',self._onkey)

        # add rtclick handlers
        fig.canvas.mpl_connect('button_press_event',self._onrtclick)
        fig.canvas.mpl_connect('button_release_event',self._onrtlift)
        fig.canvas.mpl_connect('motion_notify_event',self._onrtmove)
        

        return fig,ax

    def lnorm(self,clip_lo=0,clip_hi=100,bias=0.5,contrast=1.0):
        return DS9Normalize(stretch='linear',clip_lo=clip_lo,clip_hi=clip_hi,
                            bias=bias,contrast=contrast)


    def _onkey(self, event):
        if event.key == 'q':
            plt.close()
            exit()

        '''
        if event.key == 'z':
            if self.norm and self.normed == 'z':
                self.norm = None
                self.normed = None
            else:
                vmin,vmax = zscale(self.data)
                if np.isnan(vmin):
                    vmin = 0
                if np.isnan(vmax):
                    vmax = np.nanmedian(self.data)

                self.norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())
                self.normed = 'z'
                
            self.datahandle.remove()
            self.datahandle = self.ax.imshow(self.data,cmap=plt.cm.gist_heat,origin='lower',norm=self.norm)
            self.fig.canvas.draw()
        '''
            
        return

    def _onrtclick(self,event):
        if event.button != 3: return
        if event.inaxes != self.ax: return

        self.rtxi, self.rtyi = event.xdata, event.ydata #data coordinates
        if self.fig.canvas.manager.toolbar._active == 'ZOOM':
            self.fig.canvas.manager.toolbar.zoom()
        return

    def _onrtmove(self,event):
        if event.button != 3: return
        if event.inaxes != self.ax:
            self.rtxi = None
            self.rtyi = None
            return
                               
        if self.rtxi is None:
            return

        rtxo,rtyo = event.xdata,event.ydata
        dx = rtxo-self.rtxi
        dy = rtyo-self.rtyi

        self.norm.bias += (dx/1000.)
        if self.norm.bias < 0: self.norm.bias = 0
        if self.norm.bias > 1: self.norm.bias = 1


        self.norm.contrast += (dy/200.)
        if self.norm.contrast < -2: self.norm.contrast = -2
        if self.norm.contrast > 2: self.norm.contrast = 2

        
        #self.datahandle.remove()
        self.datahandle.set_norm(self.norm)
        # self.ax.imshow(self.data,cmap=plt.cm.gist_heat,origin='lower',norm=self.norm)
        self.fig.canvas.draw()


    def _onrtlift(self,event):
        if event.button != 3: return
        self.rtxi = None
        self.rtyi = None
        return



    def show(self):
        # initialize normalization
        '''
        vmin,vmax = np.nanmin(self.data),np.nanmax(self.data)
        self.norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
        self.normed = 'm' # minmax
        '''
        
        self.norm = self.lnorm()
        
        self.datahandle = self.ax.imshow(self.data,cmap=plt.cm.bone,origin='lower',norm=self.norm)  #gist_heat

        for ap in self.apertures:
            ap.plot(ax=self.ax,color='b',lw=2)
        for ap in self.skyapertures:
            ap.plot(ax=self.ax,color='w',lw=1.5)

        
        plt.show()

    def _get_ds9(self,d=None):
        if not d:
            try:
                d = DS9('astro')
            except:
                d = DS9()
            return d

        if isinstance(d,DS9):
            return d

        if isinstance(d,str):
            try:
                d = DS9(d)
            except:
                d = DS9()
            return d



    def show_ds9(self,regions=None,d=None):
        d = self._get_ds9(d)
        
        if self.filename:
            d.set("file %s" % filename)
        elif (self.data is not None) and (self.header is not None):
            # make hdulist
            hdu = fits.PrimaryHDU(self.data,self.header)
            hdulist = fits.HDUList([hdu])
            d.set_pyfits(hdulist)
        else:
            d.set_np2arr(self.data)

        if regions:
            if regions == True:
                regions = self.make_regions()
            d.set('regions', '\n'.join(regions))

            first_ap = self.apertures[0]
            if isinstance(first_ap.positions,SkyCoord):
                x,y = first_ap.positions.to_pixel(self.wcs)
            else:
                x,y = first_ap.positions[0]
            d.set('pan to %i %i image' % (int(x),int(y)))
        
            


    @staticmethod
    def is_pix_coord(coord):
        if isinstance(coord,SkyCoord):
            return False
        try:
            x = np.float(coord[0])
            y = np.float(coord[1])
        except ValueError:
            return False
        return True

    @staticmethod
    def is_sky_coord(coord):
        if isinstance(coord,SkyCoord):
            return True
        try:
            coo = SkyCoord(' '.join(coord),unit=(u.hourangle,u.deg),frame='icrs')
        except ValueError:
            return False
        return True

    @staticmethod
    def pix2wcs(coord,wcs):
        ncoord = wcs.all_pix2world([coord],1,ra_dec_order=True)
        ncoord = SkyCoord(ncoord,unit=(u.deg,u.deg),frame='icrs') #ra,dec
        ra,dec = ncoord.to_string('hmsdms',sep=':')[0].split()

        return ra,dec

    def make_cutout(self,coord,box_size,as_ds11=False):
        '''coord is SkyCoord, box_size is in arcsec'''
        if hasattr(box_size, 'unit'):
            xw = box_size
            yw = box_size

        else:
            xw = box_size*u.arcsec
            yw = box_size*u.arcsec
            
        head = self.header.copy()
        cd1 = head['CDELT1'] if head.get('CDELT1') else head.get('CD1_1')
        cd2 = head['CDELT2'] if head.get('CDELT2') else head.get('CD2_2')
        if cd1 is None or cd2 is None:
            raise Exception("Missing CD or CDELT keywords in header %s" % self.filelist[im])
        cd1 = np.abs(cd1) * u.deg
        cd2 = np.abs(cd2) * u.deg
        wcs = WCS(head)
        xx,yy = wcs.wcs_world2pix(coord.ra,coord.dec,0)

        #print xw,yw,xx,yy,cd1,cd2
        xmin,xmax = np.max([0,xx-xw/cd1]),np.min([head['NAXIS1'],xx+xw/cd1])
        ymin,ymax = np.max([0,yy-yw/cd2]),np.min([head['NAXIS2'],yy+yw/cd2])

        head['CRPIX1'] -= xmin
        head['CRPIX2'] -= ymin
        head['NAXIS1'] = int(xmax-xmin)
        head['NAXIS2'] = int(ymax-ymin)

        #print xmin,xmax,ymin,ymax
        data = self.data[ymin:ymax,xmin:xmax]

        if as_ds11:
            return ds11(data,head)
        else:
            return data,header



def push_frames(filelist,regions=None,pan=None):
    try:
        d = DS9('astro')
    except:
        d = DS9()

    for idx,filename in enumerate(filelist):
        frame = idx+1
        d.set('frame %i' % frame)
        d.set('file %s' % filename)
        if regions:
            d.set('regions', '\n'.join(regions))
        d.set('scale mode 99.5')
        if pan:
            d.set('pan to %i %i image' % (int(pan[0]),int(pan[1])))


    d.set('frame 1')
    d.set('lock frame image')
