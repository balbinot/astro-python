#! /usr/bin/env python
## 01/24/2017
## NOTE!  deprecated CardList from pyfits api require modification of
##   wcs_helper.py
##     line 248:  self._pywcs = pywcs.WCS(header)
from astropy.wcs import WCS
from astropy.wcs.utils import wcs_to_celestial_frame
from astropy.io import fits
import numpy as np
import pywcsgrid2
import matplotlib.pyplot as plt
from scipy import ndimage
from astropy.coordinates import SkyCoord
import astropy.units as u
from glob import glob
import os.path


class Reproject(object):
    def __init__(self,filename,xlim=None,ylim=None,fig=None,vmin=None,vmax=None,twin=None,cmap=None,right=False,aspect='auto',gs=None,angle=0.0,interp=None,shift=None,norm=None,plot=True,reproject=None,exact=False,hide_xaxis=False,hide_yaxis=False,forcemontage=False,userepr=False):
        self.filename = filename
        self.fig = fig
        self.hdu = fits.open(self.filename,ext=0)
        self.data = self.hdu[0].data
        self.h = self.hdu[0].header

        try:
            if self.data.ndim == 3:
                self.data = self.data[0]
                self.hdu[0].data = self.data
        except AttributeError:
            # PACS is ext 1.  Also, cancer
            self.hdu[0].data,self.hdu[0].header = fits.getdata(self.filename,ext=1,header=True)
            self.h = self.hdu[0].header
            self.data = self.hdu[0].data

        # For dumb-ass PACS images
        if 'CUNIT1' in self.h and self.h['CUNIT1'] in ['Degree','Degrees']:
            #self.h = fits.getheader(self.filename,ext=1)
            self.h['CUNIT1'] = ('deg','WCS unit1')
            self.h['CUNIT2'] = ('deg','WCS unit2')

            self.w = WCS(self.h)
            self.h = self.w.to_header()
            self.h['NAXIS1'] = self.data.shape[1]
            self.h['NAXIS2'] = self.data.shape[0]
            self.h['NAXIS'] = 2
            self.hdu[0].header = self.h

        elif self.h['NAXIS'] == 3:
            self.h['NAXIS'] = 2
            del self.h['NAXIS3']
            self.w = WCS(self.h,naxis=2)
            self.hdu[0].header= self.h

        else:
            self.w = WCS(self.h)

        self.w_orig = self.w.copy()
        self.xlim = xlim
        self.ylim = ylim
        self.vmin = vmin
        self.vmax = vmax
        self.cmap = cmap
        self.twin = twin
        self.right = right
        self.aspect = aspect
        self.gs = gs
        self.angle = angle
        self.interp = interp
        self.shift = shift
        self.norm = norm
        self.hide_xaxis = hide_xaxis
        self.hide_yaxis = hide_yaxis
        self.reproject = reproject  # header of reprojection
        self.exact = exact # exact reprojection or interp

        if self.reproject:
            self._transform(userepr=userepr)
        else:
            #self._rotate()
            self._rotate_WORKING(forcemontage)

        if self.shift:
            self.data = ndimage.interpolation.shift(self.data,self.shift)
            self.hdu[0].data = self.data
        
        if plot:
            self._make_axis()

    def _transform(self,userepr=False):
        from reproject import reproject_interp,reproject_exact
        if self.hdu[0].header == self.reproject:
            return

        if userepr:
            to_glob = os.path.splitext(self.filename)[0]
            reprfile = ''.join([to_glob,'.repr.fits'])
            print 'Using repr for reproj'
            self._rotate_WORKING()

        if self.exact:
            array,_= reproject_exact(self.hdu, self.reproject,hdu_in=0)
        else:
            array,_= reproject_interp(self.hdu, self.reproject,order='biquadratic')
        self.hdu[0].data = array
        self.w = WCS(self.reproject)
        self.data = self.hdu[0].data

    def plot(self):
        if hasattr(self,'ax'):
            print 'already plotted'
            return

        self._make_axis()

    def montage(self,reprfile):
        from montage_wrapper.wrappers import reproject
        reproject(str(self.filename),reprfile,system='EQUJ',north_aligned=True)
        
        
    def _rotate_WORKING(self,forcemontage=False):
        lattyp,lngtyp = self.w.wcs.lattyp,self.w.wcs.lngtyp
        if lattyp == 'GLAT':
            # FUCKED, requires reprojection
            to_glob = os.path.splitext(self.filename)[0]
            reprfile = ''.join([to_glob,'.repr.fits'])

            if forcemontage:
                print 'Generating %s...'%reprfile
                self.montage(reprfile)

            else:
                try:
                    hdu = fits.open(reprfile)
                    print 'Using %s' % reprfile
                except IOError:
                    #must reproject
                    print '%s reprfile not found. Generating...'%reprfile
                    self.montage(reprfile)
                    hdu = fits.open(reprfile)

            self.hdu = hdu
            self.filename = reprfile
            self.data = self.hdu[0].data
            self.h = self.hdu[0].header
            self.w = WCS(self.h)
            return

        else:
            # Construct header in equatorial frame
            self.h1 = self.h.copy()

            # Shift the WCS reference pixel to center of grid
            nx, ny = np.float(self.h["naxis1"]), np.float(self.h["naxis2"])
            i0, j0 = (nx + 1.0)/2, (ny + 1.0)/2  # 1 because indexing

            (ra0,dec0), = self.w.all_pix2world([[i0, j0]], 1)

            self.h1.update(crpix1=i0, crpix2=j0, crval1=ra0, crval2=dec0)

            self.h1['CUNIT1'] = ('deg','WCS unit1')
            self.h1['CUNIT2'] = ('deg','WCS unit2')

            try:
                self.h1.update(
                    cd1_1=-np.hypot(self.h["CD1_1"], self.h["CD1_2"]),
                    cd1_2=0.0,
                    cd2_1=0.0, 
                    cd2_2=np.hypot(self.h["CD2_1"], self.h["CD2_2"]),
                    orientat=0.0)
            except KeyError:
                self.h1.update(
                    cd1_1=-np.hypot(self.h["CDELT1"], self.h["CDELT2"]),
                    cd1_2=0.0,
                    cd2_1=0.0,
                    cd2_2=np.hypot(self.h["CDELT1"], self.h["CDELT2"]),
                    orientat=0.0)

            self.h = self.h1
            self.hdu[0].header = self.h
            self.w = WCS(self.h)



    def _rotate(self):
        # Construct header in equatorial frame
        self.h1 = self.h.copy()
        
        # Shift the WCS reference pixel to center of grid
        nx, ny = np.float(self.h["naxis1"]), np.float(self.h["naxis2"])
        i0, j0 = (nx + 1.0)/2, (ny + 1.0)/2  # 1 because indexing
        
        (ra0,dec0), = self.w.all_pix2world([[i0, j0]], 1)

        # ALL OF THIS IS FOR MSX BULLSHIT
        lattyp,lngtyp = self.w.wcs.lattyp,self.w.wcs.lngtyp
        if lattyp == 'GLAT':
            coord = SkyCoord(frame='galactic',b=dec0*u.deg,l=ra0*u.deg)
            coord =  coord.transform_to('fk5')
            ra0 = coord.ra.value
            dec0 = coord.dec.value
            
            typ = self.h1['CTYPE1'].split('-')[1]
            
            if typ == 'CAR':
                '''
                #force remap
                from tempfile import NamedTemporaryFile
                import subprocess
                from os import unlink
                # run wcstools remap
                with NamedTemporaryFile(delete=True) as f:
                    subprocess.call(['remap',
                                     '-w','TAN','-v',
                                     '-o',f.name,self.filename])
                    self.h1 = fits.getheader(f.name)
                '''
                pass

        self.h1.update(crpix1=i0, crpix2=j0, crval1=ra0, crval2=dec0)
        if lattyp == 'GLAT':
            typ = self.h1['CTYPE1'].split('-')[1]
            self.h1['CTYPE1'] = '---'.join(['RA',typ])
            self.h1['CTYPE2'] = '--'.join(['DEC',typ])
        self.h1['CUNIT1'] = ('deg','WCS unit1')
        self.h1['CUNIT2'] = ('deg','WCS unit2')
        
        try:
            self.h1.update(
                cd1_1=-np.hypot(self.h["CD1_1"], self.h["CD1_2"]),
                cd1_2=0.0,
                cd2_1=0.0, 
                cd2_2=np.hypot(self.h["CD2_1"], self.h["CD2_2"]),
                orientat=0.0)
        except KeyError:
            self.h1.update(
                cd1_1=-np.hypot(self.h["CDELT1"], self.h["CDELT2"]),
                cd1_2=0.0,
                cd2_1=0.0,
                cd2_2=np.hypot(self.h["CDELT1"], self.h["CDELT2"]),
                orientat=0.0)

        self.h = self.h1
        self.hdu[0].header = self.h
        self.w = WCS(self.h)

        # MORE MSX BULLSHIT
        #temph = self.w.to_header()
        typ = temph['CTYPE1'].split('-')[-1]
        if typ == 'CAR':
            # require reprojection
            pass

            
            '''
            rotate = temph['LATPOLE']
            #force remap
            from tempfile import NamedTemporaryFile
            import subprocess
            from os import unlink
            # run wcstools remap
            with NamedTemporaryFile(delete=True) as f:
                subprocess.call(['remap','-j',str(ra0),str(dec0),
                                 #'-a','%s'%(90.-rotate),
                                 '-w','TAN','-v',
                                 '-o',f.name,self.filename])
                self.h1 = fits.getheader(f.name)
                self.w = WCS(self.h1)
                print self.w
            '''
            '''
            from astropy.nddata.utils import Cutout2D
            diag = np.hypot(nx,ny)
            cutout = Cutout2D(self.data,(i0,j0),size=diag,mode='partial')
            print cutout
            print cutout.data.shape
            '''

        #fits.writeto('msx_test_image.fits',self.data,header=self.w.to_header(),overwrite=True,output_verify='silentfix')
        #fits.writeto('msx_test_image.fits',cutout.data,header=self.w.to_header(),overwrite=True,output_verify='silentfix')
        #exit()


    def _make_axis(self):
        if self.fig is None:
            self.fig = plt.figure()

        if self.cmap is None:
            self.cmap = plt.cm.gray_r
        
        if self.gs:
            gs = self.gs
        else:
            gs = 121

        self.ax = plt.subplot(gs, projection=self.w)
        self.ax.imshow(self.data, origin='lower',vmin=self.vmin,vmax=self.vmax,cmap=self.cmap,interpolation=self.interp,norm=self.norm)

        if self.xlim:
            self.ax.set_xlim(self.xlim)
        if self.ylim:
            self.ax.set_ylim(self.ylim)

        self.ax.grid(color='black', alpha=0.2, linestyle='dashed')
        #self.ax.grid(color='white')

        try:
            self.xax = self.ax.coords['ra']
            self.yax = self.ax.coords['dec']
            self.xax.set_major_formatter('hh:mm:ss')
            self.yax.set_major_formatter('dd:mm:ss')
            self.xax.set_axislabel(r"$\alpha_{2000}$",fontsize=14)
            self.yax.set_axislabel(r"$\delta_{2000}$",fontsize=14)

        except KeyError:
            self.xax = self.ax.coords['glon']
            self.yax = self.ax.coords['glat']
            self.xax.set_major_formatter('dd:mm:ss')
            self.yax.set_major_formatter('dd:mm:ss')
            #self.xax.set_axislabel(r"$\alpha_{2000}$",fontsize=14)
            #self.yax.set_axislabel(r"$\delta_{2000}$",fontsize=14)
                    
        self.xax.set_ticks(number=4)
        self.yax.set_ticks(number=4)

        if self.right:
            self.yax.set_ticks_position('r')
            self.yax.set_ticklabel_position('r')
            self.yax.set_axislabel_position('r')

        if self.hide_xaxis:
            self.xax.set_ticks_visible(False)
            self.xax.set_ticklabel_visible(False)
            self.xax.set_axislabel('')

        if self.hide_yaxis:
            self.yax.set_ticks_visible(False)
            self.yax.set_ticklabel_visible(False)
            self.yax.set_axislabel('')


    def show_xaxis(self):
        self.xax.set_ticks_visible(True)
        self.xax.set_ticklabel_visible(True)
        self.xax.set_axislabel(r"$\alpha_{2000}$",fontsize=14)
        self.xax.set_ticks_position('b')

    def show_yaxis(self,right=False):
        if right:
            self.yax.set_ticks_position('r')
            self.yax.set_ticklabel_position('r')
            self.yax.set_axislabel_position('r')

        self.yax.set_ticks_visible(True)
        self.yax.set_ticklabel_visible(True)
        self.yax.set_axislabel(r"$\delta_{2000}$",fontsize=14)

    def tight_layout(self):
        plt.tight_layout()

    def set_ax_lims(self,xlim,ylim,wcs=False):
        if not wcs:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
        else:
            bl = self.w.all_world2pix([[xlim[0],ylim[0]]],1)
            tr = self.w.all_world2pix([[xlim[1],ylim[1]]],1)
            #br = self.w.all_pix2world([[xlim[1],ylim[0]]],1)
            #tl = self.w.all_pix2world([[xlim[0],ylim[1]]],1)

            xlim = [bl[0],tr[0]]
            ylim = [bl[1],tr[1]]
            self.set_ax_lims(xlim,ylim)
            
        
    def get_ax_lims(self):
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        bl = self.w.all_pix2world([[xlim[0],ylim[0]]],1)
        tr = self.w.all_pix2world([[xlim[1],ylim[1]]],1)
        br = self.w.all_pix2world([[xlim[1],ylim[0]]],1)
        tl = self.w.all_pix2world([[xlim[0],ylim[1]]],1)
        
        return xlim,ylim,(bl,br,tl,tr)

    def transform_lims(self,wcs):
        _,_,coords = self.get_ax_lims()
        #print coords
        bl,br,tl,tr = coords
        bl1 = wcs.all_world2pix(bl,1)[0]
        tr1 = wcs.all_world2pix(tr,1)[0]
        #print bl1, tr1
        xlim = [bl1[0],tr1[0]]
        ylim = [bl1[1],tr1[1]]
        return np.array(xlim),np.array(ylim)

    def set_ticksize(self,size=8):
        try:
            self.ax.coords['ra'].set_ticks(exclude_overlapping=True,size=size)
            self.ax.coords['dec'].set_ticks(exclude_overlapping=True,size=size)

            #[tick.label.set_fontsize(size) for tick in self.ax.coords['ra'].get_major_ticks()]
            #[tick.label.set_fontsize(size) for tick in self.ax.coord['dec'].get_major_ticks()]
            self.ax.coords['ra'].set_ticklabel(size=size)
            self.ax.coords['dec'].set_ticklabel(size=size)

        except KeyError:
            self.ax.coords['glon'].set_ticks(exclude_overlapping=True,size=size)
            self.ax.coords['glat'].set_ticks(exclude_overlapping=True,size=size)
            self.ax.coords['glon'].set_ticklabel(size=size)
            self.ax.coords['glat'].set_ticklabel(size=size)


    def show(self):
        plt.show()


class North(object):
    def __init__(self,filename,xlim=None,ylim=None,fig=None,vmin=None,vmax=None,twin=None,cmap=None,right=False,aspect='auto',gs=None,angle=0.0,interp=None,shift=None,norm=None,plot=True):
        self.filename = filename
        self.fig = fig
        self.data,self.h0 = fits.getdata(self.filename,header=True,ext=0)
        if self.data.ndim == 3:
            self.data = self.data[0]

        # For dumb-ass PACS images
        if 'CUNIT1' in self.h0 and self.h0['CUNIT1'] in ['Degree','Degrees']:
            self.h0 = fits.getheader(self.filename,ext=1)
            self.h0['CUNIT1'] = ('deg','WCS unit1')
            self.h0['CUNIT2'] = ('deg','WCS unit2')

            self.w0 = WCS(self.h0)
            self.h0 = self.w0.to_header()
            self.h0['NAXIS1'] = self.data.shape[1]
            self.h0['NAXIS2'] = self.data.shape[0]
            self.h0['NAXIS'] = 2

        elif self.h0['NAXIS'] == 3:
            self.h0['NAXIS'] = 2
            del self.h0['NAXIS3']
            self.w0 = WCS(self.h0,naxis=2)

        else:
            self.w0 = WCS(self.h0)
            
        self.xlim = xlim
        self.ylim = ylim
        self.vmin = vmin
        self.vmax = vmax
        self.cmap = cmap
        self.twin = twin
        self.right = right
        self.aspect = aspect
        self.gs = gs
        self.angle = angle
        self.interp = interp
        self.shift = shift
        self.norm = norm

        self._rotate()

        if plot:
            self._make_axis()
        

    def _rotate(self):
        # Shift the WCS reference pixel to center of grid
        nx, ny = np.float(self.h0["naxis1"]), np.float(self.h0["naxis2"])
        i0, j0 = (nx + 1.0)/2, (ny + 1.0)/2  # 1 because indexing
        if self.h0['CTYPE1'] == 'GLON-CAR':
            pass
        elif 'PACS' in self.filename:
            self.h1 = self.h0.copy()
        else:
            [ra0, dec0], = self.w0.wcs_pix2world([[i0, j0]], 1)
            self.h0.update(crpix1=i0, crpix2=j0, crval1=ra0, crval2=dec0)
        self.nx = nx
        self.ny = ny

        # Construct header in equatorial frame
        self.h1 = self.h0.copy()
        
        try:
            self.h1.update(
                cd1_1=-np.hypot(self.h0["CD1_1"], self.h0["CD1_2"]), 
                cd1_2=0.0, 
                cd2_1=0.0, 
                cd2_2=np.hypot(self.h0["CD2_1"], self.h0["CD2_2"]), 
                orientat=0.0)
        except KeyError:
            self.h1.update(
                cd1_1=-np.hypot(self.h0["CDELT1"], self.h0["CDELT2"]), 
                cd1_2=0.0, 
                cd2_1=0.0, 
                cd2_2=np.hypot(self.h0["CDELT1"], self.h0["CDELT2"]), 
                orientat=0.0)            

        ### CURRENTLY BROKEN
        if self.angle == 90:
            print 'flipped'
            self.h1.update(
                cd1_1=np.hypot(self.h0["CD1_1"], self.h0["CD1_2"]),
                cd1_2=0.0,
                cd2_1=0.0,
                cd2_2=np.hypot(self.h0["CD2_1"], self.h0["CD2_2"]),
                orientat=0.0)

        if self.shift:
            self.data = ndimage.interpolation.shift(self.data,self.shift)

    def _make_axis(self):
        if self.fig is None:
            self.fig = plt.figure()

        if self.cmap is None:
            self.cmap = plt.cm.gray_r

        if self.twin and not self.right:
            if self.gs:
                gs = self.gs
            else:
                gs = 121
            self.ax = pywcsgrid2.subplot(gs, wcs=self.h1, aspect=self.aspect,adjustable='box',)
            
        elif self.twin and self.right:
            if self.gs:
                gs = self.gs
            else:
                gs = 122
            self.ax = pywcsgrid2.subplot(gs, wcs=self.h1, aspect=self.aspect,adjustable='box')
            
        elif self.gs:
            if 'PACS' in self.filename:
                self.ax = plt.subplot(self.gs,projection=self.w0)
            else:
                self.ax = pywcsgrid2.subplot(self.gs, wcs=self.h1, aspect=self.aspect,adjustable='box')

            #if self.h0['CTYPE1'] == 'GLON-CAR':
            #    self.ax.set_display_coord_system('fk5')
        else:
            self.ax = pywcsgrid2.axes(wcs=self.h1, aspect=self.aspect) #aspect=1
        self.ax.set_xlim(-self.nx/4, 5*self.nx/4)
        self.ax.set_ylim(-self.ny/4, 5*self.ny/4)

        if self.interp is None:
            if 'PACS' in self.filename:
                self.im = self.ax.imshow(self.data,origin='lower',vmin=self.vmin,vmax=self.vmax,cmap=self.cmap,norm=self.norm)
            else:
                self.im = self.ax[self.h0].imshow_affine(self.data,origin='lower',vmin=self.vmin,vmax=self.vmax,cmap=self.cmap,norm=self.norm)
        else:
            self.im = self.ax[self.h0].imshow(self.data,origin='lower',vmin=self.vmin,vmax=self.vmax,cmap=self.cmap,interpolation=self.interp,norm=self.norm)

        if self.xlim:
            self.ax.set_xlim(self.xlim)
        if self.ylim:
            self.ax.set_ylim(self.ylim)

        self.ax.grid(color='black', alpha=0.2, linestyle='dashed')

        if self.right:
            self.ax.axis['left'].major_ticklabels.set_visible(False)
            #self.ax.axis['right'].major_ticklabels.set_visible(True)
            self.ax.axis['right'].label.set_text(r"$\delta_{2000}$")
            plt.setp(self.ax.axis['right'].label,fontsize=14)
            self.ax.axis['left'].label.set_text('')

        else:
            try:
                plt.setp(self.ax.axis['left'].label,fontsize=14)
            except TypeError:
                plt.setp(self.ax.yaxis.label,fontsize=14)
        try:
            plt.setp(self.ax.axis['bottom'].label,fontsize=14)
        except TypeError:
            plt.setp(self.ax.xaxis.label,fontsize=14)
            
        return self.ax

    def show_yaxis(self,right=False):
        if right:
            self.ax.axis['left'].major_ticklabels.set_visible(False)
            self.ax.axis['right'].label.set_text(r"$\delta_{2000}$")
            plt.setp(self.ax.axis['right'].label,fontsize=14)
            self.ax.axis['left'].label.set_text('')
        else:
            self.ax.axis['left'].major_ticklabels.set_visible(True)
            self.ax.axis['left'].label.set_text(r"$\delta_{2000}$")
            plt.setp(self.ax.axis['left'].label,fontsize=14)



    def set_ticksize(self,size=8):
        [tick.label.set_fontsize(size) for tick in self.ax.xaxis.get_major_ticks()]
        [tick.label.set_fontsize(size) for tick in self.ax.yaxis.get_major_ticks()]
    
    def square(self,aspect=None):
        #self.ax.set_aspect(1./self.ax.get_data_ratio())
        if aspect:
            self.ax.set_aspect(aspect)
        else:
            self.ax.set_aspect('equal', adjustable='box')

            
    def get_aspect(self):
        return self.ax.get_data_ratio()

    def get_ax_size(self,pixels=False):
        bbox = self.ax.get_window_extent().transformed(self.fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        if pixels:
            return width*self.fig.dpi, height*self.fig.dpi
        else:
            return width, height

    def get_ax_lims(self):
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        wcs = WCS(self.h1)
        bl = wcs.all_pix2world([[xlim[0],ylim[0]]],1)
        tr = wcs.all_pix2world([[xlim[1],ylim[1]]],1)
        br = wcs.all_pix2world([[xlim[1],ylim[0]]],1)
        tl = wcs.all_pix2world([[xlim[0],ylim[1]]],1)
        
        return xlim,ylim,(bl,br,tl,tr)

    def transform_lims(self,wcs):
        _,_,coords = self.get_ax_lims()
        print coords
        bl,br,tl,tr = coords
        bl1 = wcs.wcs_world2pix(bl,1)[0]
        tr1 = wcs.wcs_world2pix(tr,1)[0]
        print bl1, tr1
        xlim = [bl1[0],tr1[0]]
        ylim = [bl1[1],tr1[1]]
        return np.array(xlim),np.array(ylim)


    def tight_layout(self):
        plt.tight_layout()

    def get_header(self):
        return self.h1
    
    def show(self):
        plt.show()
