#! /usr/bin/env python
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import pywcsgrid2
import matplotlib.pyplot as plt

class North(object):
    def __init__(self,filename,xlim=None,ylim=None,fig=None,vmin=None,vmax=None,twin=None,cmap=None,right=False,aspect='auto',gs=None):
        self.filename = filename
        self.fig = fig
        self.data,self.h0 = fits.getdata(self.filename,header=True)
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

        self._rotate()
        
        self._make_axis()
        

    def _rotate(self):
        # Shift the WCS reference pixel to center of grid
        nx, ny = np.float(self.h0["naxis1"]), np.float(self.h0["naxis2"])
        i0, j0 = (nx + 1.0)/2, (ny + 1.0)/2  # 1 because indexing
        [ra0, dec0], = self.w0.wcs_pix2world([[i0, j0]], 1)
        self.h0.update(crpix1=i0, crpix2=j0, crval1=ra0, crval2=dec0)
        self.nx = nx
        self.ny = ny

        # Construct header in equatorial frame
        self.h1 = self.h0.copy()
        self.h1.update(
            cd1_1=-np.hypot(self.h0["CD1_1"], self.h0["CD1_2"]), 
            cd1_2=0.0, 
            cd2_1=0.0, 
            cd2_2=np.hypot(self.h0["CD2_1"], self.h0["CD2_2"]), 
            orientat=0.0)

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
            
        else:
            self.ax = pywcsgrid2.axes(wcs=self.h1, aspect=self.aspect) #aspect=1
        self.ax.set_xlim(-self.nx/4, 5*self.nx/4)
        self.ax.set_ylim(-self.ny/4, 5*self.ny/4)
        self.ax[self.h0].imshow_affine(self.data,origin='lower',vmin=self.vmin,vmax=self.vmax,cmap=self.cmap)

        if self.xlim:
            self.ax.set_xlim(self.xlim)
        if self.ylim:
            self.ax.set_ylim(self.ylim)

        self.ax.grid()

        if self.right:
            self.ax.axis['left'].major_ticklabels.set_visible(False)
            #self.ax.axis['right'].major_ticklabels.set_visible(True)
            self.ax.axis['right'].label.set_text(r"$\delta_{2000}$")
            self.ax.axis['left'].label.set_text('')
        return self.ax

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

    def show(self):
        plt.show()
