#! /usr/bin/env python
from astroplan.plots import plot_finder_image
from astroplan import FixedTarget
import matplotlib.pyplot as plt
import argparse
from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord
import os.path
from matplotlib.backends.backend_pdf import PdfPages
from astropy.utils.console import ProgressBar
from astroquery.skyview import SkyView
from astropy.wcs import WCS
from astropy.visualization import PercentileInterval,LinearStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from urllib2 import HTTPError

class FindingChart(object):

    def __init__(self,coord,name=None,survey=None,radius=600,grid=False,ax=None,log=False,reticle=True,style_kwargs=None):
        self.coord = coord
        self.name = name
        self.survey = survey
        self.radius = radius
        self.target = FixedTarget(coord,name=name)
        self.grid = grid
        self.log = log
        self.reticle = reticle
        self.style_kwargs = style_kwargs
        self.hdu = None
        self.setstretch = None

        # set ax to false to suppress
        if ax != False:
            self._make_figure(ax=ax)

    def _make_figure(self,ax=None,right=False,style_kwargs=None):
        if ax:
            self.ax = ax
        else:
            self.ax = None
            plt.figure(dpi=100)

        if self.style_kwargs is None:
            self.style_kwargs = style_kwargs

        try:
            self.ax, self.hdu = plot_finder_image(self.target,
                                                  survey=self.survey,
                                                  fov_radius=self.radius*u.arcsec,
                                                  reticle=self.reticle,ax=self.ax,
                                                  grid=self.grid,log=self.log,
                                                  style_kwargs=self.style_kwargs)
        
            lon = self.ax.coords[0]
            lat = self.ax.coords[1]

            if right:
                #lat.set_axislabel('DEC',fontsize=10)
                lat.set_axislabel_position('r')
                lat.set_axislabel(r'$\delta_{2000}$',fontsize=14)
                lat.set_ticks_position('r')
                lat.set_ticklabel_position('r')
                #lon.set_axislabel('RA',fontsize=10)
                lon.set_axislabel(r'$\alpha_{2000}$',fontsize=14)

            lon.set_major_formatter('hh:mm:ss.s')
            lat.set_major_formatter('dd:mm:ss')

            if right:
                lat.ticklabels.set_fontsize(8)
                lon.ticklabels.set_fontsize(8)
                #self.ax.yaxis.set_label_position("right")
                #lat.set_ticklabel_position('r')
                #self.ax.set_ylabel(r'$\delta_{2000}$',fontsize=14)

        #except IndexError as e:
        except:
            self.ax = None
            return None


    def _get_projection(self):
        if self.hdu is None:
            try:
                hdu = SkyView.get_images(position=self.coord,
                                         coordinates='icrs',
                                         survey=self.survey,
                                         radius=self.radius*u.arcsec,
                                         grid=self.grid)[0][0]
                wcs = WCS(hdu.header)
                self.vlim = PercentileInterval(99.).get_limits(hdu.data)

            except (IndexError, HTTPError):
                hdu = SkyView.get_images(position=self.coord,
                                         coordinates='icrs',
                                         survey=self.survey,
                                         radius=self.radius*u.arcsec-1*u.arcsec,
                                         grid=self.grid)[0][0]
                wcs = WCS(hdu.header)
                self.vlim = PercentileInterval(99.).get_limits(hdu.data)
                
        else:
            wcs = WCS(self.hdu.header)
            self.vlim = PercentileInterval(99.).get_limits(self.hdu.data)
        return wcs

    def get_pix_coord(self,coord):
        return coord.to_pixel(wcs=self._get_projection())


    def show(self):
        plt.show()

    def save(self,pp):
        pp.savefig(plt.gcf())
        plt.close(plt.gcf())



def main():
    parser = argparse.ArgumentParser(description='Make finding charts')
    parser.add_argument('table',type=str,help='Input table')
    parser.add_argument('-pdf',type=str,help='If specified, output to directory')
    parser.add_argument('-fmt',type=str,default='ascii.tab',help='Specify table format for astropy (default = ascii.tab')
    parser.add_argument('-survey',type=str,default='2MASS-J',help='Pull images from this survey (default = 2MASS-J')
    parser.add_argument('-radius',type=float,default=120,help='FOV radius in arcsec (default = 120)')
    

    args = parser.parse_args()

    table = Table.read(args.table,format=args.fmt)

    if args.pdf:
        pp = PdfPages(args.pdf,keep_empty=False)
    else:
        pp = None
        
    for row in ProgressBar(table):
        coord = SkyCoord(ra=row['ra'],dec=row['dec'],unit=(u.hourangle,u.deg))
        if row['ID'] in ['S72','S52','S49','S09']:
            log = True
        else:
            log = False
        fc = FindingChart(coord,name=row['ID'],survey=args.survey,radius=args.radius,log=log)
        if pp:
            fc.save(pp)
        else:
            fc.show()

    if pp:
        print 'Writing finding charts to %s' % args.pdf
        pp.close()








        
        
if __name__ == '__main__':
    main()
