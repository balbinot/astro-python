#! /usr/bin/env python
from astropy.coordinates import SkyCoord,EarthLocation
from astropy.coordinates.errors import UnknownSiteException as USE
import astropy.units as u
from astroplan import Observer
from astropy.time import Time,TimezoneInfo
from astroplan.plots import plot_airmass
import matplotlib.pyplot as plt
import numpy as np
import pdb

class StarAlt(object):
    def __init__(self,coordinates,observatory='Mt Graham',night=None,timezone='UTC',ax=None,timedelta=None):
        if not isinstance(coordinates,SkyCoord):
            self.coord = SkyCoord(coordinates,unit=(u.hourangle,u.deg))
        else:
            self.coord = coordinates

        try:
            self.observer = Observer.at_site(observatory,timezone=timezone)
        except USE:
            print EarthLocation.get_site_names()
            raise

        if not night:
            self.time = Time.now()
        else:
            self.time = Time(night)

        if timedelta:
            self.time = timedelta

        self.ax = ax

        self._make_figure()


    def _make_figure(self,xlim=None,ylim=None,ax=None):
        if ax:
            self.ax = ax

        plot_airmass(self.coord,self.observer,self.time,ax=self.ax,
                     altitude_yaxis=True,brightness_shading=True)

        self.ax.tick_params(labelsize=8)
        #self._make_local_time()

    def _make_local_time(self,utcoffset=-7):
        xaxis = self.ax.xaxis
        new = self.ax.twiny()
        
        
        pdb.set_trace()

    def _set_ylim(self,limits):
        #get other axis
        for ax2 in self.ax.figure.axes:
            if ax2 == self.ax:
                continue
        self.ax.set_ylim(limits)

        altitude_ticks = np.array([90, 60, 50, 40, 35, 30, 20])
        airmass_ticks = 1./np.cos(np.radians(90. - altitude_ticks))

        ax2.set_yticks(airmass_ticks)
        ax2.set_yticklabels(altitude_ticks,fontsize=9)
        ax2.set_ylim(self.ax.get_ylim())

    def show(self):
        plt.tight_layout()
        plt.show()


        

def main():
    obs = StarAlt('17:51:27.7 -23:51:39.1','Subaru',night='2017-06-01',timezone='US/Hawaii')
    obs.show()
    
if __name__ == '__main__':
    main()
