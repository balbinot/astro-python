#! /usr/bin/env python
#05/21/2016
#utility functions for astronomy
import re
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from astropy.table import Table,Column,vstack
import matplotlib.pyplot as plt
from astropy.wcs.utils import proj_plane_pixel_scales
from pyraf import iraf
from astropy.io import fits
import tempfile
from collections import OrderedDict,namedtuple
import numpy as np
import os
import argparse
#from astropy.utils.console import ProgressBar

class RadProf(object):
    # Output columns from iraf.rimexam
    R_COLMAP = ['prad', 'frad', 'flux', 'bgd', 'mpeak', 'ellip', 'pa', 'mbeta', '2dfwhm', 'mfwhm', 'fwhm']

    #final table names/Profile class names
    T_COLMAP = ['xi','yi','prad', 'frad', 'flux', 'bgd', 'mpeak', 'ellip', 'pa', 'mbeta','mfwhm','mfwhm_s','fwhm','fwhm_s']

    #GKIDECODE expressions
    fitRe = re.compile('color=1\s*\n\s*(?:polyline.*\s*([\d\.\d*\s]*))*flush')
    pairRe = re.compile('(\d\.\d*\s\d\.\d*)')
    dataRe = re.compile('color=1\n(?:(?:polyline.*\n)\s(.*\n))+flush')
    extractRe = re.compile('xavg=(\d\.\d\d).*yavg=(\d\.\d\d)')
    ndcRe = re.compile('set_wcs.*\n(.*)')

    # tempfiles for pyraf
    tempdatafile = None
    tempcoofile = None
    tempradfile = None

    def __init__(self,data=None,header=None,wcs=None):
        # if data is filename, read from file
        if isinstance(data,str):
            self.filename = data
            self.data,self.header = fits.getdata(data,header=True)
            self.wcs = WCS(self.header)
            
        else:
            # data is array.  make tempfile to store data
            self.data = data
            self.header = header
            self.wcs = wcs
            hdul = fits.PrimaryHDU(self.data)
            hdul.header = self.header
            self.tempdatafile = self.__generate_temp_file(data=hdul,suffix='.fits',ref=self.tempdatafile)
            self.filename = self.tempdatafile.name

        if self.wcs:
            # read pxscale from wcs
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

        # Generate tempfiles for pyraf
        self.tempcoofile = self.__generate_temp_file(suffix='.coo',ref=self.tempcoofile)
        self.tempradfile = self.__generate_temp_file(suffix='.prf',ref=self.tempradfile)


    def __generate_temp_file(self,data=None,suffix='',ref=None):
        if ref is not None:
            try:
                os.remove(ref)
            except:
                pass
                
        with tempfile.NamedTemporaryFile(prefix='utils_tmp',suffix=suffix,delete=False) as f:
            if data:
                data.writeto(f)
        return f

    def __clear_temp_files(self):
        if self.tempdatafile is not None:
            try:
                os.remove(self.tempdatafile.name)
            except:
                pass

        if self.tempcoofile is not None:
            try:
                os.remove(self.tempcoofile.name)
            except:
                pass

        if self.tempradfile is not None:
            try:
                os.remove(self.tempradfile.name)
            except:
                pass


    def _write_coord_file(self,x,y):
        with open(self.tempcoofile.name,'w') as fcoo:
            fcoo.write(' '.join([x,y]))
        return self.tempcoofile.name

    def __del__(self):
        self.__clear_temp_files()

    @staticmethod
    def extract_gkidata(text):
        try:
            fitdata = RadProf.fitRe.findall(text)
            fitdata = RadProf.pairRe.findall(fitdata[1]) # get pair vals
            fitdata = [x.split() for x in fitdata]
            fitdata = [(np.float(x[0]),np.float(x[1])) for x in fitdata]
            xf,yf = zip(*fitdata)
        except IndexError:
            #no fit
            xf,yf = (None,None)

        try:
            data = RadProf.dataRe.search(text).group(0)
            data = RadProf.extractRe.findall(data)
            data = [(np.float(x[0]),np.float(x[1])) for x in data]
        except AttributeError:
            return None,None,None,None

        radius,pixelval = zip(*data)

        #get units
        unitline = RadProf.ndcRe.findall(text)[0]
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
        if xf:
            xf = xfit(xf)
        return radius,np.array(pixelval),np.array(xf),np.array(yf)

    def _parse_coords(self,coords):
        if isinstance(coords,SkyCoord):
            coordlist = [x.to_string('hmsdms',sep=':').split() for x in coords]
        elif isinstance(coords,np.ndarray) or len(coords) > 2:
            coordlist = []
            for x in coords:
                if isinstance(x,SkyCoord):
                    coordlist.append(x.to_string('hmsdms',sep=':').split())
                else:
                    coordlist.append((str(x[0]),str(x[1])))
        elif len(coords) == 1:
            #[(x,y)]
            for x in coords:
                if isinstance(x,SkyCoord):
                    coordlist.append(x.to_string('hmsdms',sep=':').split())
                else:
                    coordlist.append((str(x[0]),str(x[1])))
        elif len(coords) == 2:
            coords = np.array(coords)
            if coords.ndim == 1:
                #2-tuple of x,y
                coordlist = [(str(coords[0]),str(coords[1]))]
            else:
                coordlist = [(str(x[0]),str(x[1])) for x in coords]
        else:
            raise ValueError('Cannot recognize data format')

        for idx,coord in enumerate(coordlist):
            x,y = coord
            if ':' in x:
                coord = SkyCoord(ra=x,dec=y,unit=(u.hourangle,u.deg))
                coord = coord.to_pixel(self.wcs)
                coordlist[idx] = ('%.3f'%coord[0],'%.3f'%coord[1])
        return coordlist

    
    def profile(self,coords,table=True):
        coordlist = self._parse_coords(coords)
        profs = []
        #for coord in ProgressBar(coordlist):
        for coord in coordlist:
            coordfile = self._write_coord_file(*coord)
            datum = self._profiler(coordfile,coord)
            #datum = Profile(xr,yr,xf,yf,tbl)
            profs.append(datum)

        table = vstack([p.tab for p in profs])
        table.pprint()
        if table:
            return profs,table
        else:
            if len(profs) == 1:
                return profs[0]
            else:
                return profs
        

    def _profiler(self,coordfile,coord):
        # deploy task to pyraf
        results = iraf.imexamine(self.filename,image=self.filename,use_display='no',logfile='',defkey='r',imagecur=coordfile,frame=1,mode='h',StdoutG=self.tempradfile.name,Stdout=1)
        try:
            results = results[0].split()
            # prad frad flux bgd mpeak ellip pa mbeta 2dfwhm mfwhm fwhm
        except IndexError:
            # some error
            results = [np.nan]*len(self.R_COLMAP)

        # first, process stdout
        for idx,r in enumerate(results):
            try:
                results[idx] = np.float(r)
            except ValueError:
                results[idx] = np.nan
            
        results = OrderedDict(zip(self.R_COLMAP,results))
        tbl = Table([results])
        tbl = tbl[self.R_COLMAP]
        tbl['flux'].unit = u.ct
        tbl['bgd'].unit = u.ct
        tbl['mpeak'].unit = u.ct
        tbl['pa'].unit = u.deg
        #tbl['2dfwhm'].unit = u.pix
        tbl['mfwhm'].unit = u.pix
        tbl['fwhm'].unit = u.pix

        # add xi,yi
        xi,yi = (np.float(i) for i in coord)
        tbl.add_column(Column([xi],name='xi',format='%.2f'),0)
        tbl.add_column(Column([yi],name='yi',format='%.2f'),1)

        #remove 2dfwhm cuz wut?
        tbl.remove_column('2dfwhm')

        #add arcsec columns
        for fwhmcol in tbl.colnames[-2:]:
            try:
                val = tbl[fwhmcol]*self.pxscale
            except TypeError:
                val = [np.nan]
            tbl.add_column(Column(val,name='%s_s'%fwhmcol,unit=u.arcsec,format='%.3f'),tbl.index_column(fwhmcol)+1)

        # process profile graphics
        gki = iraf.gkidecode(self.tempradfile.name,Stdout=1,verbose=True)
        gki = '\n'.join(gki)
        xr,yr,xf,yf = self.extract_gkidata(gki)

        vals = [tbl[x][0] for x in tbl.colnames] + [xr,yr,xf,yf,tbl]
        p = Profile(*vals)

        return p


class Profile(namedtuple('Profile',RadProf.T_COLMAP+['rad','int','radf','intf','tab'])):
    __slots__ = ()
    def __repr__(self):
        #rep = ','.join(['%s=%.1f' % (x,self.__dict__[x]) for x in RadProf.T_COLMAP])
        return 'Profile(xi=%.2f,yi=%.2f,fwhm=%r)' % (self.xi,self.yi,self.fwhm)

    

# testing
def main():
    parser = argparse.ArgumentParser(description = 'Astro utilities')
    parser.add_argument('filename',type=str,help='Input file for analysis')
    parser.add_argument('-coords',nargs=2,help='X Y coordinates')
    parser.add_argument('-task',default='radprof',help='Specify task (default=radprof')

    args = parser.parse_args()

    if args.task == 'radprof':
        r = RadProf(args.filename)
        #coords = [(1390,523),(1361,670)]
        #r.profile(coords)
        #coords = [(1390,523)]
        #r.profile(coords)
        #coords = (1390,670)
        #r.profile(coords)
        coords = [(1390,523),(1361,670),(1361,670),('9:13:22.286','76:29:07.48'),(1390,523),(1361,670),(1361,670),('9:13:22.286','76:29:07.48'),(1000,1000),SkyCoord('9:13:22.286','76:29:07.48',unit=(u.hourangle,u.deg))]
        profs,tbl = r.profile(coords)

if __name__ == '__main__':
    main()
