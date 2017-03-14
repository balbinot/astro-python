#! /usr/bin/env python
# version 05/03/2016
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from cardelli import ext_correction,ext_correction_indebetouw
from mag2flux import get_flux,get_wave,get_flux_single,photcols,flux_to_mag
from astropy.table import Table,Column,MaskedColumn
from scipy.integrate import quad,simps
from scipy.interpolate import interp1d
import astropy.units as u
import argparse
#from scipy.integrate import simps

h = 6.626e-27  # erg*s
c = 3.0e10     # cm/s
k = 1.381e-16  # erg/K


class SED(object):
    def __init__(self,row,cols=None,xlim=None,ylim=None,fit=None,name=None,Av=0.0,extlaw='cardelli',scale=None,ax=None,title=True,suppress_plot=False,ms=5):
        self.row = Table(row)
        self.origRow = self.row.copy()
        self.name = name
        self.Av = Av
        self.extlaw = extlaw
        self.scale = scale
        self.xlim = xlim
        self.ylim = ylim
        self.fit = fit
        self.ms = ms
        
        if title == True:
            self.title = name
        elif title:
            self.title = title
        else:
            self.title = False
            
        if cols is not None:
            self.photcols = cols
        else:
            self.photcols = photcols

        self._add_phot_cols()
        if Av and not np.isnan(Av):
            self._add_ext_cols()
        else:
            self.Av = 0

        self.Teff = None
            
        if not suppress_plot:
            self._make_figure(xlim,ylim,fit,ax)

    @classmethod
    def from_object_table(cls,table,ax=None,**kwargs):
        name = table['SED_name']
        Av = table['SED_Av']
        extlaw = table['SED_extlaw']
        scale = table['SED_scale']
        xlim = [table['SED_xmin'],table['SED_xmax']]
        ylim = [table['SED_ymin'],table['SED_ymax']]
        fit = table['SED_fit']
        if 'fit' in kwargs:
            fit = kwargs.pop('fit')
        if 'xlim' in kwargs:
            xlim = kwargs.pop('xlim')
        if 'ylim' in kwargs:
            ylim = kwargs.pop('ylim')
        if 'title' in kwargs:
            title = kwargs.pop('title')
        return cls(table,name=name,Av=Av,extlaw=extlaw,scale=scale,
                   xlim=xlim,ylim=ylim,fit=fit,ax=ax,title=title,**kwargs)

    @classmethod
    def from_astropy_table(cls,table,row=0):
        return cls(table[row])

    def make_object_table(self):
        #make new table so object can be restored
        ## caution:  cols and lims not added
        namecol = Column(data=[self.name],name='SED_name')
        avcol =   Column(data=[self.Av],name='SED_Av',unit=u.mag)
        extcol =  Column(data=[self.extlaw],name='SED_extlaw')
        if self.scale:
            scalecol =Column(data=[self.scale.value],name='SED_scale',unit=self.scale.unit,dtype=float)
        else:
            scalecol =MaskedColumn(data=[np.nan],name='SED_scale',dtype=float,mask=[True],fill_value=np.nan)
        xmincol = Column(data=[self.xlim[0]],name='SED_xmin')
        xmaxcol = Column(data=[self.xlim[1]],name='SED_xmax')
        ymincol = Column(data=[self.ylim[0]],name='SED_ymin')
        ymaxcol = Column(data=[self.ylim[1]],name='SED_ymax')
        if not self.fit:
            self.fit = -1
        fitcol = Column(data=[self.fit],name='SED_fit',dtype=float)
                        
        newtable = self.origRow.copy()
        newtable.add_columns([namecol,avcol,extcol,scalecol,
                              xmincol,xmaxcol,ymincol,ymaxcol,fitcol])
        if not 'Teff' in newtable.colnames:
            if self.Teff:
                newtable.add_column(Column([self.Teff],name='Teff',unit=u.K))

        if not 'Lum' in newtable.colnames:
            if hasattr(self,'Lum'):
                newtable.add_column(Column([self.Lum],name='Lum'))
        return newtable

        
        

    def _add_phot_cols(self):
        self.waves = []
        self.energy = []
        self.filters = []
        
        for col in self.photcols:
            if col not in self.row.colnames:
                continue

            if self.row[col].unit == u.Jy:
                Jy_col,F_lam,lam_F_lam = get_flux(self.row[col],Jy=True)
            elif self.row[col].unit == u.mJy:
                self.row[col] /= 1000
                self.row[col].unit = u.Jy
                Jy_col,F_lam,lam_F_lam = get_flux(self.row[col],Jy=True)
            elif self.row[col].unit == u.W/u.cm**2:
                ####NOT IMPLEMENTED YET
                Jy_col = Column([None]*len(self.row[col]),name='F_Jy_bad_%s'%col)
                F_lam = Column([None]*len(self.row[col]),name='F_lam_bad_%s'%col)
                lam_F_lam = Column(self.row[col].to(u.erg/u.cm**2/u.s).value,unit=self.row[col].to(u.erg/u.cm**2/u.s),name='lam_F_lam_good_%s'%col)
            
            else:
                Jy_col,F_lam,lam_F_lam = get_flux(self.row[col])

            try:
                self.row.add_columns([Jy_col,F_lam,lam_F_lam])
            except ValueError as e:
                #print 'ValueError: _add_phot_cols'
                #print e
                pass

            val = lam_F_lam
            energy = float(val) if (val < 70) and (val > 1e-40) else np.nan
            if not np.isnan(energy):
                self.waves.append(get_wave(col))
                self.energy.append(energy)
                self.filters.append(col)
        self.waves *= u.micron
        self.energy *= u.erg/u.s/u.cm**2

    def _add_ext_cols(self):
        self.extwaves = []
        self.extenergy = []
        self.extfilters = []
        if self.extlaw == 'indebetouw':
            extdict = ext_correction_indebetouw(self.Av)
        else:
            extdict = ext_correction(self.Av)
        for k,v in extdict.iteritems():
            if k in self.filters:
                if self.row[k].unit == u.Jy:
                    # have to quickly convert to mag and back to extinct
                    mag = flux_to_mag(self.row[k],k,Jy=True)
                    mag -= v
                    self.row[k],_,_ = get_flux_single(k,mag)
                else:
                    self.row[k] -= v

        for col in self.filters:
            if self.row[col].unit == u.Jy:
                _,_,lam_F_lam = get_flux_single(col,self.row[col],Jy=True)
            else:
                _,_,lam_F_lam = get_flux_single(col,self.row[col])

            val = lam_F_lam
            energy = float(val) if (val < 70) and (val > 1e-40) else np.nan
            if not np.isnan(energy):
                self.extwaves.append(get_wave(col))
                self.extenergy.append(energy)
                self.extfilters.append(col)
        self.extwaves *= u.micron
        self.extenergy *= u.erg/u.s/u.cm**2


    def observed_lum(self):
        ###make sure sorted
        #### try np.trapz() instead!!!
        flux = self.energy / self.waves
        tflux = simps(flux,self.waves)
        tflux *= flux.unit*self.waves.unit
        return np.abs(tflux)


    def get_mag_at(self,wave):
        idx = np.argmin(np.abs(wave*u.micron-self.waves))
        energy = self.energy[idx]
        flux = energy/self.waves[idx]
        filt = self.filters[idx]
        mag = flux_to_mag(flux.value,filt)
        
        return mag

    def get_projected_mag(self,band,scale=0):
        if self.popt is None:
            return None
        
        wave = get_wave(band)
        #fit = interp1d(self.extwaves,self.extenergy)
        #print blackbody(wave,*self.popt)
        #energy = fit(wave)
        energy = blackbody(wave,*self.popt)
        energy = float(energy)*u.erg/u.cm**2/u.s
        '''
        print energy
        print self.scale,'scale'
        if self.scale and not np.isnan(self.scale):
            waved = np.linspace(0.01,100,num=5000)
            bb = blackbody(waved,*self.popt)
            factor = bb.max()
            energy = (energy/factor) * self.scale
            if not hasattr(energy,'unit'):
                energy *= u.erg/u.cm**2/u.s
        '''
        #print energy,'final'
        wave *= u.micron
        flux = (energy/wave)
        #print flux
        flux = flux.to(u.Jy,equivalencies=u.spectral_density(wave))
        #print flux
        mag = flux_to_mag(flux.value,band,Jy=True)

        if self.extlaw == 'indebetouw':
            Ak = self.Av
            if band == 'I':
                Av = Ak/0.479
            elif band == 'R':
                Av = Ak/0.751
            else:
                Av = Ak
            ext = ext_correction(Av,filter=band)
            '''
            filter = band if wave.value >= get_wave('J') else 'J'
            ext = ext_correction_indebetouw(self.Av,filter=filter)
            if filter == 'R':
                ext *= 2.5
            elif filter == 'I':
                ext *= 1.7
            '''
        else:
            ext = ext_correction(self.Av,filter=band)
        print mag,ext
        _,_,newenergy = get_flux_single(band,mag+ext+scale)
        #print wave,newenergy
        return (mag+ext+scale, newenergy)
            

            

    def _make_figure(self,xlim=None,ylim=None,fit=None,ax=None):
        if ax:
            self.ax = ax
        else:
            self.ax = plt.figure(dpi=100).gca()

        if self.Av:
            self.exthandle = self.ax.loglog(self.extwaves,self.extenergy,'s',markersize=self.ms,color='#707070',markeredgewidth=2,mec='#707070',zorder=2)
        else:
            self.extwaves = self.waves
            self.extenergy = self.energy
        self.energyhandle = self.ax.loglog(self.waves,self.energy,'ko',
                       markersize=self.ms,markeredgewidth=2,markeredgecolor='k',zorder=9)
        
        if fit is None:
            fit = self.fit
        if fit:
            if fit < 0:
                fit = True
            if fit == True:
                try:
                    popt,pcov = curve_fit(blackbody,self.extwaves,self.extenergy,p0=[4000,np.max(self.extenergy).value])
                except RuntimeError:
                    print 'No fit'
                    popt = None
            else:
                # fit is value, so assign temp
                try:
                    bounds = [(fit-10,-np.inf),(fit+10,np.inf)]
                    popt,pcov = curve_fit(blackbody,self.extwaves,self.extenergy,p0=[fit,np.max(self.extenergy).value],bounds=bounds)
                except RuntimeError:
                    print 'No fit'
                    popt = None

            if popt is not None:
                self.popt = popt
                self.Teff = popt[0]
                if xlim:
                    waved = np.linspace(xlim[0],xlim[1],num=1000)
                else:
                    waved = np.linspace(np.min(self.waves).value/10.,np.max(self.waves).value*10.,num=1000)

                bb = blackbody(waved,*popt)
                if self.scale and not np.isnan(self.scale):
                    bb /= bb.max()
                    bb *= self.scale
                self.Flux = 1.3586*bb.max()
                self.Flux *= u.erg/u.s/(u.cm**2)

                self.fithandle = self.ax.loglog(waved,bb,linestyle='dashed',color='#707070',lw=1,zorder=2)
            else:
                self.popt = None
                self.Teff = None
                self.Flux = None
            
        self.ax.set_ylabel(r'$\lambda\,F_{\lambda}\,\left(\mathrm{ergs/s/cm}^2\right)$',fontsize=14)
        self.ax.set_xlabel(r'$\lambda\,\left(\mu\mathrm{m}\right)$',fontsize=14)

        if xlim:
            self.ax.set_xlim(xlim)
        else:
            if self.xlim:
                self.ax.set_xlim(self.xlim)
                   
        if ylim:
            self.ax.set_ylim(ylim)
        elif self.ylim:
            self.ax.set_ylim(self.ylim)
        else:
            self.ax.set_ylim([np.min(self.extenergy).value/10.,np.max(self.extenergy).value*10])

        self.xlim = self.ax.get_xlim()
        self.ylim = self.ax.get_ylim()
            
        if self.title:
            self.ax.set_title(self.title)

        self.fix_axes()

    def fix_axes(self):
        #self.ax.tick_params(axis='x',which='both',pad=5)
        majorticks = self.ax.xaxis.get_majorticklocs()
        majortickvals = [major_ticks(x) for x in majorticks]
        self.ax.set_xticklabels(majortickvals)

        minorticks = self.ax.xaxis.get_minorticklocs()
        minortickvals = [minor_ticks(x) for x in minorticks]
        self.ax.set_xticklabels(minortickvals,minor=True)

    def show(self):
        plt.show()

    def save(self,pp):
        pp.savefig(plt.gcf())
        plt.close(plt.gcf())
    def save_png(self,filename):
        plt.savefig(filename)


    def as_table(self):
        outtable = Table([self.waves,self.energy],names=['lambda','lambda_F_lam'])
        outtable.sort(keys='lambda')
        return outtable

    def as_dict(self):
        return {x:z for x,z in zip(self.filters,self.energy)}

    def __str__(self):
        return str([(z,'%0.3f %s'%(x.value,x.unit),'%0.3g %s'%(y.value,y.unit)) for z,x,y in zip(self.filters,self.waves,self.energy)])
    

def blackbody(waves, T, scale=np.pi):
    # waves in microns

    # convert to cm
    wavesCM = waves*1e-4

    # flux in ergs/s/cm^2/cm/sr
    flux = 2*h*c**2/wavesCM**5 * 1.0/(np.exp(h*c/(k*wavesCM*T)) - 1.0)

    # rescale for fitting
    flux = scale*flux

    # return as intensity, ergs/s/cm^2
    return flux * wavesCM



def correct_phot(phot,Av,Rv=3.2):
    ext_dicts = [ext_correction(Av,Rv=Rv)]
    for i,z in enumerate(zip(phot,ext_dicts)):
        row,ext_dict = z
        for k,v in ext_dict.iteritems():
            if v is None:
                continue
            if (k in row.colnames) and (row[k] or row[k] != 0) and (row[k] != '--') and (float(row[k]) < 70) and (not np.isnan(float(row[k]))):
                phot[k][i] = float(phot[k][i]) - v

    return phot

def major_ticks(value):
    if value < 1:
        return ''
    else:
        return '{0:d}'.format(int(value))
    
def minor_ticks(value):
    if value not in [0.2,0.5,2,5,20]:
        return ''
    if value == 0.5:
        return '{0:.1f}'.format(0.5)
    if value in [2,5,20]:
        return '{0:d}'.format(int(value))

    exp = np.floor(np.log10(value))
    base = value/10**exp
    if exp == 0 or exp == 1:   
        return '${0:d}$'.format(int(value))
    if exp == -1:
        return '${0:.1f}$'.format(value)
    else:
        return '${0:d}\\times10^{{{1:d}}}$'.format(int(base), int(exp))

def add_flux_cols(phot,cols,ext=''):
    for col in cols:
        Jy_col,F_lam,lam_F_lam = get_flux(phot[col],name_ext=ext)
        phot.add_columns([Jy_col,F_lam,lam_F_lam])
    return phot


def main():
    parser = argparse.ArgumentParser(description='Test SED object')
    parser.add_argument('store',type=str,help='SED object FITS data store')
    args = parser.parse_args()

    table = Table.read(args.store)
    print table
    print table.colnames
    for row in table:
        sed = SED.from_object_table(row)
        sed.show()


if __name__ == '__main__':
    main()
