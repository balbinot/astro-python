#! /usr/bin/env python
# version 05/03/2016
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from cardelli import ext_correction,ext_correction_indebetouw
from mag2flux import get_flux,get_wave,get_flux_single,photcols,flux_to_mag
from astropy.table import Table,Column
from scipy.integrate import quad,simps
from scipy.interpolate import interp1d
import astropy.units as u

h = 6.626e-27  # erg*s
c = 3.0e10     # cm/s
k = 1.381e-16  # erg/K


class SED(object):
    def __init__(self,row,cols=None,xlim=None,ylim=None,fit=False,name=None,Av=0.0,extlaw='cardelli',scale=None):
        self.row = Table(row)
        self.name = name
        self.Av = Av
        self.extlaw = extlaw
        self.scale = scale
        
        if cols is not None:
            self.photcols = cols
        else:
            self.photcols = photcols

        self._add_phot_cols()
        if Av:
            self._add_ext_cols()
            
        self._make_figure(xlim,ylim,fit)

    def _add_phot_cols(self):
        self.waves = []
        self.energy = []
        self.filters = []
        
        for col in self.photcols:
            if col not in self.row.colnames:
                continue

            if self.row[col].unit == u.Jy:
                Jy_col,F_lam,lam_F_lam = get_flux(self.row[col],Jy=True)
            else:
                Jy_col,F_lam,lam_F_lam = get_flux(self.row[col])
                
            self.row.add_columns([Jy_col,F_lam,lam_F_lam])

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
        

    def _make_figure(self,xlim=None,ylim=None,fit=False):
        self.fig = plt.figure(dpi=100)

        if self.Av:
            plt.loglog(self.extwaves,self.extenergy,'s',markersize=5,
                       color='#707070',markeredgewidth=2,mec='#707070',zorder=2)
        else:
            self.extwaves = self.extwaves
            self.extenergy = self.energy
        plt.loglog(self.waves,self.energy,'ko',
                   markersize=8,markeredgewidth=2,markeredgecolor='k',zorder=9)

        if fit:
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
                if xlim:
                    waved = np.linspace(xlim[0],xlim[1],num=1000)
                else:
                    waved = np.linspace(np.min(self.waves).value/10.,np.max(self.waves).value*10.,num=1000)

                bb = blackbody(waved,*popt)
                if self.scale:
                    bb /= bb.max()
                    bb *= self.scale

                plt.loglog(waved,bb,linestyle='dashed',color='#707070',lw=1,zorder=2)
            
        plt.ylabel(r'$\lambda\,F_{\lambda}\,\left(ergs/s/cm^2\right)$',fontsize=14)
        plt.xlabel(r'$\lambda\,(\mu m)$',fontsize=14)

        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        else:
            plt.ylim([np.min(self.extenergy).value/10.,np.max(self.extenergy).value*10])
        if self.name:
            plt.title(self.name)
        majorticks = plt.gca().xaxis.get_majorticklocs()
        majortickvals = [major_ticks(x) for x in majorticks]
        plt.gca().set_xticklabels(majortickvals)

        minorticks = plt.gca().xaxis.get_minorticklocs()
        minortickvals = [minor_ticks(x) for x in minorticks]
        plt.gca().set_xticklabels(minortickvals,minor=True)

    def show(self):
        plt.show()

    def save(self,pp):
        pp.savefig(self.fig)
        plt.close(self.fig)
            

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
        return '${0:d}$'.format(int(value))
    
def minor_ticks(value):
    if value not in [0.2,0.5,2,5,20]:
        return ''

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
