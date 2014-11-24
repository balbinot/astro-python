#! /usr/bin/env python
from astropy.table import Table,Column
import numpy as np

## http://ssc.spitzer.caltech.edu/warmmission/propkit/pet/magtojy/ref.html
zp = {'U':1823,'B':4130,'V':3781,'R':2941,'I':2635,'3.6':277.5,'4.5':179.5,'5.8':116.6,'8.0':63.1,'J':1594,'H':1024,'K':666.7,'W1':309.540,'W2':171.787,'W3':31.674,'W4':8.363,'8':63.1}  #Jy

wave = {'U':0.36,'B':0.44,'V':0.55,'R':0.71,'I':0.97,'3.6':3.55,'4.5':4.439,'5.8':5.731,'8.0':7.872,'J':1.235,'H':1.662,'K':2.159,'W1':3.3526,'W2':4.6028,'W3':11.5608,'W4':22.0883,'8':7.872}  #microns

def get_wave(filter):
    return wave.get(filter)

def get_flux_single(name,val):
    try: 
        val = np.float(val)
    except:
        return (np.nan,np.nan,np.nan)
        
    if val and val < 90:
        Jy = np.power(10.0,-val/2.5)*zp[name]
        F_lam = 3.0e-9 * Jy / (wave[name]**2)
        lam_F_lam = wave[name] * F_lam
        return  (Jy,F_lam,lam_F_lam)

    return (np.nan,np.nan,np.nan)


def get_flux(tablecol,name_ext = ''):
    '''Returns [F_Jy,F_lambda,lam_F_lam]'''
    name = tablecol.name

    Jy_col = [np.power(10.0,-val/2.5)*zp[name] if val else np.nan for val in tablecol.data]
    #Jy_col = [np.float(x) if x else np.nan for x in Jy_col]
   
    # convert to F_lambda, erg/s/cm^2/micron
    F_lam = [3.0e-9 * val / (wave[name]**2) if val else np.nan for val in Jy_col]
    #F_lam = [np.float(x) if x else None for x in F_lam]

    
    lam_F_lam = [val * wave[name] if val else np.nan for val in F_lam]
    #lam_F_lam = [np.float(x) if x else 99.99 for x in lam_F_lam]

    Jy_col = [x if not np.isnan(x) else 99.99 for x in Jy_col]
    F_lam = [x if not np.isnan(x) else 99.99 for x in F_lam]
    lam_F_lam = [x if not np.isnan(x) else 99.99 for x in lam_F_lam]
    
    Jy_col = Column(Jy_col,name='F_%s_Jy%s' % (name,name_ext),description='Zeropoint: %i Jy' % zp[name],unit='Jy')
    F_lam = Column(F_lam,name='F_%.2f_um%s' % (wave[name],name_ext),unit='erg*s^-1*cm^-2*micron^-1',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[name],wave[name]))
    lam_F_lam = Column(lam_F_lam,name='lam_F_%.2f_um%s' % (wave[name],name_ext),unit='erg*s^-1*cm^-2',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[name],wave[name]))

        
    return [Jy_col,F_lam,lam_F_lam]



