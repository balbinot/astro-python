#! /usr/bin/env python
from astropy.table import Table,Column
import numpy as np
import astropy.units as u

photcols = ['U','B','V','R','Rc','I','Ic','3.6','4.5','5.8','8.0','J','H','K','L','W1','W2','W3','W4','u','g','r','i','z','450w','814w','IRAS12','IRAS25','IRAS60','IRAS100','MIPS24','MIPS70','MIPS100','AKARIS9W','AKARIL18','AKARIN60','AKARIWIDES','AKARIWIDEL','MSXA','MSXB1','MSXB2','MSXC','MSXD','MSXE','FOR_F197','FOR_F253','FOR_F315','FOR_F348','FOR_F371','19.7', '25.358', '31.5', '31.5', '34.8', '34.8', '37.1','PACS70','PACS100','PACS160','SPIRE250','SPIRE350','SPIRE500','LABOCA870']

## http://ssc.spitzer.caltech.edu/warmmission/propkit/pet/magtojy/ref.html
##http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=CFHT&gname2=MegaCam
zp = {'U':1823,'B':4130,'V':3781,'R':2941,'Rc':3028,'I':2635,'Ic':2458.3,'L':278,'3.6':277.5,'4.5':179.5,'5.8':116.6,'8.0':63.1,'J':1594,'H':1024,'K':666.7,'W1':309.540,'W2':171.787,'W3':31.674,'W4':8.363,'8':63.1,'u':2669.7,'g':3920.0,'r':3126.8,'i':2581.8,'z':2255.8,'450w':3930.0,'814w':2372.1,'IRAS12':30.9,'IRAS25':7.3,'IRAS60':1.1,'IRAS100':0.4,'MIPS24':7.140,'MIPS70':0.775,'MIPS100':0.159,'AKARIS9W':47.5,'AKARIL18':10.7,'AKARIN60':0.9,'AKARIWIDES':0.5,'AKARIWIDEL':0.2,'MSXA':55.810,'MSXB1':193.613,'MSXB2':189.003,'MSXC':26.423,'MSXD':18.267,'MSXE':8.666,'FOR_F197':0,'FOR_F253':0,'FOR_F315':0,'FOR_F348':0,'FOR_F371':0,'PACS70':0.8,'PACS100':0.4,'PACS160':0.1,'SPIRE250':0.1,'SPIRE350':0.0,'SPIRE500':0,'LABOCA870':0}  #Jy

wave = {'U':0.36,'B':0.44,'V':0.55,'R':0.71,'Rc':0.6358,'I':0.97,'Ic':0.78292,'L':3.4500,'3.6':3.55,'4.5':4.439,'5.8':5.731,'8.0':7.872,'J':1.235,'H':1.662,'K':2.159,'W1':3.3526,'W2':4.6028,'W3':11.5608,'W4':22.0883,'8':7.872,'u':0.38816,'g':0.47670,'r':0.61917,'i':0.74674,'z':0.88240,'450w':0.44441,'814w':0.81864, '24':24.,'19.7':19.7,'25.358':25.358,'31.5':31.5,'37.1':37.1,'70':70,'100':100,'IRAS12':10.14646,'IRAS25':21.72655,'IRAS60':51.98873,'IRAS100':95.29710,'MIPS24':23.680,'MIPS70':71.420,'MIPS100':155.900,'AKARIS9W':8.22836,'AKARIL18':17.60949,'AKARIN60':62.95067,'AKARIWIDES':76.90352,'AKARIWIDEL':140.85611,'MSXA':8.280,'MSXB1':4.290,'MSXB2':4.350,'MSXC':12.130,'MSXD':14.650,'MSXE':21.340,'FOR_F197':19.7,'FOR_F253':25.358,'FOR_F315':31.5,'FOR_F348':34.8,'FOR_F371':37.1,'PACS70':68.92474,'PACS100':97.90361,'PACS160':153.94392,'SPIRE250':247.12451,'SPIRE350':346.71804,'SPIRE500':496.10677,'LABOCA870':870.0}  #microns

def get_waves(filters):
    return [wave.get(f) for f in list(filters)]

def get_wave(filter):
    return wave.get(filter)

def get_flux_single(name,val,Jy=False):
    try: 
        val = np.float(val)
    except:
        return (np.nan,np.nan,np.nan)
        
    if val and val < 90:
        if Jy:
            Jyval = val
        else:
            Jyval = np.power(10.0,-val/2.5)*zp[name]
        F_lam = 3.0e-9 * Jyval / (wave[name]**2)
        lam_F_lam = wave[name] * F_lam
        return  (Jyval,F_lam,lam_F_lam)

    return (np.nan,np.nan,np.nan)


def flux_to_mag(F_lam,name,Jy=False):
    #convert to Jy
    if not Jy:
        Jy = wave[name]**2 * F_lam / 3.0e-9
    else:
        Jy = F_lam
        
    #convert to mag
    if hasattr(Jy,'unit'):
        try:
            val = -2.5 * np.log10(Jy/(zp[name]*u.Jy))
        except TypeError:
            dat = np.array([float(x) if x and x != '--' else np.nan for x in Jy])
            val = -2.5 * np.log10(dat/zp[name])
    else:
        val = -2.5 * np.log10(Jy/zp[name])
    return val

def get_flux(tablecol,name_ext = '',name=None,Jy=False):
    '''Returns [F_Jy,F_lambda,lam_F_lam]'''

    if name is None:
        name = tablecol.name

    data = [x for x in tablecol.data]
    for idx,val in enumerate(data):
        try:
            data[idx] = np.float(val)
        except:
            data[idx] = None

    # if data already in Jy, just clean it up
    if Jy:
        Jy_col = [np.float(val) if val else np.nan for val in data]
    else:
        Jy_col = [np.power(10.0,-val/2.5)*zp[name] if val else np.nan for val in data]
   
    # convert to F_lambda, erg/s/cm^2/micron
    F_lam = [3.0e-9 * val / (wave[name]**2) if val else np.nan for val in Jy_col]
    #F_lam = [np.float(x) if x else None for x in F_lam]

    
    lam_F_lam = [val * wave[name] if val else np.nan for val in F_lam]
    #lam_F_lam = [np.float(x) if x else 99.99 for x in lam_F_lam]

    Jy_col = [x if not np.isnan(x) else 99.99 for x in Jy_col]
    F_lam = [x if not np.isnan(x) else 99.99 for x in F_lam]
    lam_F_lam = [x if (not np.isnan(x)) and (x > 1e-30) else 99.99 for x in lam_F_lam]
    
    Jy_col = Column(Jy_col,name='F_%s_Jy%s' % (wave[name],name_ext),description='Zeropoint: %i Jy' % zp[name],unit='Jy')
    #    Jy_col = Column(Jy_col,name='F_%s_Jy%s' % (name,name_ext),description='Zeropoint: %i Jy' % zp[name],unit='Jy')
    F_lam = Column(F_lam,name='F_%.2f_um%s' % (wave[name],name_ext),unit='erg*s^-1*cm^-2*micron^-1',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[name],wave[name]))
    lam_F_lam = Column(lam_F_lam,name='lam_F_%.2f_um%s' % (wave[name],name_ext),unit='erg*s^-1*cm^-2',description='Zeropoint: %i Jy, Eff_wave: %f' %(zp[name],wave[name]))

        
    return [Jy_col,F_lam,lam_F_lam]



