#! /usr/bin/env python
#version 12/19/2016
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table,hstack,Column
from astropy.coordinates import SkyCoord
import astropy.units as u
from photutils.centroids import centroid_com,centroid_1dg,centroid_2dg
from photutils.morphology import data_properties
from photutils.utils import cutout_footprint
from photutils import properties_table, CircularAperture, CircularAnnulus, aperture_photometry, SkyCircularAnnulus,SkyCircularAperture
import numpy as np

def is_pix_coord(coord):
    try:
        x = np.float(coord[0])
        y = np.float(coord[1])
    except ValueError:
        return False
    return True

def is_sky_coord(coord):
    try:
        coo = SkyCoord(' '.join(coord),unit=(u.hourangle,u.deg),frame='icrs')
    except ValueError:
        return False
    return True

def parse_wcs(header):
    try:
        wcs = WCS(header,naxis=2)
    except:
        print 'Failed to parse wcs'
        return None
    return wcs

def pix2wcs(coord,wcs):
    ncoord = wcs.all_pix2world([coord],1,ra_dec_order=True)
    ncoord = SkyCoord(ncoord,unit=(u.deg,u.deg),frame='icrs') #ra,dec
    ra,dec = ncoord.to_string('hmsdms',sep=':')[0].split()

    return ra,dec

def wcs2pix(coord,wcs):
    if not isinstance(coord,SkyCoord):
        coord = SkyCoord(' '.join(coord),unit=(u.hourangle,u.deg),frame='icrs')
    ncoord = wcs.all_world2pix(coord.ra,coord.dec,1,ra_dec_order=True)
    x,y = ncoord[0],ncoord[1]
    return x,y


def centroid(data,coord,box_size=30,wcs=None):
    if is_pix_coord(coord):
        coord = [np.float(coord[0]),np.float(coord[1])]
        print 'Centroiding at (%.2f,%.2f)' % (coord[0],coord[1]),
    elif is_sky_coord(coord):
        if wcs:
            try:
                x,y = wcs2pix(coord,wcs)
                print 'Centroiding at (%s,%s) -> (%.2f,%.2f)' % (coord[0],coord[1],x,y)
                coord = (x,y)
            except:
                raise ValueError('Cannot parse coordinates with input wcs')
        else:
            raise InputError("Must provide wcs object to parse sky coordinates")
    else:
        raise ValueError('Cannot parse coordinates')       


    ## make cutout mask
    #print coord
    _,_,_,slices = cutout_footprint(data,coord,box_size=box_size)
    mask = np.ones_like(data,dtype=bool)
    mask[slices] = False    
    xc,yc = centroid_2dg(data,mask=mask)
    print ' -> (%f, %f)' % (xc,yc),

    if wcs:
        ra,dec = pix2wcs([xc,yc],wcs)
        print ' [%s, %s]' % (ra,dec)
        tbl = Table([[xc],[yc],[ra],[dec]],names=['xcen','ycen','ra','dec'])
        return tbl
    else:
        print
        tbl = Table([[xc],[yc]],names=['xcen','ycen'])
        return tbl

    
def qphot(data,coord,rad,skyradin,skyradout,wcs=None,calfctr=None,skycoord=None,unit='Jy',error=None,filter=None):
    if is_pix_coord(coord):
        coord = [np.float(coord[0]),np.float(coord[1])]
        aperture = CircularAperture(coord,r=rad)
    else:
        coord = SkyCoord(' '.join(coord),unit=(u.hourangle,u.deg),frame='icrs')
        aperture = SkyCircularAperture(coord,r=rad*u.arcsec)

    if skycoord:
        if is_pix_coord(skycoord):
            scoord = [np.float(skycoord[0]),np.float(skycoord[1])]
            annulus = CircularAnnulus(scoord,skyradin,skyradout)
        else:
            scoord = SkyCoord(' '.join(skycoord),unit=(u.hourangle,u.deg),frame='icrs')
            annulus = SkyCircularAnnulus(scoord,skyradin*u.arcsec,skyradout*u.arcsec)
    else:
        if isinstance(coord,SkyCoord):
            annulus = SkyCircularAnnulus(coord,skyradin*u.arcsec,skyradout*u.arcsec)
        else:
            annulus = CircularAnnulus(coord,skyradin,skyradout)

    #mask out nans or negative
    mask = np.logical_or(np.isnan(data),data<0.)
    
    apflux = aperture_photometry(data,aperture,wcs=wcs,error=error,mask=mask)
    skyflux = aperture_photometry(data,annulus,wcs=wcs,error=error,mask=mask)
    phot_table = hstack([apflux, skyflux], table_names=['src', 'sky'])

    #calculate mean local background in annulus
    if isinstance(annulus,SkyCircularAnnulus):
        sky_area = annulus.to_pixel(wcs).area()
    else:
        sky_area = annulus.area()

    if isinstance(aperture,SkyCircularAperture):
        src_area = aperture.to_pixel(wcs).area()
    else:
        src_area = aperture.area()
    sky_mean = phot_table['aperture_sum_sky'] / sky_area
    sky_sum = sky_mean * src_area
    
    final_sum = phot_table['aperture_sum_src'] - sky_sum
    phot_table['residual_aperture_sum'] = final_sum
    phot_table['residual_aperture_sum'].unit = unit
    phot_table['aperture_area_sky'] = sky_area
    phot_table['aperture_area_src'] = src_area
    phot_table['aperture_mean_sky'] = sky_mean
    phot_table['aperture_mean_sky'].unit = unit
    phot_table['aperture_rad_src'] = rad
    phot_table['aperture_irad_sky'] = skyradin
    phot_table['aperture_orad_sky'] = skyradout
    phot_table['aperture_area_sky'].unit = u.pix**2
    phot_table['aperture_area_src'].unit = u.pix**2

    if error is not None:
        src_err = phot_table['aperture_sum_err_src']
        sky_err = phot_table['aperture_sum_err_sky']
        src_var = src_err/phot_table['residual_aperture_sum']
        sky_var = sky_err/sky_sum
        color_err = 0. # 10 percent photoerr
        color_var = color_err * phot_table['residual_aperture_sum']
        phot_table['residual_aperture_err'] = np.sqrt(src_var**2+sky_var**2+color_var**2)
    
    phot_table.remove_columns(['xcenter_src','ycenter_src','xcenter_sky','ycenter_sky'])
    if 'center_input' in phot_table.colnames:
        phot_table.remove_columns('center_input')
    if 'center_input_src' in phot_table.colnames:
        phot_table.remove_columns('center_input_src')
    if 'center_input_sky' in phot_table.colnames:
        phot_table.remove_columns('center_input_sky')
    
    if isinstance(coord,SkyCoord):
        phot_table['xcen'],phot_table['ycen'] = wcs2pix(coord,wcs)
        phot_table['ra'],phot_table['dec'] = coord.to_string('hmsdms',sep=':').split()
    else:
        phot_table['xcen'] = coord[0]
        phot_table['ycen'] = coord[1]

    if skycoord:
        if isinstance(scoord,SkyCoord):
            phot_table['xcensky'],phot_table['ycensky'] = wcs2pix(scoord,wcs)
            phot_table['rasky'],phot_table['decsky'] = scoord.to_string('hmsdms',sep=':').split()
        else:
            phot_table['xcensky'] = scoord[0]
            phot_table['ycensky'] = scoord[1]

    if wcs:
        if 'ra' not in phot_table.colnames:
            ra,dec = pix2wcs(coord,wcs)
            phot_table['ra'] = ra
            phot_table['dec'] = dec

        if skycoord:
            if 'rasky' not in phot_table.colnames:
                skyra,skydec = pix2wcs(scoord,wcs)
                phot_table['rasky'] = skyra
                phot_table['decsky'] = skydec


    if calfctr:
        flux = phot_table['residual_aperture_sum']/calfctr
        phot_table.add_column(Column(flux,name='ap_flux',unit=unit))

        if error is not None:
            fluxerr = phot_table['residual_aperture_err']/calfctr
            phot_table.add_column(Column(fluxerr,name='ap_err',unit=unit))

    if filter:
        phot_table['filter'] = str(filter)

    return phot_table
