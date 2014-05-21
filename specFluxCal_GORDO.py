#! /usr/bin/env python
broken_mods = []
import argparse
import sys

try: import numpy as np
except: broken_mods.append('numpy')
try: import matplotlib.pyplot as plt
except: broken_mods.append('matplotlib')
try: import pyfits
except: broken_mods.append('pyfits')
try: from fitting import *
except: broken_mods.append('fitting.py')
from fitting import *

if broken_mods:
  errM = 'The following mods need to be installed to run specFluxCal.py:'
  broken_mods = '\n'.join([line for line in broken_mods])
  sys.exit('%s\n%s' % (errM, broken_mods))

Zp = 48.60
speed_of_light = 2.99792458e18  


def get_extinction(filename):
  f = open(filename, 'r+')
  
  w = []
  ext = []
  for line in f:
    wc, extc = line.split()
    w.append(np.float(wc))
    ext.append(np.float(extc))

  f.close()
  return (w, ext)

def get_spectrum(filename):
  """Read spectrum using pyfits.  Returns spectrum and header.
  """
  spectrum, header = pyfits.getdata(filename, 0, header=True)
  spectrum = spectrum[:1][0]

  if type(spectrum) is not 'list': spectrum = list(spectrum)
  if len(spectrum) != int(header['naxis1']): spectrum = spectrum[0]

  
  lambda1 = header['crval1']
  lambda2 = np.round(len(spectrum)*header['cdelt1']+header['crval1'])

  wave = np.arange(lambda1, lambda2, header['cdelt1'])
  wave = wave[:len(spectrum)]
 
  return (header, wave, spectrum)


def get_std_data(filename):
  """Read in the calibrated standard star data.
  """  
  f = open(filename, 'r+')

  w = []
  m = []
  bin = []
  
  for line in f:
    wc, mc, binc = line.split()
    w.append(np.float(wc))
    m.append(np.float(mc))
    bin.append(np.float(binc))

  f.close()
  return (w, m, bin)

def get_ABmag(wavelength, flux):

  ABmag = -2.5*np.log10(flux*wavelength**2/speed_of_light) - Zp

  return ABmag

def get_flux(wavelength, magnitudes):
  """Convert magnitudes to fluxes in ergs/s/cm**2/Hz.
  """

  flux = np.power(10,(-Zp - magnitudes)/2.5) *speed_of_light / wavelength**2.

  return flux


def get_closest_wave(cal_wave, wavelist):
  """Find the closest wavelength in an observed spectrum to the calibrated 
     standard spectrum.
  """
  diff = np.abs(cal_wave-wavelist)
  
  return wavelist[np.argmin(diff)]


def integrate(start, finish, wavelengths, spectrum):
  """Integrate the values of spectrum between the start and finish wavelengths.
     The wavelengths argument contains the corresponding wavelength values for 
     the input spectrum. 
  """  

  sumlist = [y for x,y in zip(wavelengths, spectrum) if x >= start and x < finish]
  
  return np.sum(sumlist)/(finish-start)


def exptime_normalize(spectrum, header):
  """Normalize the observed spectrum by its exposure time which is fetched from
     the header card header['exptime'].
  """

  spectrum = np.array(spectrum)

  return spectrum/header['exptime']


def get_summed_counts(cal_wave, cal_mag, cal_bin, std_wave, std_spectrum):
  """Takes the calibrated and observed spectrum values as arguments.  
     The observed spectrum is then summed over the same bins and bin sizes as
     the calibrated standard spectrum.  
     The left wavelength edge, the integrated counts at that wavelength, and the
     corresponding calibrated spectrum mags at that wavelength are returned.
  """
  
  data = []

  for w,b,m in zip(cal_wave, cal_bin, cal_mag):
    start_wave = get_closest_wave(w, std_wave)
    end_wave = get_closest_wave(w+b, std_wave)

    if start_wave != end_wave: 
      counts = integrate(start_wave, end_wave, std_wave, std_spectrum)
      data.append((start_wave, counts, m))

  edge_wave, edge_counts, edge_mags = zip(*data)
  
  return (np.array(edge_wave),np.array(edge_counts),np.array(edge_mags))


def get_response(counts, cal_flux):

  resp_func = cal_flux / counts

  return resp_func


def fit_response(resp_func, edge_wave, wave, weights=None, *args, **kwargs):
  
  f = interp1d(edge_wave, resp_func, kind='cubic')
  
  interp_resp_func = np.zeros(wave.shape)
  interp_resp_func[:-1] = f(wave[:-1])[:]
  interp_resp_func[-1] = interp_resp_func[-1]
 
  return interp_resp_func


def flux_cal(spectrum, resp_func):
  return spectrum * resp_func


def airmass_correction(mags, airmass, extinction):
  return mags - airmass*extinction


def smooth(spectrum, kernel):
  return medfilt(spectrum, kernel_size = kernel)

def main():

  parser = argparse.ArgumentParser(description='Extract the one dimensional spectra from mult-extension fits files.')
  parser.add_argument('spec',metavar='spectrum', type=str, 
                      help='Filename for single spectrum.')
  parser.add_argument('obstd',metavar='standard', type=str, 
                      help='Filename for standard spectrum.') 
  parser.add_argument('caldat',metavar='calibration-data', type=str, 
                      help='Filename for standard star mags.')
  parser.add_argument('--l', dest='list',action='store_true', 
                      help='List of sources to flux calibrate.')
  parser.add_argument('extinction',metavar='extinction',
                      help='Filename of extinction coefficients.')
  args = parser.parse_args()
  
  # Read data from fits files
  header, wave, spectrum = get_spectrum(args.spec)
  cal_wave, cal_mag, cal_bin = get_std_data(args.caldat)
  std_header, std_wave, std_spectrum = get_spectrum(args.obstd)

  # Read extinction
  ext_wave, extinction = get_extinction(args.extinction)

  # Normalize by exposure time
  std_spectrum = exptime_normalize(std_spectrum, std_header)
  spectrum = exptime_normalize(spectrum,header)

  # Bin standard star spectrum to match calibration bins
  edge_wave, edge_counts, edge_mags = get_summed_counts(cal_wave, cal_mag, cal_bin, std_wave, std_spectrum)

  # Convert calibration magnitudes to flux units
  cal_flux = get_flux(edge_wave, edge_mags)

  # Spline fit calibration flux  
  calfit, calresid = spline_fit(cal_flux,edge_wave,edge_wave,weights=None)

  # Divide calibrated spectrum into standard star spectrum.
  resp_func = get_response(edge_counts, calfit)
  plt.plot(edge_wave,resp_func)
  plt.show()
  exit()

  ####### OLD SHITTT
  '''
  # Fit calibrated spectrum to blackbody. Residuals will then be used as
  #  weights for fitting the response function.
  bbfit, bbresiduals = blackbody_fit(cal_flux, edge_wave, edge_wave)
  
  # Divide calibrated spectrum into standard star spectrum.
  resp_func = get_response(edge_counts, cal_flux)
  #resp_func = get_response(edge_counts, bbfit)

  # Interpolate and fit response function to match wavelengths of target spectrum
  interp_resp_func, interp_residuals = chebyshev_fit(resp_func, edge_wave, std_wave, weights=1/bbresiduals,order=7)
  #interp_resp_func = fit_response(resp_func, edge_wave, std_wave, weights=None)
  
  # Multiply standard star spectrum by response (sanity check)
  std_spectrum_flux = flux_cal(std_spectrum, interp_resp_func)

  # Convert standard star to AB mag  
  std_ABmag = get_ABmag(std_wave, std_spectrum_flux)

  # Fit curve to extinction and interpolate to target wavelengths
  interp_ext, ext_residuals = legendre_fit(extinction, ext_wave, std_wave,weights=None,order=10)

  # Airmass correction
  std_airmass = np.float(std_header['airmass'])
  std_ABmag_corrected = airmass_correction(std_ABmag, std_airmass, interp_ext)
  '''


  
  '''
  # Wait, let's try this again from the beginning
  # C = 2.5 log10(O/(TBF)) + A*E
  #   Need E interp'd to edge_wave
  E,er = legendre_fit(extinction, ext_wave, edge_wave,weights=None, order=12)

  C =  2.5 * np.log10(edge_counts/bbfit) + std_airmass* E

  # Convert back to flux
  C_resp_func = get_flux(edge_wave,C)

  # Interpolate resp_func to std waves
  interp_resp_func, interp_residuals = chebyshev_fit(resp_func, edge_wave, std_wave, weights=1/bbresiduals,order=7)

  # Multiply standard star spectrum by response (sanity check)
  std_spectrum_flux = flux_cal(std_spectrum, interp_resp_func)

  # Convert standard star to AB mag  
  std_ABmag = get_ABmag(std_wave, std_spectrum_flux)

  # Fit curve to extinction and interpolate to target wavelengths
  interp_ext, ext_residuals = legendre_fit(extinction, ext_wave, std_wave,weights=None,order=12)

  # Airmass correction
  std_airmass = np.float(std_header['airmass'])
  std_ABmag_corrected = airmass_correction(std_ABmag, std_airmass, interp_ext)
  '''
  
  '''
  plt.plot(std_wave, std_ABmag_corrected, 'b')
  plt.plot(cal_wave, cal_mag, 'r.')
  plt.gca().invert_yaxis()
  plt.xlim(3500, 6000)
  plt.show()
  #f = open('spec_Feige-66.specFluxCal.dat','w')
  #for x,y in zip(std_wave,std_ABmag_corrected):
  #  f.write('{0}\t{1}\n'.format(x,y))
  #f.close()
  '''


  
if __name__ == '__main__':
  main()

