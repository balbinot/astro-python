#! /usr/bin/env python

failed_modules = []

import sys
import argparse

try: import pyfits
except: failed_modules.append('pyfits')

try: import numpy as np
except: failed_modules.append('numpy')

try: import matplotlib.pyplot as plt
except: failed_modules.append('matplotlib')

if failed_modules:
  for module in failed_modules:  
    print "'%s' module required." % module
  sys.exit('Dependencies required')


def get_spectrum(filename):
  '''Read spectrum using pyfits.  Returns spectrum and header.
  '''
  spectrum,header = pyfits.getdata(filename, 0, header=True)

  return (header,spectrum[:1][0])


def crop_spectrum(spectrum, npix):
  '''If called, will crop a spectrum to the specified number of pixels'
  '''
  spectrum = spectrum[:npix]

  return spectrum


def extract_spectrum(spectrum, header, append=False, extract=False, useastrom=False):
  '''Takes an input spectrum and header and will write it to a file.  If
     useastrom is specified, the filename will contain the RA and DEC,
     otherwise the object name is used.
  '''

  if useastrom is True:
    outputfile = 'spec_'+','.join([header['RA'], header['DEC']])+'_1d-counts.fits' 
  else: outputfile = 'spec_'+header['object']+'_1d-counts.fits'

  if extract is True and append is False: 
    pyfits.writeto(outputfile,
                   spectrum, header=header, clobber=True)
  if extract is True and append is True:
    pyfits.append('spec_all_1d-counts.fits',
                   spectrum, header=header, clobber=True)


def plot_spectrum(spectrum, header, show=False, useastrom=False):

  if type(spectrum) is not 'list': spectrum = list(spectrum) #make sure spectrum is a list

  if len(spectrum) != int(header['naxis1']): spectrum = spectrum[0]
  
  lambda1 = header['crval1']
  
  lambda2 = np.round(len(spectrum)*np.round(header['cdelt1'],2)+header['crval1'])

  wave = np.arange(lambda1, lambda2, np.round(header['cdelt1'],2))
  
  fig, ax = plt.subplots(figsize=(14,8), dpi=72)
  fig.subplots_adjust(wspace=0.25, left=0.1, right=0.95,
                    bottom=0.125, top=0.925)

  plt.plot(wave,spectrum, ls='-', color='r')
  plt.xlim(lambda1, lambda2)
  plt.ylim(min(spectrum), max(spectrum))
  plt.xlabel('Wavelength (Angstroms)', fontsize=24)
  plt.ylabel('Counts', fontsize=24)
  
  if useastrom is True: 
    plt.title(','.join([header['RA'], header['DEC']]), fontsize=30)
  else: 
    plt.title(header['object'], fontsize=30)
 
  if show is True: plt.show()
  elif show is False and useastrom is False: 
    plt.savefig('spec_'+header['object']+'_1d-counts.pdf', 
              dpi=None, facecolor='w', edgecolor='w',
              transparent=True, bbox_inches=None, pad_inches=0.1,
              frameon=None)
  elif show is False and useastrom is True: 
    plt.savefig('spec_'+','.join([header['RA'], header['DEC']])+'_1d-counts.pdf', 
              dpi=None, facecolor='w', edgecolor='w',
              transparent=True, bbox_inches=None, pad_inches=0.1,
              frameon=None)


def sum_spectra(spec_list):
  '''Input is a list of spectra to sum.  The first spectrum of the list will
     be read in and act as the reference spectrum.  All subsequent spectra
     will be added to the reference spectrum.  If the spectra are of different
     lengths, then all the spectra will be cropped to the shortest spectrum.
     
     A HISTORY card will be added to the header to denote the filenames of the
     summed spectra.
     The 'EXPTIME' card will reflect the summed exposure times of all the
     summed spectra.
     The NAXIS1 card will reflect the number of pixels of the summed spectrum,
     which will be the value of the shortest spectrum.
  '''

  header0, spec0 = get_spectrum(spec_list[0])
  spec0 = spec0[0]
  
  exptime = header0['exptime']

  for line in spec_list[1:]:
    header, spec = get_spectrum(line)

    spec = spec[0]

    if header0['naxis1'] > header['naxis1']:
      spec0 = crop_spectrum(spec0, int(header['naxis1']))
    elif header['naxis1'] > header0['naxis1']:
      spec = crop_spectrum(spec, int(header0['naxis1']))

    spec0 += spec

    exptime += header['exptime']
    
  header0.update('exptime', int(exptime))
  header0.update('naxis1', int(len(spec0)))
  header0.add_history('specExtract.py: Summed ' + 
                      ','.join([x for x in spec_list]))

  return (header0,np.array([spec0]))

def main():
  '''Functionality:
      -show: a single spectrum or list of spectra.
      -sum: given a list of spectra, sum the matches (note, the summed spectra 
       may be matched by their header OBJECT cards *or* their RA and DEC cards.
      -extract: extract spectrum to a new filename.  Uses either the OBJECT
       header card for the filename *or* the RA and DEC cards if --useastrom is
       specified.
  '''

  parser = argparse.ArgumentParser(description='Extract the one dimensional \
                                   spectra from mult-extension fits files.')
  parser.add_argument('spec',metavar='spectrum', type=str, 
                      help='Filename for single spectrum.') 
  parser.add_argument('--show', action='store_true', 
                      help='If specified, display spectrum to screen, \
                      else print to pdf.')
  parser.add_argument('--list', action='store_true', help='List of sources to\
                      run through specExtract.py.')
  parser.add_argument('--sum', action='store_true', help='If set, will iterate\
                      through list to find sources with matching identifiers \
                      and then sum the spectra.')
  parser.add_argument('--extract', action='store_true', 
                      help='If set, extract to 1-d fits file.')
  parser.add_argument('--useastrom', action='store_true', 
                      help='If set, match sources using astrometry rather than header object cards.')
  parser.add_argument('--parse', action='store_true', 
                      help='If set,user may scroll through a list of spectra using N(ext) and B(ack).')
  args = parser.parse_args()
  
  if args.list is True and args.parse is not True:
    
    listfile = open(args.spec, 'r+')

    filelist = [x.strip('\n') for x in listfile if x[0] != '#']
   
    listfile.close()

    if args.useastrom is True: 
      unique_ids = set([(pyfits.getval(x, 'RA'), pyfits.getval(x, 'DEC'))
                        for x in filelist])
    else: unique_ids = set([pyfits.getval(x, 'object') for x in filelist])
    
    for line in unique_ids:
     
      if args.useastrom is True:
         matches = [x for x in filelist 
                    if (pyfits.getval(x, 'RA'), pyfits.getval(x, 'DEC')) == line]
      else: 
        matches = [x for x in filelist if pyfits.getval(x, 'object') == line]
      
      if len(matches) > 1 and args.sum is True:
        print 'Summing '+ ','.join([x for x in matches])+'.'
        try:        
          header, summed_spec = sum_spectra(matches)
       
          extract_spectrum(summed_spec, header, extract=args.extract, useastrom=args.useastrom)
          extract_spectrum(summed_spec, header, extract=args.extract, append=True)
          plot_spectrum(summed_spec, header, show=args.show, useastrom=args.useastrom)
        except: 
          print 'Unable to sum' + ','.join([x for x in matches])+'.'
          continue
      
      else:
        try: 
          header, spectrum = get_spectrum(matches[0])
          extract_spectrum(spectrum, header, extract=args.extract, useastrom=args.useastrom)
          extract_spectrum(spectrum, header, extract=args.extract, append=True)
          plot_spectrum(spectrum, header, show=args.show, useastrom=args.useastrom)
        except: print 'Error getting '+matches[0]+'.'
  
  elif args.list is True and args.parse is True:
    
    listfile = open(args.spec, 'r+')

    filelist = [x.strip('\n') for x in listfile if x[0] != '#']
   
    listfile.close()
    todo = raw_input('Please enter N(ext), B(ack), Q(uit):')
    index = 0
    while todo.lower()[0] is not 'q':
       header, spectrum = get_spectrum(filelist[index])
       plot_spectrum(spectrum, header, show=True)
       todo = raw_input('Please enter N(ext), B(ack), Q(uit)')       
       if todo.lower()[0] is 'n': index =+ 1
       if todo.lower()[0] is 'b': index =+ -1
    
  else:
    try:
      header, spectrum = get_spectrum(args.spec)
      extract_spectrum(spectrum, header, extract=args.extract)
      plot_spectrum(spectrum, header, show=args.show)
    except: print 'Error getting '+args.spec+'.'

if __name__ == '__main__':
  main()
