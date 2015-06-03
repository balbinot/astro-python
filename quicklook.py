#! /usr/bin/env python
#mod date: 2015-02-26
import argparse
import pyfits
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d


def plot_spectrum(spectrum, wave, header, color='k',title=None,xlim=None,ylim=None):
 
  plt.plot(wave,spectrum, ls='-', lw=1.5, color=color)#, label=header['object'])
  if xlim is not None:
    plt.xlim(xlim)
  #else:
  #  plt.xlim(min(wave),max(wave))
  if ylim is not None:
    plt.ylim(ylim)
  #else:
  #  plt.ylim(min(spectrum), max(spectrum))
    
  plt.xlabel('$\lambda$ ($\AA$)', fontsize=30)
  plt.tick_params(axis='x', which='both', labelsize=22)
  plt.tick_params(axis='y', which='both', left='off',right='off', labelleft='off')

  #plt.legend(frameon=False, loc='upper right', markerscale=None, fontsize=22)
  plt.gca().xaxis.labelpad = 10
  if title:
    plt.title(title)


def get_spectrum(filename):
  '''Read spectrum using pyfits.  Returns spectrum and header.
  '''
  spectrum, header = pyfits.getdata(filename, 0, header=True)

  if len(spectrum) != int(header['naxis1']):
    spectrum = spectrum[0]

  
  lambda1 = header['crval1']
  lambda2 = float(header['naxis1'])*header['cdelt1']+float(header['crval1'])

  wave = np.arange(lambda1, lambda2, header['cdelt1'])
  wave = wave[:len(spectrum)]
 
  return (header, wave, spectrum)


def onkey(event):
  """Allows the user to scroll through a list of spectra using the left and
     right keys.  If the current index is equal to the end of the list, it 
     simply sets the index back to 0 to allow for endless scrolling fun.  Press
     q to quit out.  It's just like IRAF's splot, just not horrible.
  """

  global current
  global s

  if event.key == 'right':
    fig.clf()
    if current == len(fileslist[:-1]): 
      current = 0
    else: 
      current += 1
    header, wave, spectrum = get_spectrum(fileslist[current])
    if s != 0:
      spectrum = smooth(spectrum, kernel=s)
    plot_spectrum(spectrum, wave, header,title=fileslist[current])
    fig.canvas.draw()


  if event.key == 'left':
    fig.clf()
    if current == 0: 
      current = len(fileslist[:-1])
    else: 
      current -= 1
    header, wave, spectrum = get_spectrum(fileslist[current])

    if s != 0:
      spectrum = smooth(spectrum, kernel=s)
    plot_spectrum(spectrum, wave, header,title=fileslist[current])
    fig.canvas.draw()

  if event.key == 'x':
    s += 0.5
    xlim = plt.gca().xaxis.get_view_interval()
    ylim = plt.gca().yaxis.get_view_interval()
    fig.clf()
    header, wave, spectrum = get_spectrum(fileslist[current])
    if s != 0:
      spectrum = smooth(spectrum, kernel=s)
    plot_spectrum(spectrum, wave, header,title=fileslist[current],xlim=xlim,ylim=ylim)
    fig.canvas.draw()

  if event.key == 'z':
    if s == 0:
      return
    s -= 0.5
    xlim = plt.gca().xaxis.get_view_interval()
    ylim = plt.gca().yaxis.get_view_interval()
    fig.clf()
    header, wave, spectrum = get_spectrum(fileslist[current])
    if s != 0:
      spectrum = smooth(spectrum, kernel=s)
    plot_spectrum(spectrum, wave, header,title=fileslist[current],xlim=xlim,ylim=ylim)
    fig.canvas.draw()

  if event.key == 'r':
    fig.clf()
    header, wave, spectrum = get_spectrum(fileslist[current])
    spectrum = smooth(spectrum, kernel=s)
    plot_spectrum(spectrum, wave, header,title=fileslist[current])
    fig.canvas.draw()

    
 
  if event.key == 'q': exit()


def smooth(spectrum, kernel=0):
  """Smooths the input spectrum using a user-specified Gaussian kernel.
  """

  spectrum = gaussian_filter1d(spectrum, sigma=kernel)

  return spectrum


parser = argparse.ArgumentParser(description='An easy way to view a 1D, single extension, spectrum.  Will not work on .ms.fits files.')
parser.add_argument('input', metavar='input', nargs='+', type=str, help='Filename of single spectrum or list of spectra.')
parser.add_argument('--s', type=float, default=0, help='Kernel size (in pix) to use for Gaussian smoothing.')
parser.add_argument('--a',nargs=1,type=int,help='Specify starting fileno')

args = parser.parse_args()
fileslist = args.input

fig, ax = plt.subplots(figsize=(15,7), dpi=72) 
fig.subplots_adjust(wspace=0.25, left=0.05, right=0.95,
                    bottom=0.125, top=0.9)


cid = fig.canvas.mpl_connect('key_press_event', onkey)  #allows key presses to be read within the plotting window.

if args.a:
  split = np.array([int(x.split('.')[0]) for x in fileslist])
  current = np.where(split == args.a[0])
  current = current[0][0]
  print fileslist[current]
  header0, wave0, spectrum0 = get_spectrum(fileslist[current])

else:
  header0, wave0, spectrum0 = get_spectrum(fileslist[0])
  current = 0

if args.s != 0:
  spectrum0 = smooth(spectrum0, kernel=args.s)
s = args.s
  

plot_spectrum(spectrum0, wave0, header0,title=fileslist[current])
plt.show()
