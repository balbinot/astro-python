#! /usr/bin/env python
import numpy as np
import pyfits
import argparse
#import ds9
import matplotlib.pyplot as plt

current = 0

def get_spectrum(filename,getheader=False):
    """Read spectrum using pyfits.  Returns spectrum and header.
    """
    spectrum, header = pyfits.getdata(filename, 0, header=True)
    #spectrum = spectrum[:1][0]
    
    if type(spectrum) is not 'list': spectrum = list(spectrum)
    if len(spectrum) != int(header['naxis1']): spectrum = spectrum[0]
    
    
    lambda1 = header['crval1']
    lambda2 = np.round(len(spectrum)*header['cdelt1']+header['crval1'])

    wave = np.arange(lambda1, lambda2, header['cdelt1'])
    wave = wave[:len(spectrum)]
    
    if getheader:
        return (wave,spectrum,header)

    else:
        return (wave, spectrum)
  




def onkey(event):
    global current
    global fig
    global xs
    global ys
    global header
    global keys
    
    if event.key == 'right':
        if current >= len(xs)-1:
            return
        fig.clf()
        current = current + 1
        plt.plot(xs[current],ys[current])
        if keys:
            plt.title(keys[current])
        fig.canvas.draw()
        #print current

    elif event.key == 'left':
        if current == 0:
            return
        fig.clf()
        current = current - 1
        plt.plot(xs[current],ys[current])
        if keys:
            plt.title(keys[current])
        fig.canvas.draw()
        #print current

    elif event.key == 'q':
        plt.close(fig)



parser = argparse.ArgumentParser(description='Python version of splot.  [Left] and [Right] to scroll through spectra.  [Q] to quit.')
    
parser.add_argument('filelist',nargs='+',help='File(s) with 1D spectra.')

args = parser.parse_args()

keys = args.filelist

data = [get_spectrum(f) for f in keys]

xs, ys = zip(*data)

fig = plt.figure()
cid = fig.canvas.mpl_connect('key_press_event', onkey)
plt.plot(xs[0],ys[0])
plt.show()

'''
fig = plt.figure()
cid = fig.canvas.mpl_connect('key_press_event', onkey)
    plt.plot(data)
    plt.show()
    exit()
else:
    rows = [x[:] for x in data[:]]
    
fig = plt.figure()
cid = fig.canvas.mpl_connect('key_press_event', onkey)

plt.plot(rows[0])

if keys:
    plt.title(header[keys[0]])
plt.show()
'''
