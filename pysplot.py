#! /usr/bin/env python
import numpy as np
import pyfits
import argparse
import ds9
import matplotlib.pyplot as plt

current = 0

def onkey(event):
    global current
    global fig
    global rows
    global header
    global keys
    
    if event.key == 'right':
        if current >= len(rows)-1:
            return
        fig.clf()
        current = current + 1
        plt.plot(rows[current])
        if keys:
            plt.title(header[keys[current]])
        fig.canvas.draw()
        #print current

    elif event.key == 'left':
        if current == 0:
            return
        fig.clf()
        current = current - 1
        plt.plot(rows[current])
        if keys:
            plt.title(header[keys[current]])
        fig.canvas.draw()
        #print current

    elif event.key == 'q':
        plt.close(fig)



parser = argparse.ArgumentParser(description='Python version of splot.  [Left] and [Right] to scroll through spectra.  [Q] to quit.')
    
parser.add_argument('file',metavar='file',help='File with 1D spectra.')

args = parser.parse_args()

f = pyfits.open(args.file)

data = f[0].data
header = f[0].header

keys = [x for x in header.keys() if 'APID' in x]

if data.ndim <= 1:  # 1D array
    rows = data
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

