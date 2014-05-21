#! /usr/bin/env python
import numpy as np
#import pyfits
import argparse
#import ds9
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
    
parser.add_argument('file',metavar='file',nargs='*',help='File(s) with 1D spectra.')

args = parser.parse_args()

data = []
for f in args.file:
    data.append([np.genfromtxt(f,names=['x','y']),f])

fig = plt.figure()
cid = fig.canvas.mpl_connect('key_press_event', onkey)
plt.plot([x[0]['x'] for x in data],[x[0]['y'] for x in data])
plt.title([x[1] for x in data])
plt.show()
exit()
#else:
#    rows = [x[:] for x in data[:]]
    
fig = plt.figure()
cid = fig.canvas.mpl_connect('key_press_event', onkey)

plt.plot(rows[0])

if keys:
    plt.title(header[keys[0]])
plt.show()

