#! /usr/bin/env python
import argparse
import matplotlib.pyplot as plt
import matplotlib.patheffects as PE
from scipy.ndimage.interpolation import shift
import pyfits
import os.path
import numpy as np

class Plotter(object):
    def __init__(self,filelist,step,outdir):
        self.filelist = filelist
        self.step = step
        self.outdir = outdir
        self.current = 0
        self.orig_data = pyfits.getdata(filelist[0])
        self.active_data = self.orig_data.copy()
        # store total offsets
        self.offsets = np.zeros((len(filelist),2))
        

        # set up subparser
        self.subparser=argparse.ArgumentParser(description='Parse window text')
        self.subparser.add_argument('-s',type=float,help='Step size, default=%.1f' % self.step)
        self.subparser.add_argument('-q',action='store_true',help='Quit')
        self.subparser.add_argument('-r',action='store_true',help='Restore original')
        self.fig = plt.figure()
        self.fig.canvas.mpl_disconnect(self.fig.canvas.manager.key_press_handler_id)
        self.keycid = self.fig.canvas.mpl_connect('key_press_event',self.onkey)
        self.pausetext = '-'
        self.pid = None

        #Show initial
        self.display(self.active_data)
        self.displaytext('[0, 0], s=%.2f'%self.step,x=0.60)
        plt.show()
        

    def display(self, data):
        self.fig.clf()
        plt.imshow(data,origin='lower',interpolation='none')
        plt.title(self.filelist[self.current])
        self.fig.canvas.draw()

    def displaytext(self,text,x=0.05,y=0.05,remove=None):
        if remove:
            remove.remove()
        pid = plt.text(x,y,text,color='w',
                       horizontalalignment='left',
                       verticalalignment='bottom',
                       transform=plt.gca().transAxes,
                       path_effects=[PE.withStroke(linewidth=2,foreground='k')])
        self.fig.canvas.draw()
        return pid
        


    def parsetext(self,text):
        args = None
        try:
            args = self.subparser.parse_args(text.split())
        except SystemExit:
            return
            
        if not args:
            return
            
        if args.q:
            plt.close()
            exit()

        if args.s:
            self.step = args.s

        if args.r:
            self.active_data = self.orig_data.copy()
            self.offsets[self.current][0] = 0.0
            self.offsets[self.current][1] = 0.0

        
    def pausekey(self,event):
        if event.key == 'enter':
            self.fig.canvas.mpl_disconnect(self.keycid)
            self.keycid = self.fig.canvas.mpl_connect('key_press_event',self.onkey)
            self.parsetext(self.pausetext)
            self.pausetext = '-'
            self.display(self.active_data)
            self.displaytext('[%.2f, %.2f] s=%.2f'%
                             (self.offsets[self.current][0],
                              self.offsets[self.current][1],
                              self.step),
                             x=0.60)
            return

        elif event.key == 'backspace':
            self.pausetext = self.pausetext[0:-1]
            
        elif len(event.key) > 1:
            return
            
        else:
            self.pausetext = ''.join([self.pausetext,event.key])

        self.pid = self.displaytext(self.pausetext,remove=self.pid)
        return
        
    
    def onkey(self, event):
        if event.key in ['.','>']:
            if self.current >= len(self.filelist)-1:
                return
            self.current += 1
            self.orig_data = pyfits.getdata(self.filelist[self.current])
            self.active_data = self.orig_data.copy()
            
        elif event.key in [',','<']:
            if self.current == 0:
                return
            self.current -= 1
            self.orig_data = pyfits.getdata(self.filelist[self.current])
            self.active_data = self.orig_data.copy()

        elif event.key == '-':
             self.fig.canvas.mpl_disconnect(self.keycid)
             self.pausetext = '-'
             self.pid = self.displaytext(self.pausetext)
             self.keycid = self.fig.canvas.mpl_connect('key_press_event',self.pausekey)
             return

        elif event.key == 'left':
            if self.active_data is None:
                return
            self.active_data = shift(self.active_data,[0,-self.step])
            self.offsets[self.current][0] -= self.step

        elif event.key == 'right':
            if self.active_data is None:
                return
            self.active_data = shift(self.active_data,[0,self.step])
            self.offsets[self.current][0] += self.step

        elif event.key == 'down':
            if self.active_data is None:
                return
            self.active_data = shift(self.active_data,[-self.step,0])
            self.offsets[self.current][1] -= self.step

        elif event.key == 'up':
            if self.active_data is None:
                return
            self.active_data = shift(self.active_data,[self.step,0])
            self.offsets[self.current][1] += self.step

        self.display(self.active_data)
        self.displaytext('[%.2f, %.2f] s=%.2f'%
                         (self.offsets[self.current][0],
                          self.offsets[self.current][1],
                          self.step),
                         x=0.60)
        

def main():
    parser = argparse.ArgumentParser(description='Nudge an image, or series of images, by a given step size')
    parser.add_argument('filelist',nargs='+',help='List of input FITS files to be nudged.')
    parser.add_argument('-s',type=float,default=1.0,help='Specify initial step size (default=1.0).')
    parser.add_argument('-o',type=str,default='.',help="Output directory (default='.').")

    args = parser.parse_args()
    plotter = Plotter(args.filelist,args.s,args.o)


if __name__ == '__main__':
    main()







