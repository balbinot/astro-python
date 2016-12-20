#! /usr/bin/env python
import argparse
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u

class HRD(object):
    def __init__(self,table,tcol='Teff',lcol='Lum',xlim=None,ylim=None,ax=None):
        self.table = table
        self.tcol = tcol
        self.lcol = lcol
        self.xlim = xlim
        self.ylim = ylim

        # set ax to false to suppress
        if ax != False:
            self._make_figure(xlim,ylim,ax=ax)

    def _make_figure(self,xlim=None,ylim=None,ms=30,ax=None,fontsize=12,nformatter=None):
        if ax:
            self.ax = ax
        else:
            self.ax = plt.figure(dpi=100).gca()

        self.ax.scatter(self.table[self.tcol],self.table[self.lcol],s=ms)
        if xlim:
            self.ax.set_xlim(xlim)
        elif self.xlim:
            self.ax.set_xlim(self.xlim)
        else:
            self.ax.invert_xaxis()

        if ylim:
            self.ax.set_ylim(ylim)
        elif self.ylim:
            self.ax.set_ylim(self.ylim)


        self.ax.set_xlabel(r'$T_{eff}\,\left(K\right)$',fontsize=fontsize)
        self.ax.set_ylabel(r'$\log\,\left(L/L_{\odot}\right)$',fontsize=fontsize)
        if nformatter:
            for label in self.ax.get_xticklabels()[::nformatter]:
                label.set_visible(False)

        self.ax.tick_params(labelsize=8)

    def show(self):
        plt.show()

    def save(self,filename):
        plt.savefig(filename)
        plt.close(plt.gcf())
            
        
    


def main():
    parser = argparse.ArgumentParser(description='Plot stars on HR diagram')
    parser.add_argument('table',type=str,help='Input table')
    parser.add_argument('-o',type=str,metavar='outfile',help='Save fig')
    parser.add_argument('-fmt',default='ascii.tab',help='Specify table format for astropy (default = ascii.tab')
    parser.add_argument('-tcol',type=str,default='Teff',help='Temp col name in table (default = Teff)')
    parser.add_argument('-lcol',type=str,default='Lum',help='Lum col name in table (default = Lum)')
    parser.add_argument('-xlim',nargs=2,type=float,help='Temp limits')
    parser.add_argument('-ylim',nargs=2,type=float,help='Lum limits')

    args = parser.parse_args()

    table = Table.read(args.table,format=args.fmt)
    hrd = HRD(table,args.tcol,args.lcol,args.xlim,args.ylim)

    labels = table['ID']
    xoffsets = [-20.]*len(labels)
    yoffsets = [0.]*len(labels)
    if labels[0][0] == 'F':
        xoffsets[-1] = -100

    for i,v in enumerate(zip(labels,xoffsets,yoffsets)):
        l,xo,yo = v
        x,y = table[args.tcol][i],table[args.lcol][i]
        if l == 'F15':
            x = 4450
        if l == 'F08':
            yo = 0.02
        if l == 'F07':
            yo = -0.02
        if l == 'S26':
            yo = -0.05
        if l == 'S06':
            yo = -0.02


        plt.annotate(l,xy=(x,y),xytext=(x+xo,y+yo),arrowprops=dict(arrowstyle="-"),fontsize=12,va='center')
        

    if args.o:
        hrd.save(args.o)
    else:
        hrd.show()



if __name__ == '__main__':
    main()
