#! /usr/bin/env python
import argparse
import numpy as np
import bokeh.plotting as blt

def bokeh_plot(data,outfile):

    table = np.genfromtxt(data,names=['x','y'])
    blt.output_file(outfile)

    blt.hold()
    blt.figure(tools="pan,wheel_zoom,box_zoom,reset,previewsave,select")
    blt.scatter(table['x'], table['y'])
    blt.show()

    

def main():
    parser = argparse.ArgumentParser(description='Output data to HTML plot')
    parser.add_argument('data',type=str,help='Input data file')
    parser.add_argument('-o',type=str,default='output.html',dest='outfile',help='Specify output file (default = "output.html")')

    args = parser.parse_args()

    bokeh_plot(args.data,args.outfile)




if __name__ == '__main__':
    main()
