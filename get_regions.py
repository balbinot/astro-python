#! /usr/bin/env python
import pyregion
import argparse
import ds9


def main():
    parser = argparse.ArgumentParser(description='Generate regions from ds9')
    parser.add_argument('-o',dest='outfile',type=str, default=None,help='Output file.')

    args = parser.parse_args()
    
    d = ds9.ds9()

    frames = d.get('frame all').split()

    files = []
    coords = []
    print 'Loading regions from %i frames' % len(frames)
    for frame in frames:
        d.set('frame %s' % frame)
        regions = d.get('regions')
        try:
            r = pyregion.parse(regions)[0].coord_list
        except:
            continue
        files.append(d.get('file'))
        coords.append((r[0],r[1],r[2]))

    if len(coords) == 0:
        print 'No regions found'
        return 1

    if args.outfile == None:
        for fo,coo in zip(files,coords):
            print '%s\t%f\t%f\t%f' % (fo,coo[0],coo[1],coo[2])
        
    else:
        f = open(args.outfile,'w')
        for fo,coo in zip(files,coords):
            f.write('%s\t%f\t%f\t%f\n' % (fo,coo[0],coo[1],coo[2]))

        f.close()
    return 0


if __name__ == '__main__':
    main()
