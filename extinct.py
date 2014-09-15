#! /usr/bin/env python
import argparse
from astropy.table import Table,Column,join
import numpy as np
import os

def get_intrinsic(specTable,refDict):
    intrinsic = []
    for row in specTable:
        if row and row['Lum'] != 'V':
            if row['SpTy'][0:2] == 'F2':
                intrinsic.append(np.mean([refDict['F0'],refDict['F5']]))
            elif row['SpTy'][0:2] == 'A7':
                intrinsic.append(np.mean([refDict['A8'],refDict['A5']]))
            elif row['SpTy'][0:2] == 'A0':
                intrinsic.append(np.mean([refDict['B9'],refDict['A2']]))
            elif row['SpTy'][0:2] == 'B2':
                intrinsic.append(np.mean([refDict['B1'],refDict['B3']]))
            else:
                intrinsic.append(refDict.get(row['SpTy'][0:2]))
        else:
            intrinsic.append(None)
    return intrinsic


def get_observed(specTable,photDict):
    observed = []
    for row in specTable:
        observed.append(photDict.get(row['oID']))

    return observed
    

def get_Av(observed,intrinsic,R=3.1):
    EBV = [BV-BV0 if BV and BV0 else None for BV,BV0 in zip(observed,intrinsic) ]

    Av = [x*R if x else 0.0 for x in EBV]

    return Av

def main():
    parser = argparse.ArgumentParser(description='Calculate Av via observed color')
    parser.add_argument('specfile',help='TSV file with observed spectral types.')
    parser.add_argument('reffile',help='Reference file with intrinsic colors.')
    parser.add_argument('photfile',help='File with photometry')
    ##reffile = ~/Software/python/colorExcess/supergiants.tsv
    parser.add_argument('-R',type=float,default=3.1,help='Rv. Default=3.1')

    args = parser.parse_args()

    specTable = Table.read(args.specfile,format='ascii.tab')
    refTable = Table.read(args.reffile,format='ascii.tab')
    photTable = Table.read(args.photfile)

    # Make dictionary of spectral type to color
    refDict = {MK[0:2]:np.float(color) for MK,color in zip(refTable['MK'],refTable['B-V']) if MK}

    # Make dictionary of oID to observed color
    photDict = {oID:np.float(x-y) for oID,x,y in zip(photTable['oID'],photTable['B'],photTable['V'])}
    
    intrinsic = get_intrinsic(specTable,refDict)
    observed = get_observed(specTable,photDict)

    Av = get_Av(observed,intrinsic,args.R)
    #for x,y,z in zip(observed,intrinsic,Av,photTable):
    #    print x,y,z
    #exit()
    
    # update input list
    specTable['SpTy'] = [x if x else '' for x in specTable['SpTy']]
    specTable['Lum'] = [x if x else '' for x in specTable['Lum']]
    specTable['Notes'] = [x if x else '' for x in specTable['Notes']]
    updateTable = zip(specTable['oID'],specTable['File'],specTable['SpTy'],specTable['Lum'],Av,specTable['Notes'])
    updateTable = Table(rows=updateTable,names=('oID','SpecFile','SpTy','Lum','Av','Notes'))


    # join tables
    photTable = join(photTable,updateTable,join_type='left')
    outFITS = '../Drout_list_ZOMG.fits'
    outTSV = os.path.splitext(outFITS)[0] + '.tsv'

    #photTable.write(outFITS)
    photTable.write(outTSV,format='ascii.tab')
    
    


if __name__ == '__main__':
    main()
