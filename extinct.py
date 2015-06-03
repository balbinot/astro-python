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


def isValid(cell):
    try:
        cell = float(cell)
    except:
        return False
    if cell == 0:
        return True
    if not cell:
        return False
    if cell == '':
        return False
    if cell > 80:
        return False
    if cell == '--':
        return False

    return True

def get_Av(observed,intrinsic,R=3.1):
    EBV = [BV-BV0 if BV and BV0 else None for BV,BV0 in zip(observed,intrinsic) ]

    Av = [x*R if x else 0.0 for x in EBV]

    return Av

def main():
    parser = argparse.ArgumentParser(description='Calculate Av via observed color')
    parser.add_argument('specfile',help='TSV file with observed spectral types.')
    parser.add_argument('reffile',help='Reference file with intrinsic colors.')
    parser.add_argument('photfile',help='FITS file with photometry')
    ##reffile = ~/Software/python/colorExcess/supergiants.tsv
    parser.add_argument('-o',metavar='outfile',dest='outfile',default=False,help="Output table (Default=photfile).")
    parser.add_argument('-R',type=float,default=3.2,help='Rv. Default=3.2')

    args = parser.parse_args()

    specTable = Table.read(args.specfile,format='ascii.tab')
    refTable = Table.read(args.reffile,format='ascii.tab')


    if os.path.splitext(args.photfile)[1] == '.tsv':
        photTable = Table.read(args.photfile,format='ascii.tab')
    else:
        photTable = Table.read(args.photfile)

    # Make dictionary of spectral type to color
    refDict = {MK[0:2]:np.float(color) for MK,color in zip(refTable['MK'],refTable['B-V']) if MK}

    # Make dictionary of oID to observed color
    BV = [float(x)-float(y) if isValid(x) and isValid(y) else None for x,y in zip(photTable['B'],photTable['V'])]
    photDict = {oID:color for oID,color in zip(photTable['oID'],BV)}
    
    intrinsic = get_intrinsic(specTable,refDict)
    observed = get_observed(specTable,photDict)

    Av = get_Av(observed,intrinsic,args.R)
    for x in Av:
        if x <= 0.0:
            print '--'
        else:
            print x
    exit()
    
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
    photTable = join(photTable,updateTable,join_type='left',keys='oID')
    #print len(photTable)
    #print photTable.colnames
    #exit()

    # Resolve conflicts by choosing just new values
    #  WILL OVERWRITE ANY EXISTING ONES
    for coln in updateTable.colnames:
        if coln+'_2' in photTable.colnames:
            photTable.remove_column(coln+'_1')
            photTable.rename_column(coln+'_2',coln)


    # Write out
    if args.outfile:
        outFITS = args.outfile
    else:
        outFITS = args.photfile
    #outFITS = '../Drout_ZOMG_w_M31B.fits'#'../Drout_list_ZOMG.fits'
    outTSV = os.path.splitext(outFITS)[0] + '.tsv'
    

    #photTable.write(outFITS)
    ##print 'Writing to %s' % outFITS
    #photTable.write(outTSV,format='ascii.tab')
    for x in photTable['Av']:
        if x <= 0.0:
            print '--'
        else:
            print x
    ###print 'writing to %s' % outTSV
    
    


if __name__ == '__main__':
    main()









