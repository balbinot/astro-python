#! /usr/bin/env python
from astropy.table import Table
import numpy as np
from scipy.interpolate import interp1d

def get(SpTy, temp = None):
    table = Table.read('/home/gordon/Software/python/colorExcess/supergiants.tsv',format='ascii.tab')

    if SpTy == 'FeII':
        return np.nan
    
    if temp:
        # interpolate
        f = interp1d(table['Teff'],table['BC'],kind='slinear',bounds_error=False)
        return f(temp)
    
    MK = SpTy[0:2]

    if MK in table['MK']:
        idx = np.where(table['MK'] == MK)
        return np.float(table['BC'][idx])

    for idx,x in enumerate(table['MK']):
        if not x:
            continue
        if len(x) < 2:
            table['MK'][idx] = x+'0'
    

    tM,tK = zip(*[(x[0],x[1]) if x else (None,None) for x in table['MK']])

    if MK[0] not in tM:
        return np.nan

    tM = np.array(tM)
    tK = np.array(tK)
    
    tMidx = np.where(tM == MK[0])  # M's are equal
    tKA = np.array([float(x) for x in tK[tMidx]])

    match = np.searchsorted(tKA,float(MK[1]))
    tMidx = tMidx[0]
    if match == 0:
        # return first
        return table['BC'][tMidx[0]]

    elif match == len(tKA)-1:
        # return last
        return table['BC'][tMidx[-1]]

    else:
        prev = table['BC'][tMidx[match-1]]
        next = table['BC'][tMidx[match]]

        return np.mean([prev,next])

    
