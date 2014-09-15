#! /usr/bin/env python
import pyfits
import numpy as np
import argparse

_MaxIter = 500
_GridSize = 3

def ijIndFromXYPos(xyPos):
    """Return the integer index of the pixel whose center is nearest the specified position."""
    return [int(round(xyPos[ii] - 0.5)) for ii in (1, 0)]


def centroid(data,xyguess,rad):
    if isinstance(data,str):
        data = pyfits.getdata(data)
        
    if len(xyguess) != 2:
        raise ValueError("Initial guess=%r must have 2 elements" % (xyGuess,))

    # Get index of pixel closest to xyguess
    ijguess = ijIndFromXYPos(xyguess)

    # Use this as first guess at maximum.  Extract radial profiles in
    # a 3x3 gridlet around this point, and walk to find minimum fitting error
    maxi, maxj = ijguess
        
    asymmArr = np.zeros([_GridSize,_GridSize], float)
    totPtsArr = np.zeros([_GridSize,_GridSize], int)
    totCountsArr = np.zeros([_GridSize,_GridSize], float)

    niter = 0
    while True:
        niter += 1
        if niter > _MaxIter:
            raise RuntimeError("Could not find object in %i iterations" % _MaxIter)

        for i in range(_GridSize):
            ii = maxi + i - 1
            for j in range(_GridSize):
                kk = maxj + j - 1
                if totPtsArr[i,j] != 0:
                    continue
                asymmArr[i,j],totCountsArr[i,j],totPtsArr[i,j] = radProf
