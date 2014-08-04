#! /usr/bin/env python
import argparse
from astropy.table import Table
import sys
import os
import numpy as np

def get_Av(UB, BV):
    if len(UB) != len(BV):
        sys.exit('UB and BV lists must be same length')

    Q = np.array(UB) - 0.72 * np.array(BV)

    BV0 = 0.332*Q
    EBV = BV - BV0
    Av = 3.2 * EBV

    return Av

        
        

def main():
    parser = argparse.ArgumentParser(description='Calculate Av via Qmethod')
    parser.add_argument('-UB',required=True,type=float,nargs='+',help='U-B color')
    parser.add_argument('-BV',required=True,type=float,nargs='+',help='B-V color')

    args = parser.parse_args()
    if len(args.UB) != len(args.BV):
        sys.exit('UB and BV lists must be same length')
    

    print get_Av(args.UB,args.BV)


if __name__ == '__main__':
    main()
