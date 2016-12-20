#! /usr/bin/env python
# 05/03/16
import argparse
from collections import OrderedDict
import numpy as np
from mag2flux import get_wave


def ext_correction_indebetouw(Ak,filter = None):
    #indebetouw 2005
    filters = ['R','I','Rc','Ic','J','H','K','3.6','4.5',
               '5.8','8.0','W1','W2','W3','W4',
               'MIPS24','MSXA','MSXB1','MSXB2','MSXC','MSXD','MSXE',
               'IRAS12', 'IRAS25', 'IRAS60', 'IRAS100',
               'AKARIS9W', 'AKARIN60', 'AKARIWIDES']
    waves = [get_wave(x) for x in filters]
    l_Alam_Ak = 0.61 - 2.22*(np.log10(waves)) + 1.21*(np.log10(waves)**2)
    Alam_Ak = 10.0**l_Alam_Ak
    Alam = Alam_Ak * np.float(Ak)
    Alam[np.where(Alam < 0)] = np.min(np.abs(Alam))  # smallest number
    
    dictak = OrderedDict(zip(filters,Alam))
    if filter:
        return dictak[filter]
    else:
        return dictak

def ext_correction_multiple(Av,Rv=3.2):
    # all values from Cardelli 1989
    filters = ['U','B','V','R','Rc','I','Ic','J','H','K','L']
    x = [2.78,2.27,1.82,1.43,1.43,1.11,1.11,0.80,0.63,0.46,0.29]
    ax = [0.9530,0.9982,1.0000,0.8686,0.8686,0.6800,0.6800,0.4008,0.2693,0.1615,0.0800]
    bx = [1.9090,1.0495,0.0,-0.3660,-0.3660,-0.6239,-0.6239,-0.3679,-0.2473,-0.1483,-0.0734]

    ext_law = [a+(b/Rv) for a,b in zip(ax,bx)]
    ext_dicts = []

    for A in list(Av):
        print A
        if (A or A != 0) and (A != '--'):
            al = [float(A)*ext for ext in ext_law]
        else:
            al = [None for ext in ext_law]
            
        dict_av = OrderedDict(zip(filters,al))
        ext_dicts.append(dict_av)

    return ext_dicts

def ext_cardelli_smooth(waves,Av,Rv=3.2):
    # waves in microns
    opt_idx = np.where(waves < 1.1)

    a = np.zeros(len(waves))
    b = np.zeros(len(waves))

    a[opt_idx] = 0.574*np.power(1./waves[opt_idx],1.61)
    b[opt_idx] = -0.527*np.power(1./waves[opt_idx],1.61)

    nir_idx = np.where(waves >= 1.1)
    y = (1./waves[nir_idx]) - 1.82
    a[nir_idx] = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
    b[nir_idx] = 1.41338*y + 2.28305*y**2 + 1.07223*y**3 - 5.38434*y**4 -0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7

    ext = Av*(a+(b/Rv))
    return ext
    
    
    

def ext_correction(Av,Rv=3.2,filter = None):
    # all values from Cardelli 1989
    filters = ['U','B','V','R','Rc','Ic','I','J','H','K','L']
    x = [2.78,2.27,1.82,1.43,1.43,1.11,0.80,0.63,0.46,0.29]
    ax = [0.9530,0.9982,1.0000,0.8686,0.8686,0.6800,0.6800,0.4008,0.2693,0.1615,0.0800]
    bx = [1.9090,1.0495,0.0,-0.3660,-0.3660,-0.6239,-0.6239,-0.3679,-0.2473,-0.1483,-0.0734]

    ext_law = [a+(b/Rv) for a,b in zip(ax,bx)]
    A = [float(Av)*ext for ext in ext_law]

    dict_av = dict(zip(filters,A))
    if filter:
        if filter in dict_av:
            return dict_av[filter]
        else:
            return 0.0
    else:
        return dict(zip(filters,A))




def main():
    parser = argparse.ArgumentParser(description='Given Av, return all A_lambda from Cardelli et al 1989')

    parser.add_argument('Av',type=float,help='Extinction in visual')
    parser.add_argument('--Rv',type=float, default=3.1, help='ISM parameter, default=3.1')

    args = parser.parse_args()

    ext_dict = ext_correction(args.Av,args.Rv)
    print ext_dict
    return 0


if __name__ == '__main__':
    main()
