#! /usr/bin/env python
from astropy.table import Table,MaskedColumn,Column,hstack,vstack
from astroquery.vizier import Vizier
from astroquery.irsa import Irsa
from functools import partial
import astropy.units as u
import numpy as np

def query_cat(coords,catalog,radius=0.5*u.arcsec,cols=None,fill_val=-99.99,full=False,IRSA=False):
    #one-to-one query
    if IRSA:
        results = Irsa.query_region(coords,catalog=catalog,radius=radius)
    else:
        results = Vizier.query_region(coords,catalog=catalog,radius=radius)
    if len(results) == 0:
        return None

    if full:
        return results

    results = results[0]

    if cols is not None:
        # if dict, remap colnames
        if isinstance(cols,dict):
            for k,v in cols.iteritems():
                results.rename_column(k,v)
            names = cols.values()
        else:
            names = cols
    else:
        names = results.colnames

    # make new columns one-to-one with coords
    newtable = Table(masked=True)
    for col in names:
        oldcol = results[col]
        newcol = MaskedColumn(data=np.zeros(len(coords),dtype=oldcol.dtype),unit=oldcol.unit,name=col,mask=np.ones(len(coords),dtype=bool),fill_value=fill_val)

        # copy data from results
        for row in results:
            if not row[col]:
                continue
            # _q IS 1-BASED INDEXING?!
            newcol[row['_q']-1] = row[col]
            newcol.mask[row['_q']-1] = False

        newtable.add_column(newcol)

    return newtable


query_irsa = partial(query_cat,IRSA=True)
query_vizier = partial(query_cat,IRSA=False)

