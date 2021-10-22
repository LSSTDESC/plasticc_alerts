#!/usr/bin/env python
# coding: utf-8

# Simple alert to run on cori: will ingest a dictionary rather than an SNANA file and write out the alerts based on that object.

import numpy as np
import matplotlib.pyplot as plt


import lsst.alert.packet
from pathlib import Path
import shutil
from fastavro import writer, reader

import os
from copy import copy
import json
from collections import OrderedDict as odict
import tarfile
import glob

from astropy.table import Table, vstack
import pandas as pd
from astropy.io import fits
from astropy.time import Time


# taken from sncosmo
def read_snana_fits(head_file, phot_file, snids=None, n=None):
    """Read the SNANA FITS format: two FITS files jointly representing
    metadata and photometry for a set of SNe.

    Parameters
    ----------
    head_file : str
        Filename of "HEAD" ("header") FITS file.
    phot_file : str
        Filename of "PHOT" ("photometry") FITS file.
    snids : list of str, optional
        If given, only return the single entry with the matching SNIDs.
    n : int
        If given, only return the first `n` entries.

    Returns
    -------
    sne : list of `~astropy.table.Table`
        Each item in the list is an astropy Table instance.

    Notes
    -----
    If `head_file` contains a column 'SNID' containing strings, leading and
    trailing whitespace is stripped from all the values in that column.

    If `phot_file` contains a column 'FLT', leading and trailing whitespace
    is stripped from all the values in that column.

    Examples
    --------
    >>> sne = read_snana_fits('HEAD.fits', 'PHOT.fits')
    >>> for sn in sne:
    ...     sn.meta  # Metadata in an OrderedDict.
    ...     sn['MJD']  # MJD column

    """

    
    memmap = (snids is not None or n is not None)
    # Get metadata for all the SNe
    head_data = fits.getdata(head_file, 1, view=np.ndarray)
    phot_data = fits.getdata(phot_file, 1, view=np.ndarray, memmap=memmap)

    # Strip trailing whitespace characters from SNID.
    if 'SNID' in head_data.dtype.names:
        try:
            head_data['SNID'][:] = np.char.strip(head_data['SNID'])
        except TypeError:
            pass

    # Check which indicies to return.
    if snids is None and n is None:
        idx = range(len(head_data))
    elif n is None:
        if 'SNID' not in head_data.dtype.names:
            raise RuntimeError('Specific snids requested, but head file does'
                               ' not contain SNID column')
        idx = []
        for snid in snids:
            i = np.flatnonzero(head_data['SNID'] == snid)
            if len(i) != 1:
                raise RuntimeError('Unique snid requested, but there are '
                                   '{0:d} matching entries'.format(len(i)))
            idx.append(i[0])
    elif snids is None:
        idx = range(n)
    else:
        raise ValueError("cannot specify both 'snids' and 'n' arguments")

    # Loop over SNe in HEAD file
    sne = []
    for i in idx:
        meta = odict(zip(head_data.dtype.names, head_data[i]))

        j0 = head_data['PTROBS_MIN'][i] - 1
        j1 = head_data['PTROBS_MAX'][i]
        data = phot_data[j0:j1]
        if 'BAND' in data.dtype.names:
            data['BAND'][:] = np.char.strip(data['BAND'])
        sne.append(Table(data, meta=meta, copy=False))

    return sne


def import_diasrc_from_fits(sn, my_diasrc,i=0):
        
        my_diasrc['diaObjectId'] = 5 #sn['SNID'][i]
        my_diasrc['filterName'] = sn['BAND'][i]     
        my_diasrc['apFlux'] = sn['FLUXCAL'][i]
        my_diasrc['apFluxErr'] = sn['FLUXCALERR'][i]
        my_diasrc['snr'] = sn['FLUXCAL'][i]/sn['FLUXCALERR'][i]
        my_diasrc['midPointTai'] = sn['MJD'][i]
    
        # General properties
        my_diasrc['ra'] = sn.meta['RA']
        my_diasrc['decl'] = sn.meta['DEC']
        my_diasrc['nobs'] = sn.meta['NOBS']
        my_diasrc['mwebv'] = sn.meta['MWEBV']
        my_diasrc['mwebv_err'] = sn.meta['MWEBV_ERR']
        my_diasrc['z_final'] = sn.meta['REDSHIFT_FINAL']
        my_diasrc['z_final_err'] = sn.meta['REDSHIFT_FINAL_ERR']
        my_diasrc['hostgal_ra'] = sn.meta['HOSTGAL_RA']
        my_diasrc['hostgal_dec'] = sn.meta['HOSTGAL_DEC']
        my_diasrc['hostgal_snsep']= sn.meta['HOSTGAL_SNSEP']
        my_diasrc['hostgal_z']=sn.meta['HOSTGAL_SPECZ']
        my_diasrc['hostgal_z_err']=sn.meta['HOSTGAL_SPECZ_ERR']   
        my_diasrc['hostgal_mag_u']= sn.meta['HOSTGAL_MAG_u']
        my_diasrc['hostgal_mag_g']= sn.meta['HOSTGAL_MAG_g']
        my_diasrc['hostgal_mag_r']= sn.meta['HOSTGAL_MAG_r']
        my_diasrc['hostgal_mag_i']= sn.meta['HOSTGAL_MAG_i']
        my_diasrc['hostgal_mag_z']= sn.meta['HOSTGAL_MAG_z']
        my_diasrc['hostgal_mag_Y']= sn.meta['HOSTGAL_MAG_Y']
        my_diasrc['hostgal_magerr_u']= sn.meta['HOSTGAL_MAGERR_u']
        my_diasrc['hostgal_magerr_g']= sn.meta['HOSTGAL_MAGERR_g']
        my_diasrc['hostgal_magerr_r']= sn.meta['HOSTGAL_MAGERR_r']
        my_diasrc['hostgal_magerr_i']= sn.meta['HOSTGAL_MAGERR_i']
        my_diasrc['hostgal_magerr_z']= sn.meta['HOSTGAL_MAGERR_z']
        my_diasrc['hostgal_ellipticity'] = 5 #sn.meta['HOSTGAL_ELLIPTICITY']
        my_diasrc['hostgal_sqradius'] = 5 # sn.meta['HOSTGAL_SQRADIUS']
        my_diasrc['hostgal2_ellipticity'] = 5 # sn.meta['HOSTGAL2_ELLIPTICITY']
        my_diasrc['hostgal2_size'] = 5 #sn.meta['HOSTGAL2_SQRADIUS']

        return my_diasrc

## Data-frame specific code

def make_dataframe(packet):
    df = pd.DataFrame(packet['diaSource'], index=[0])
    df_prv = pd.DataFrame(packet['prvDiaSources'])
    return pd.concat([df,df_prv], ignore_index=True)


####################################

rundir='./'
schemadir = rundir+'/plasticc_schema/'
modelname='MLAG_HOSTLIB_ALLCLASS_NONIaMODEL0'


schema = lsst.alert.packet.Schema.from_file(filename=schemadir+'lsst.v4_1.alert.avsc')
# Load an example json alert, and clear the numberical input
print(schemadir+'sample_data/plasticc.json')

with open(schemadir+'sample_data/plasticc.json') as f:
    alert_data = json.load(f)
    
alert_data_orig = alert_data.copy()
#%json.load(f)
    
diasrc = alert_data['prvDiaSources'][0]

num_sne =1 #how many supernova alerts you want to make
num_mjd = 1 # how many mjd you want to plot for this object
print(f"Running alerts for the model list: {modelname}")


name_head=modelname+'/'+'%s-0010_HEAD.FITS.gz'%modelname
name_phot=modelname+'/'+'%s-0010_PHOT.FITS.gz'%modelname
sne = read_snana_fits(name_head, name_phot)


for iterNum,sn in enumerate(sne[0:num_sne]):
    ra = sn.meta['RA']
    mjd = sn['MJD']
    
    alert = copy(alert_data_orig)  
    diasrc = alert_data_orig['prvDiaSources'][0]
    my_diasrc = copy(diasrc)
    alert = copy(alert_data_orig)
    alert['prvDiaSources'].clear()
    
    # Make a dummy source ID and import additional information to the diaSource
    my_diasrc['diaSourceId'] = np.int64(28132306237521+1000*np.random.uniform())
    my_diasrc = import_diasrc_from_fits(sn,my_diasrc,i=0)
    alert['diaSource'] = my_diasrc
    
    
    for j, day in enumerate(mjd[0:num_mjd]):
        print(mjd[j])
        my_diasrc = copy(diasrc)
        my_diasrc['diaSourceId'] = my_diasrc['diaSourceId']+j+1
        my_diasrc = import_diasrc_from_fits(sn,my_diasrc,i=j)
        alert['prvDiaSources'].append(alert['diaSource'])
        
        # serialize the alert    
        avro_bytes = schema.serialize(alert)
        messg = schema.deserialize(avro_bytes)
        
        mjdlabel=mjd[j]
        with open("%.1f_%i_%i.avro"%(mjdlabel,my_diasrc['diaObjectId'],my_diasrc['diaSourceId']), "wb") as f:
            schema.store_alerts(f, [alert])
            
    print(f"Done saving models for {modelname}")

