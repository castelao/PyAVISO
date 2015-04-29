#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 expandtab


"""
"""

import os
import os.path
from datetime import datetime, timedelta
import time
import re
import logging
import logging.handlers

import numpy
import numpy as np
from numpy import ma

try:
    import netCDF4
    from netCDF4 import date2num
except:
    import pupynere
from pydap.client import open_url

try:
    from rings import okuboweiss
except:
    pass


def products(ncfile):
    """ Calculate some products from the downloaded data

        Think about, should I filter the data here or during
          the analysis? Probably over there, but if I do here
          I'll need to include some metadata, like filter type
          and half power cut freq.
    """
    assert os.path.isfile(ncfile)

    nc = netCDF4.Dataset(ncfile, 'a', format='NETCDF4')

    W = nc.createVariable('W', 'f4', ('time', 'latitude', 'longitude'),
            fill_value=netCDF4.default_fillvals['f4'])
    zeta = nc.createVariable('zeta', 'f4', ('time', 'latitude', 'longitude'),
            fill_value=netCDF4.default_fillvals['f4'])

    for tn in range(nc.variables['time'].size):
        data = {'Lat': nc.variables['Lat'][:],
            'Lon': nc.variables['Lon'][:],
            'u': nc.variables['u'][tn],
            'v': nc.variables['v'][tn],
            }

        products = okuboweiss.okuboweiss(data)

        W[tn] = products['W']
        zeta[tn] = products['zeta']

    nc.close()


def mask_shallow(ncfile, zlimit=-150, zfile=None):
    """ Mask all variables in vars @ gridpoints shallower than mindepth

         Use http://opendap.ccst.inpe.br/Misc/etopo2/ETOPO2v2c_f4.nc

         In the future move it out of here, into a different support
           command line
    """
    if zfile is None:
        zfile = "http://opendap.ccst.inpe.br/Misc/etopo2/ETOPO2v2c_f4.nc"

    nc = netCDF4.Dataset(ncfile, 'a')
    Lat = nc.variables['latitude'][:]
    Lon = nc.variables['longitude'][:]

    #Lon,Lat = numpy.meshgrid(self.data['lon'],self.data['lat'])
    #from fluid.common.common import get_bathymery
    #z = get_bathymery(Lat, Lon, etopo_file=zfile)
    # ========================================================================
    # Not cute, but works for now.
    ncz = netCDF4.Dataset(zfile, 'r')
    # -180:180
    lon_z = ncz.variables['x'][:]
    lat_z = ncz.variables['y'][:]

    # If input 0:360
    ind = Lon > 180
    Lon[ind] = Lon[ind] - 360

    I = Lat.size
    J = Lon.size
    z = np.empty((I, J))
    for i in range(I):
        for j in range(J):
            x_ind = (np.absolute(Lon[j] - lon_z)).argmin()
            y_ind = (np.absolute(Lat[i] - lat_z)).argmin()
            z[i, j] = ncz.variables['z'][y_ind, x_ind]

    # ========================================================================
    ind = z > zlimit
    for v in nc.variables.keys():
        if nc.variables[v].dimensions == (u'time', u'latitude', u'longitude'):
            if nc.variables[v].shape[1:] != ind.shape:
                return
            I, J = np.nonzero(ind)
            for i, j in zip(I, J):
                nc.variables[v][:,i,j] = nc.variables[v]._FillValue
    nc.sync()
    nc.close()

def mask(self):
    """ Improve it. Make it more flexible
    """
    #
    print "Masking data shallower then: %s" % self.metadata['mask_shallow']
    bath_mask = self.data['z']>self.metadata['mask_shallow']
    #
    Lon,Lat = numpy.meshgrid(self.data['lon'],self.data['lat'])
    equator_mask = (Lat>-2.5) & (Lat<2.5)
    #
    for i in range(len(self.data['datetime'])):
        self.data['u'][i,:,:]=ma.masked_array(self.data['u'][i,:,:].data,mask=(self.data['u'][i,:,:].mask) | bath_mask | equator_mask)
        self.data['v'][i,:,:]=ma.masked_array(self.data['v'][i,:,:].data,mask=(self.data['v'][i,:,:].mask) | bath_mask | equator_mask)
        self.data['h'][i,:,:]=ma.masked_array(self.data['h'][i,:,:].data,mask=(self.data['h'][i,:,:].mask) | bath_mask)


def eke(cutperiod=360, dt=7, verbose=False):
    """
        Include the possibility to do with a different dataset, like anomaly or ref.

        ATENTION, need to move user and password out of here.
    """
    from maud import window_1Dmean
    l = cutperiod*24    # From days to hours. Aviso time is on hours.

    #self.metadata['urlbase'] = "http://%s:%s@opendap.aviso.oceanobs.com/thredds/dodsC" % (self.metadata['username'], self.metadata['password'])

    url_uv = "http://aviso-users:grid2010@opendap.aviso.oceanobs.com/thredds/dodsC/dataset-duacs-dt-upd-global-merged-madt-uv-daily"
    dataset = open_url(url_uv)

    T, I, J = dataset.Grid_0001.shape
    eke = ma.masked_all((I,J))

    I,J = numpy.nonzero(ma.masked_values(dataset.Grid_0001.Grid_0001[-300::60,:,:], dataset.Grid_0001.attributes['_FillValue']).max(axis=0))
    t = ma.array(dataset.time[::dt])
    if verbose:
        from progressbar import ProgressBar
        pbar = ProgressBar(maxval=I.shape[0]).start()
        n=-1
    for i, j in zip(I,J):
        if verbose:
            n+=1
            pbar.update(n)
        doit = True
        while doit:
            try:
                u = ma.masked_values(dataset.Grid_0001.Grid_0001[::dt,i,j], dataset.Grid_0001.attributes['_FillValue'])*1e-2
                v = ma.masked_values(dataset.Grid_0002.Grid_0002[::dt,i,j], dataset.Grid_0002.attributes['_FillValue'])*1e-2
                u_prime = u-window_1Dmean(u, l=l, t=t, axis=0)
                v_prime = v-window_1Dmean(v, l=l, t=t, axis=0)
                eke[i,j] = (u_prime**2+v_prime**2).mean()/2.
                doit=False
            except:
                print "I had some trouble. I'll wait a litte bit and try again"
                time.sleep(10)
    if verbose:
        pbar.finish()

    return eke
