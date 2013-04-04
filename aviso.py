#!/usr/bin/env python
# -*- coding: Latin-1 -*-
# vim: tabstop=4 shiftwidth=4 expandtab


""" Module to deal with AVISO datasets

    Migrate from pickle to shove
    (https://bitbucket.org/luizirber/extract/src/0403757ce04d/tests/extract_simple.py#cl-37)
    I'm not sure if that's a good idea. I had several problems with the shove in other project.
"""

import os
import os.path
from datetime import datetime, timedelta
import time
from UserDict import UserDict
import re

import numpy
from numpy import ma

import pupynere
from pydap.client import open_url
from maud import window_1Dmean


class Aviso_map(object):
    """ Class to get the maps of h and uv from Aviso

        ATENTION, should improve it. Should include a check on the file if
          is the same time and area coverage.
    """
    def __init__(self, metadata=None, auto=True):
        """
        """
        self.metadata = metadata or {}
        #
        self.set_default()
        if self.metadata['map'] == 'madt+msla':
            print "Map: madt+msla"
            print "First getting madt"
            m = self.metadata.copy()
            m['map'] = 'madt'
            data_madt = Aviso_map(metadata=m)
            print "Now getting msla"
            m = self.metadata.copy()
            m['map'] = 'msla'
            data_msla = Aviso_map(metadata=m)
            if (data_madt['Lat']!=data_msla['Lat']).any() | \
                    (data_madt['Lon']!=data_msla['Lon']).any() | \
                    (data_madt['datetime']!=data_msla['datetime']).any():
                print("PROBLEMS!!!!! data_madt and data_msla are not equivalent!")
                return
            else:
                self.data = data_madt.data
                self.data['ssh'] = self.data['h']
                del(self.data['h'])
                self.data['eta'] = data_msla.data['h']
                self.data['u_anom'] = data_msla.data['u']
                self.data['v_anom'] = data_msla.data['v']
                #del(data_madt)
                #del(data_msla)
                #self.data = data
        elif (self.metadata['map'] == 'madt') or (self.metadata['map'] == 'msla'):
            self.set_source_filename()
            if auto==True:
                self.get_data()

    def __getitem__(self, key):
        """
        """
        if type(self.data[key]) == 'numpy.ma.core.MaskedArray':
            return self.data[key]
        elif hasattr(self.data[key], 'missing_value'):
            return ma.masked_values(self.data[key][:], getattr(self.data[key], 'missing_value'))
        return ma.array(self.data[key])

    def __setitem__(self, key, value):
        """
        """
        self.data[key] = value

    def keys(self):
        return self.data.keys()

    def set_default(self):
        """
        """
        if ('username' not in self.metadata) | ('password' not in self.metadata):
            print("Aviso DAP server requires a registered username and password, sorry.")
            return
        if 'type' not in self.metadata:
            self.metadata['type'] = "upd"       # upd, ref
        if 'map' not in self.metadata:
            self.metadata['map'] = "madt"       # madt, msla
        if 'limits' not in self.metadata:
            self.metadata['limits'] = {'latini': 0, 'latfin': 15, 'lonini': 296, 'lonfin': 317}
        if 'mask_shallow' not in self.metadata:
            self.metadata['mask_shallow'] = -250
        if 'datadir' not in self.metadata:
            self.metadata['datadir'] = "../data"
        if 'urlbase' not in self.metadata:
            if ('username' in self.metadata) & ('password' in self.metadata):
                self.metadata['urlbase'] = "http://%s:%s@opendap.aviso.oceanobs.com/thredds/dodsC" % (self.metadata['username'], self.metadata['password'])
            else:
                self.metadata['urlbase'] = "http://opendap.aviso.oceanobs.com/thredds/dodsC"
        if 'force_download' not in self.metadata:
            self.metadata['force_download'] = False
        return

    def set_source_filename(self):
        """
        """
        self.metadata['source_filename'] = "dataset-duacs-dt-%s-global-merged-%s" % (self.metadata['type'],self.metadata['map'])
        return

    def get_data(self):
        """
        """
        #file = os.path.join(self.metadata['datadir'],self.metadata['source_filename']+".pkl.bz2")
        self.file = os.path.join(self.metadata['datadir'],self.metadata['source_filename']+".nc")
        #self.nc = pupynere.netcdf_file(file, 'w')
        #self.download()


        if self.metadata['force_download'] == True:
            self.download()
        else:
            if os.path.isfile(file):
                print "I already downloaded that. I'm just reloading"
            else:
                print "I need to download it. This might take a while. (%s)" % datetime.now()
                self.download()
                if os.path.isdir(self.metadata['datadir']) == False:
                    print(" There is no data directory: %s" & self.metadata['datadir'])
                try:
                    pass
                except:
                    print "Couldn't save the data on pickle"
        self.set_z()
        #self.mask()
        return

    def download(self):
        """

            Migrate it to use np.lib.arrayterator.Arrayterator
        """
        url_h = "%s/%s-h-daily" % (self.metadata['urlbase'], self.metadata['source_filename'])
        print url_h
        dataset_h = open_url(url_h)
        url_uv = "%s/%s-uv-daily" % (self.metadata['urlbase'], self.metadata['source_filename'])
        dataset_uv = open_url(url_uv)
        # ----
        if 't_ini' not in self.metadata['limits']:
            self.metadata['limits']['t_ini'] = 0
        if 't_fin' not in self.metadata['limits']:
            self.metadata['limits']['t_fin'] = dataset_h['time'].shape[0]
        if 't_step' not in self.metadata['limits']:
            self.metadata['limits']['t_step'] = 0
        else:
            print "Atention!! t_step set to: %s" % self.metadata['limits']['t_step']
        t_ini = self.metadata['limits']['t_ini']
        t_fin = self.metadata['limits']['t_fin']
        t_step = self.metadata['limits']['t_step']
        # ----
        data={}
        #
        #from coards import from_udunits
        t0=datetime(1950,1,1)
        #if (re.match('^hours since \d{4}-\d{2}-\d{2}$',dataset_h['time'].attributes['units'])):
        if (re.match('^hours since 1950-01-01',dataset_h['time'].attributes['units'])):
            data['datetime']=numpy.array([t0+timedelta(hours=h) for h in dataset_h['time'][t_ini:t_fin:t_step].tolist()])
        else:
            print "Problems interpreting the time"
            return

        #self.nc.createDimension('time', len(range(t_ini,t_fin,t_step)))
        #time = self.nc.createVariable('time', 'i', ('time',))
        #time[:] = dataset_h['time'][t_ini:t_fin:t_step]
        #time.units = dataset_h['time'].attributes['units']
        #data['time'] = time
        #
        limits=self.metadata['limits']
        Lat=dataset_h['NbLatitudes']
        Lon=dataset_h['NbLongitudes']
        Latlimits=numpy.arange(Lat.shape[0])[(Lat[:]>=limits["latini"]) & (Lat[:]<=limits["latfin"])]
        Latlimits=[Latlimits[0],Latlimits[-1]]
        Lonlimits=numpy.arange(Lon.shape[0])[(Lon[:]>=limits["lonini"]) & (Lon[:]<=limits["lonfin"])]
        Lonlimits=[Lonlimits[0],Lonlimits[-1]]

        #data['Lat'], data['Lon'] = numpy.meshgrid( (Lat[Latlimits[0]:Latlimits[-1]]), (Lon[Lonlimits[0]:Lonlimits[-1]]))
        data['Lon'], data['Lat'] = numpy.meshgrid( (Lon[Lonlimits[0]:Lonlimits[-1]]), (Lat[Latlimits[0]:Latlimits[-1]]) )


        #------
        self.data = data
        #Arrayterator = numpy.lib.arrayterator.Arrayterator
        #dataset = dataset_h['Grid_0001']['Grid_0001']
        #ssh = Arrayterator(dataset)[t_ini:t_fin:t_step]

        #blocks = 1e4


        file = os.path.join(self.metadata['datadir'],self.metadata['source_filename']+".nc")
        nc = pupynere.netcdf_file(file,'w')
        nc.createDimension('time', len(range(t_ini,t_fin,t_step)))
        nc.createDimension('lon', (Lonlimits[-1]-Lonlimits[0]))
        nc.createDimension('lat', (Latlimits[-1]-Latlimits[0]))

        dblocks = max(1,int(1e5/((Lonlimits[-1]-Lonlimits[0])*(Latlimits[-1]-Latlimits[0]))))

        ti = numpy.arange(t_ini, t_fin, t_step)
        blocks = ti[::dblocks]
        if ti[-1] not in blocks:
            #blocks.append(t_fin)
            blocks = numpy.append(blocks,t_fin)

        ntries = 40
        #------
        for v, dataset, missing_value in zip(['h','u','v'], [dataset_h['Grid_0001']['Grid_0001'], dataset_uv['Grid_0001']['Grid_0001'], dataset_uv['Grid_0002']['Grid_0002']], [dataset_h['Grid_0001']._FillValue, dataset_uv['Grid_0001']._FillValue, dataset_uv['Grid_0002']._FillValue]):

            print "Getting %s" % v
            #data['h'] = ma.masked_all((len(ti),Lonlimits[-1]-Lonlimits[0], Latlimits[-1]-Latlimits[0]), dtype=numpy.float64)
            #self.data[v] = nc.createVariable(v, 'f', ('time', 'lon', 'lat'))
            self.data[v] = nc.createVariable(v, 'f', ('time', 'lat', 'lon'))
            self.data[v].missing_value = missing_value
            for b1, b2 in zip(blocks[:-1], blocks[1:]):
                print "From %s to %s of %s" % (b1, b2, blocks[-1])
                ind = numpy.nonzero((ti>=b1) & (ti<b2))
                for i in range(ntries):
                    print "Try n: %s" % i
                    try:
                        #self.data[v][ind] = dataset[b1:b2:t_step, Lonlimits[0]:Lonlimits[-1],Latlimits[0]:Latlimits[-1]]
                        self.data[v][ind] = dataset[b1:b2:t_step, Lonlimits[0]:Lonlimits[-1],Latlimits[0]:Latlimits[-1]].swapaxes(1,2).astype('f')
                        break
                    except:
                        waitingtime = 30+i*20
                        print "Failed to download. I'll try again in %ss" % waitingtime
                        time.sleep(waitingtime)
                        #ssh = dataset_h['Grid_0001']['Grid_0001']
            #data['h'] = 1e-2*data['h'].swapaxes(1,2)



    def set_z(self):
        """

             Use http://opendap.ccst.inpe.br/Misc/etopo2/ETOPO2v2c_f4.nc
        """
        #Lon,Lat = numpy.meshgrid(self.data['lon'],self.data['lat'])
        from fluid.common.common import get_bathymery
        self.data['z'] = get_bathymery(self.data['Lat'], self.data['Lon'], etopo_file="/Users/castelao/work/misc/ETOPO2v2c_f4.nc")

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


#class MAFromDAP(object):
#    """
#    """
#    def __init__(self, dataset, var):
#        self.dataset = dataset
#        self.var = var
#    def find_missingvalue(self):
#        """ Extract the missing value from the dataset
#        """
#        self.missing_value = getattr(self.dataset, '_FillValue')
#    def __getitem__(self, key):
#        if type(key) == slice:
#
#x = MAFromDAP(dataset_h['Grid_0001'], 'Grid_0001')
#print dir(x)
#print x[::10]
