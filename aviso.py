#!/usr/bin/env python
# -*- coding: Latin-1 -*-
# vim: tabstop=4 shiftwidth=4 expandtab


""" Module to deal with AVISO datasets
"""

import os
import os.path
from datetime import datetime, timedelta
import time
import re
import logging
import logging.handlers

import numpy
from numpy import ma

import pupynere
from pydap.client import open_url

# Time to refactor and change somethings. Once class to download, and another
#   on the top of that to handle the file like it was a MA

class AVISO_fetch(object):
    """ Class to fetch maps from AVISO

        - Deal with the file to save the data
        - Download
            - Download LatLon
                - Adjust -180/180 <-> 0/360
                - First to second limit
            - Download time
                - Parse data input into time index
            - Define the size of the blocks (number of time snapshots per download)
            - Download data in blocks
            - Think about how to save the data in a generic way? This class should worry only about download and save it. Another Class should be created on the top of this to use it and offer the downloaded data as a MA.
    """

    def __init__(self, cfg):
        """
        """
        self.cfg = cfg

        self.set_logger()
        self.logger.info("Initializing AVISO_fetch class")
        self.logger.debug("cfg: %s" % cfg)

        if ('username' not in self.cfg) | ('password' not in self.cfg):
            self.logger.error("Aviso DAP server requires a registered username and password. I'll abort.")
            return

        if 'urlbase' not in self.cfg:
            self.cfg['urlbase'] = \
                    "http://%s:%s@opendap.aviso.oceanobs.com/thredds/dodsC" % \
                    (self.cfg['username'], self.cfg['password'])
        self.logger.debug("urlbase: %s" % self.cfg['urlbase'])

        if 'force_download' not in self.cfg:
            self.cfg['force_download'] = False

        self.set_source_filename()
        self.set_dataset()

        self.file = os.path.join(self.cfg['datadir'], self.cfg['filename'])
        self.nc = pupynere.netcdf_file(self.file,'w')

        self.download_time()
        self.download_LonLat()
        self.download_data()


    def set_logger(self):
        """
        """
        # Creating another log level
        logging.VERBOSE = logging.DEBUG - 1
        logging.addLevelName(logging.VERBOSE, 'VERBOSE')

        #create logger
        logger = logging.getLogger("AVISO fetch")
        logger.setLevel(logging.VERBOSE)

        #create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(logging.WARN)
        #create formatter
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s -  %(message)s")
        #add formatter to ch
        ch.setFormatter(formatter)
        #add ch to logger
        logger.addHandler(ch)

        if 'logfile' in  self.cfg:
            #create a rotate file handler
            fh = logging.handlers.RotatingFileHandler(
                    self.cfg['logfile'],
                    mode='a', maxBytes=1000000, backupCount=10)
            fh.setLevel(logging.DEBUG)
            #create formatter
            formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s -  %(message)s")
            #add formatter to ch
            fh.setFormatter(formatter)
            #add ch to logger
            logger.addHandler(fh)

        #"application" code
        #logger.debug("debug message")
        #logger.info("info message")
        #logger.warn("warn message")
        #logger.error("error message")
        #logger.critical("critical message")

        self.logger = logger


    def set_source_filename(self):
        """
        """
        self.cfg['source_filename'] = "dataset-duacs-dt-%s-global-merged-%s" % (self.cfg['type'], self.cfg['map'])
        self.logger.debug("source_filename: %s" % self.cfg['source_filename'])

    def set_dataset(self):
        """
        """
        self.logger.debug("Setting the dataset on the DAP server")

        self.url = {'h': "%s/%s-h-daily" % (self.cfg['urlbase'], self.cfg['source_filename']), 
                'uv': "%s/%s-uv-daily" % (self.cfg['urlbase'], self.cfg['source_filename'])}
        self.dataset = {'h': open_url(self.url['h']), 'uv': open_url(self.url['uv'])}

        self.logger.debug("Dataset on the DAP server is ready")

    def download_time(self):
        """
        """
        self.logger.debug("Downloading time")
        if 't_ini' not in self.cfg['limits']:
            self.cfg['limits']['t_ini'] = 0
            self.logger.debug("Setting t_ini: %s" % self.cfg['limits']['t_ini'])

        if 't_step' not in self.cfg['limits']:
            self.cfg['limits']['t_step'] = 1
            self.logger.debug("Setting t_step: %s" % self.cfg['limits']['t_step'])

        if 't_fin' not in self.cfg['limits']:
            self.cfg['limits']['t_fin'] = self.dataset['h']['time'].shape[0]
            self.logger.debug("Setting t_fin: %s" % self.cfg['limits']['t_fin'])

        t_ini = self.cfg['limits']['t_ini']
        t_fin = self.cfg['limits']['t_fin']
        t_step = self.cfg['limits']['t_step']
        # ----
        data={}
        #
        #from coards import from_udunits
        #t0=datetime(1950,1,1)
        #if (re.match('^hours since \d{4}-\d{2}-\d{2}$',dataset_h['time'].attributes['units'])):
        #if (re.match('^hours since 1950-01-01',self.dataset['h']['time'].attributes['units'])):
        #    t = self.dataset['h']['time'][t_ini:t_fin:t_step].tolist()
        #    data['datetime'] = numpy.array([t0+timedelta(hours=h) for h in t])
        #else:
        #    self.logger.error("Problems interpreting the time")

        t = self.dataset['h']['time'][t_ini:t_fin:t_step].tolist()

        self.nc.createDimension('time', len(range(t_ini,t_fin,t_step)))
        nct = self.nc.createVariable('time', 'f', ('time', ))
        nct[:] = t
        nct.units = self.dataset['h']['time'].attributes['units']

    def download_LonLat(self):
        """ Download the Lon x Lat coordinates
        """
        self.logger.debug("Downloading LonLat")
        data = {}
        limits = self.cfg['limits']
        Lat = self.dataset['h']['NbLatitudes']
        Lon = self.dataset['h']['NbLongitudes']

        Latlimits = numpy.arange(Lat.shape[0])[(Lat[:]>=limits["latini"]) & (Lat[:]<=limits["latfin"])]
        Latlimits = [Latlimits[0],Latlimits[-1]]

        Lonlimits = numpy.arange(Lon.shape[0])[(Lon[:]>=limits["lonini"]) & (Lon[:]<=limits["lonfin"])]
        Lonlimits=[Lonlimits[0],Lonlimits[-1]]

        self.cfg['limits']['Latlimits'] = Latlimits
        self.cfg['limits']['Lonlimits'] = Lonlimits

        Lon, Lat = numpy.meshgrid( (Lon[Lonlimits[0]:Lonlimits[-1]]), (Lat[Latlimits[0]:Latlimits[-1]]) )


        self.slice_size = (Lonlimits[-1]-Lonlimits[0])* \
                    (Latlimits[-1]-Latlimits[0])

        # ========
        self.nc.createDimension('lat', (Latlimits[-1]-Latlimits[0]))
        self.nc.createDimension('lon', (Lonlimits[-1]-Lonlimits[0]))

        ncLat = self.nc.createVariable('Lat', 'f', ('lat', 'lon'))
        ncLon = self.nc.createVariable('Lon', 'f', ('lat', 'lon'))

        ncLat[:] = Lat
        ncLon[:] = Lon

    def download_data(self):
        """ Download h and uv in blocks
        """

        # Will download blocks of at most 5MB
        #   i.e. 4e7 floats of 32bits.
        dblocks = max(1, int(4e7/self.slice_size))

        ti = numpy.arange(self.cfg['limits']['t_ini'], 
                self.cfg['limits']['t_fin'], 
                self.cfg['limits']['t_step'])

        blocks = ti[::dblocks]
        if ti[-1] not in blocks:
            blocks = numpy.append(blocks, self.cfg['limits']['t_fin'])

        ntries = 40
        #------
        data = {}
        for v, dataset, missing_value in zip(['h','u','v'], 
                [self.dataset['h']['Grid_0001']['Grid_0001'], 
                    self.dataset['uv']['Grid_0001']['Grid_0001'], 
                    self.dataset['uv']['Grid_0002']['Grid_0002']], 
                [self.dataset['h']['Grid_0001']._FillValue, 
                    self.dataset['uv']['Grid_0001']._FillValue, 
                    self.dataset['uv']['Grid_0002']._FillValue]):

            print "Getting %s" % v
            #data['h'] = ma.masked_all((len(ti),Lonlimits[-1]-Lonlimits[0], Latlimits[-1]-Latlimits[0]), dtype=numpy.float64)
            data[v] = self.nc.createVariable(v, 'f', ('time', 'lat', 'lon'))
            data[v].missing_value = missing_value

            # Work on these limits. Must have a better way to handle it
            Lonlimits = self.cfg['limits']['Lonlimits']
            Latlimits = self.cfg['limits']['Latlimits']

            for b1, b2 in zip(blocks[:-1], blocks[1:]):
                print "From %s to %s of %s" % (b1, b2, blocks[-1])
                ind = numpy.nonzero((ti>=b1) & (ti<b2))
                for i in range(ntries):
                    print "Try n: %s" % i
                    try:
                        data[v][ind] = dataset[b1:b2:self.cfg['limits']['t_step'], Lonlimits[0]:Lonlimits[-1],Latlimits[0]:Latlimits[-1]].swapaxes(1,2).astype('f')
                        break
                    except:
                        waitingtime = 30+i*20
                        print "Failed to download. I'll try again in %ss" % waitingtime
                        time.sleep(waitingtime)
            #data['h'] = 1e-2*data['h'].swapaxes(1,2)





class Aviso_map(object):
    """ Class to get the maps of h and uv from Aviso

        ATENTION, should improve it. Should include a check on the file if
          is the same time and area coverage.
    """
    def __init__(self, metadata=None, auto=True):
        """
        """
        self.metadata = metadata or {}

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
        if 'datadir' not in self.metadata:
            self.metadata['datadir'] = "./" #"../data"
        if 'urlbase' not in self.metadata:
            # Double check this. I believe now it is required to have a password,
            #   so I should clean the option with no use/pass
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
        return

    def download(self):
        """

            Migrate it to use np.lib.arrayterator.Arrayterator
        """
        url_h = "%s/%s-h-daily" % (self.metadata['urlbase'], self.metadata['source_filename'])
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
            blocks = numpy.append(blocks,t_fin)

        ntries = 40
        #------
        for v, dataset, missing_value in zip(['h','u','v'], [dataset_h['Grid_0001']['Grid_0001'], dataset_uv['Grid_0001']['Grid_0001'], dataset_uv['Grid_0002']['Grid_0002']], [dataset_h['Grid_0001']._FillValue, dataset_uv['Grid_0001']._FillValue, dataset_uv['Grid_0002']._FillValue]):

            print "Getting %s" % v
            #data['h'] = ma.masked_all((len(ti),Lonlimits[-1]-Lonlimits[0], Latlimits[-1]-Latlimits[0]), dtype=numpy.float64)
            self.data[v] = nc.createVariable(v, 'f', ('time', 'lat', 'lon'))
            self.data[v].missing_value = missing_value
            for b1, b2 in zip(blocks[:-1], blocks[1:]):
                print "From %s to %s of %s" % (b1, b2, blocks[-1])
                ind = numpy.nonzero((ti>=b1) & (ti<b2))
                for i in range(ntries):
                    print "Try n: %s" % i
                    try:
                        self.data[v][ind] = dataset[b1:b2:t_step, Lonlimits[0]:Lonlimits[-1],Latlimits[0]:Latlimits[-1]].swapaxes(1,2).astype('f')
                        break
                    except:
                        waitingtime = 30+i*20
                        print "Failed to download. I'll try again in %ss" % waitingtime
                        time.sleep(waitingtime)
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
