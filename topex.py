#!/usr/bin/env python
# -*- coding: Latin-1 -*-
# vim: tabstop=4 shiftwidth=4 expandtab

"""
"""

__author__ = ["Guilherme Castelao <guilherme@castelao.net>", "Roberto De Almeida <rob@pydap.org>"]

import os
import itertools
from datetime import datetime, timedelta
from UserDict import IterableUserDict
try:
    import cPickle as pickle
except ImportError:
    import pickle

from numpy import ma
import dap.client

def topex_time_table(dt_days,dt_seconds,dt_microseconds,base_date=None):
    """
    """
    if base_date is None:
        base_date=datetime(year=1950,month=01,day=01,hour=0,minute=0,second=0)
    t=[]
    for d, s, ms in itertools.izip(dt_days.compressed(),dt_seconds.compressed(),dt_microseconds.compressed()):
	    dt=timedelta(days=int(d),seconds=int(s),microseconds=int(ms))
            t.append(base_date+dt)
    t=ma.masked_equal(t,-1)
    return t

def topex_track_table(ndata,tracks,cycles):
    """
    """
    track_list=[]
    cycle_list=[]
    for track, n, cycle in itertools.izip(tracks.compressed(),ndata.compressed(),cycles.compressed()):
        for i in range(n):
            track_list.append(track)
            cycle_list.append(cycle)
    track_list=ma.masked_equal(track_list,-1)
    cycle_list=ma.masked_equal(cycle_list,-1)
    return cycle_list,track_list

def read_file(filename,vars=['CorSSH','MSS','Bathy']):
    """
    """
    import dap.client
    try:
        dataset = dap.client.open(filename)
    except:
        return
    cycles=ma.masked_equal(dataset['Cycles'][:,0],dataset['Cycles']._FillValue)
    cycles.set_fill_value(dataset['Cycles']._FillValue)
    tracks=ma.masked_equal(dataset['Tracks'][:],dataset['Tracks']._FillValue)
    tracks.set_fill_value(dataset['Tracks']._FillValue)
    ndata=ma.masked_equal(dataset['NbPoints'][:],dataset['NbPoints']._FillValue)
    ndata.set_fill_value(dataset['NbPoints']._FillValue)
    [cycle_list,track_list]=topex_track_table(ndata,tracks,cycles)
    # Time related
    TimeDay=ma.masked_equal(dataset['TimeDay'][:,0],dataset['TimeDay']._FillValue)
    TimeDay.set_fill_value(dataset['TimeDay']._FillValue)
    TimeSec=ma.masked_equal(dataset['TimeSec'][:,0],dataset['TimeSec']._FillValue)
    TimeSec.set_fill_value(dataset['TimeSec']._FillValue)
    TimeMicroSec=ma.masked_equal(dataset['TimeMicroSec'][:,0],dataset['TimeMicroSec']._FillValue)
    TimeMicroSec.set_fill_value(dataset['TimeMicroSec']._FillValue)
    # Improve and include the check of the BeginDates
    time_list=topex_time_table(TimeDay,TimeSec,TimeMicroSec)
    # Position related
    lat=ma.masked_equal(dataset['Latitudes'][:],dataset['Latitudes']._FillValue)*dataset['Latitudes'].scale_factor
    lat.set_fill_value(dataset['Latitudes']._FillValue)
    lon=ma.masked_equal(dataset['Longitudes'][:],dataset['Longitudes']._FillValue)*dataset['Longitudes'].scale_factor
    lon.set_fill_value(dataset['Longitudes']._FillValue)
    #
    #data={'Cycles':cycles,'Tracks':tracks,'NbPoints':ndata,'Tracks Table':tracks_list,'TimeDay':TimeDay,'TimeSec':TimeSec,'TimeMicroSec':TimeMicroSec,'Time Table':time_list,'CorSSH':CorSSH,'Latitudes':lat,'Longitudes':lon,'MSS':MSS}
    #data={'Cycles':cycle_list,'Tracks':track_list,'TimeDay':TimeDay,'TimeSec':TimeSec,'TimeMicroSec':TimeMicroSec,'Datetime':time_list,'CorSSH':CorSSH,'Latitudes':lat,'Longitudes':lon,'MSS':MSS}
    data={'Cycles':cycle_list,'Tracks':track_list,'TimeDay':TimeDay,'TimeSec':TimeSec,'TimeMicroSec':TimeMicroSec,'Datetime':time_list,'Latitude':lat,'Longitude':lon}
    #
    for var in vars:
        tmp=ma.masked_equal(dataset[var][:,0],dataset[var]._FillValue)*dataset[var].scale_factor
        tmp.set_fill_value(dataset[var]._FillValue)
	data[var]=tmp
    return data

def make_SSHA(data):
    """
    """
    for c in data:
        for t in data[c]:
	    data[c][t]['SSHA'] = data[c][t]['CorSSH']-data[c][t]['MSS']
    return

def filter(data,var,limits):
    """

    ATENTION, change it to cond instead of limits, so give complete freedom for
      the conditions, like choose >= instead of > or only a lower limit.

      In work
    """
    index=(data[var].data>=limits[0])&(data[var].data<=limits[1])
    data_out={}
    for key in data:
        data_out[key]=data[key][index]
    return data_out

def load_TP_dataset(files,filtercond=None,data=None):
    """
    """
    if data is None:
        data={}
    elif type(data)!=dict:
        print "data should be a dictionary, and it's %s" % type(data)
        return
    if type(files)==str:
        fileslist = [files]
    else:
        fileslist = files
    i=0
    for file in fileslist:
        print "File: %s" % file
        try:
	    data_in = read_file(file)
	    if filtercond is not None:
	        for var in filtercond:
		  data_in=filter(data_in,var,filtercond[var])
            #
            for c in set(data_in['Cycles']):
	        print "Doing cycle: %s" % c
                if c not in data:
	            data[c]={}
                index_c = (data_in['Cycles'].data==c)
                for tck in set(data_in['Tracks'][index_c]):
		    #print "Doing track: %s" % tck
	            #if tck not in data_out[c].keys():
	            #    data_out[c][tck]={}
                    index_tck = index_c & (data_in['Tracks'].data==tck)
		    # Change it for a generic all keys
                    data[c][tck]={'Datetime':data_in['Datetime'][index_tck],'Latitude':data_in['Latitude'][index_tck],'Longitude':data_in['Longitude'][index_tck],'CorSSH':data_in['CorSSH'][index_tck],'MSS':data_in['MSS'][index_tck],'Bathy':data_in['Bathy'][index_tck]}
	    if i<=25:
	        i+=1
            else:
	        i=0
	        save_dataset(data,'load_TP_dataset.tmp')
        except:
            pass
    #
    return data

def load_from_path(path,filtercond=None):
    """

    Improve it to accept a URL too, in the case of a path for a DODS server.
    Maybe a regex too to restrict to nc files? Maybe to pattern of names.
    """
    import os
    filenames=os.listdir(path)
    filenames.sort()
    files=[os.path.join(path,filename) for filename in filenames]
    data=load_TP_dataset(files,filtercond)
    return data

def load_from_aviso(urlbase='ftp://ftp.cls.fr/pub/oceano/AVISO/SSH/monomission/dt/corssh/ref/j1/',filtercond=None):
    """ Load the data from aviso
    """
    import urllib
    import re
    import urlparse
    import StringIO
    import gzip
    import pupynere
    import tempfile

    f = urllib.urlopen(urlbase)
    content = f.read()
    filesnames = re.findall('CorSSH_Ref_\w{2}_Cycle\d{1,3}\.nc\.gz',content)
    data = {}
    for filename in filesnames[:3]:
        f = urllib.urlopen(urlparse.urljoin(urlbase,filename))
        #content = f.read()
        #f=StringIO.StringIO(content)
        f=StringIO.StringIO(f.read())
        zf = gzip.GzipFile(fileobj=f)
        f=open('tmp.nc','w')
        f.write(zf.read())
        f.close()
        data=topex.load_TP_dataset(['tmp.nc'],filtercond=filtercond,data=data)
        #x=NetCDFFile(tmp)
        #ncf=tempfile.mkstemp(text=zf.read())
        #unzf=tempfile.TemporaryFile()
        #unzf.write(zf.read())
        #unzf.seek(0)
        #ncf=pupy.NetCDFFile(fileobj=unzf)
        #print ncf.attributes
        #print ncf.variables.keys()

    return data
    

def save_dataset(data,filename):
    """
    """
    import pickle
    output = open(filename,'wb')
    pickle.dump(data, output)
    output.close()
    return

def load_dataset(filename):
    """
    """
    import pickle
    pkl_file = open(filename, 'rb')
    data = pickle.load(pkl_file)
    pkl_file.close()
    return data

def join_cycles(data):
    """Join all cycles, so that from data[c][t][vars] return data[t][vars] with all cycles
    """
    import numpy
    vars=data[data.keys()[0]][data[data.keys()[0]].keys()[0]].keys()
    data_out={}
    for t in invert_keys(data):
        data_out[t]={}
        for var in vars:
            data_out[t][var] = ma.masked_array([])

    for c in data:
        for t in data[c]:
            for var in data[c][t]:
	        data_out[t][var]=numpy.ma.concatenate((data_out[t][var],data[c][t][var]))

    return data_out


def invert_keys(data):
    """ Invert the hirerachy of the first 2 level keys in the dictionary.

        This is usable to group the data in tracks instead of cycles, like
	  data[tracks][cycles] = invert_keys(data[cycles][tracks])
    """
    data_out={}
    for c in data:
        for t in data[c]:
            if t not in data_out:
                data_out[t]={}
            if c not in data_out[t]:
                data_out[t][c]={}
            data_out[t][c] = data[c][t]
    return data_out

##############################################################################
#### Extras
##############################################################################

def make_L(data,direction='S',z=None,):
    """ Define the along track distance from one reference

        direction define the cardinal direction priority (N,S,W or E).
         S means that the reference will be the southern most point

        z define the bathymetry, if defined, the closest point to that
         bathymetry will be the reference. In case of cross this bathymetry
         more than once, the direction criteria is used to distinguish.
    """
    from fluid.common.distance import distance
    all_cycles_data = join_cycles(data)

    if z==None:
        import rpy
        #for t in topex.invert_keys(data):
        for t in all_cycles_data:
            rpy.set_default_mode(rpy.NO_CONVERSION)
            linear_model = rpy.r.lm(rpy.r("y ~ x"), data = rpy.r.data_frame(x=all_cycles_data[t]['Longitude'], y=all_cycles_data[t]['Latitude']))
            rpy.set_default_mode(rpy.BASIC_CONVERSION)
            coef=rpy.r.coef(linear_model)
            if direction=='S':
                lat0=all_cycles_data[t]['Latitude'].min()-1
                lon0 = (lat0-coef['(Intercept)'])/coef['x']
                L_correction = distance(all_cycles_data[t]['Latitude'],all_cycles_data[t]['Longitude'],lat0,lon0).min()
            for c in invert_keys(data)[t]:
                data[c][t]['L'] = distance(data[c][t]['Latitude'],data[c][t]['Longitude'],lat0,lon0)- L_correction
    # This bathymetric method was only copied from an old code. This should be atleast
    #  changed, if not removed.
    elif method=='bathymetric':
        import rpy
        for t in all_cycles_data:
            # First define the near coast values.
            idSouth=numpy.argmin(all_cycles_data[t]['Latitude'])
            L_tmp = distance(all_cycles_data[t]['Latitude'],all_cycles_data[t]['Longitude'],all_cycles_data[t]['Latitude'][idSouth],all_cycles_data[t]['Longitude'][idSouth])
            idNearCoast = L_tmp.data<400e3
            if min(all_cycles_data[t]['Bathy'][idNearCoast]) > -z:
                idNearCoast = L_tmp.data<600e3
            # Then calculate the distance to a reference
            rpy.set_default_mode(rpy.NO_CONVERSION)
            linear_model = rpy.r.lm(rpy.r("y ~ x"), data = rpy.r.data_frame(x=all_cycles_data[t]['Longitude'], y=all_cycles_data[t]['Latitude']))
            rpy.set_default_mode(rpy.BASIC_CONVERSION)
            coef=rpy.r.coef(linear_model)
            lat0 = all_cycles_data[t]['Latitude'].min()-1
            lon0 = (lat0-coef['(Intercept)'])/coef['x']
            #L = distance(,lon,lat0,lon0)
            #
            #id0 = numpy.argmin(numpy.absolute(all_cycles_data[t]['Bathy'][idNearCoast]))
            idref=numpy.argmin(numpy.absolute(all_cycles_data[t]['Bathy'][idNearCoast]+z))
            #L_correction = distance(all_cycles_data[t]['Latitude'][idNearCoast][idref],all_cycles_data[t]['Longitude'][idNearCoast][idref],all_cycles_data[t]['Latitude'][idNearCoast][idref],all_cycles_data[t]['Longitude'][idNearCoast][idref])
            L_correction = distance(all_cycles_data[t]['Latitude'][idNearCoast][idref],all_cycles_data[t]['Longitude'][idNearCoast][idref],lat0,lon0)
            for c in topex.invert_keys(data)[t]:
                #data[c][t]['L'] = distance(data[c][t]['Latitude'],data[c][t]['Longitude'],all_cycles_data[t]['Latitude'][idNearCoast][id0],all_cycles_data[t]['Longitude'][idNearCoast][id0]) - L_correction
                data[c][t]['L'] = distance(data[c][t]['Latitude'],data[c][t]['Longitude'],lat0,lon0) - L_correction
    #
    return


##############################################################################
#### 
##############################################################################

##############################################################################
class TOPEX(IterableUserDict):
    """
    """
    def __init__(self,tracks=None,nbpoints=None,Cycles=None):
        """
        """
        self.data={1:{11:[1,2,3],12:[1.1,1.2,1.3],'teste':['orange','apple','pear']},2:[10,20,30],3:[100,200,300]}
        return
    def __getitem__(self, key):
        print "Chave tipo %s" % type(key).__name__
        if isinstance(key, basestring):
            if key=='join':
	        print "key is join"
	    data_out={}
	    for k in self.data:
                if k not in data_out:
                    pass
        if isinstance(key, slice):
	    print "It's a slice"
	    if (key.start==0) & (key.stop>max(self.data.keys())+1) & (key.step==None):
	      return self.data
	    print "I'm not ready for that. Only full return, like x[:]"
	    return
        if key not in self.data:
	    print "%s is not a valid key" % key
	    return
            print "key: %s" % key
            return self.data[key]
