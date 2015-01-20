.. PyAVISO documentation master file, created by
   sphinx-quickstart on Thu Nov 27 14:00:23 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyAVISO's documentation!
===================================

Contents:

.. toctree::
   :maxdepth: 2

aviso
-----

A small library to download data from AVISO. It uses the AVISO's DAP server, to
optimize the use of the net, downloading only the selected subset. For that,
you need to register a username and password at AVISO. It's free. The data is
downloaded in blocks, so don't overload AVISO's server, and run safe in unstable
networks.

Fast aviso howto
----------------
Run this in the shell, outside Python, to install the library:

pip install aviso

Now try this to download some data into a NetCDF file::

    AVISO_download -R -5/310/15/350 -D 1999-01-01/1999-06-01 \
    -u user_at_AVISO_DAP -p password_at_AVISO_DAP --map='madt+msla' \
    --timestep=30 -o aviso.nc

The longitudes can be defined as [0, 360] or [-180, 180], like::

    -R -5/-50/15/-10 ...

If -D is not defined, it will download all available timeseries.

It's not recommended to use shallow water altimetric data without a special treatment, so one can add into the command above the options::

    --maskshallowerthan=150 --etopofile=/data/ETOPO2v2c_f4.nc

Which means that all gridpoints shallower than 150m depth will be masked. 
The depth is defined from the ETOPO bathymetry. 
On this example, the file "ETOPO2v2c_f4.nc" should be at the directory: "/data/".

You can request a username to access AVISO's OpenDAP server on the website `http://www.aviso.altimetry.fr/en/data/data-access/aviso-opendap.html <http://www.aviso.altimetry.fr/en/data/data-access/aviso-opendap.html>`_

topex
-----
A very simple library to read the NetCDF altimetry data of TOPEX and JASON-1. The objective here is to automate the procedure of making the datasets ready to use. This include reading how all the required NetCDF files, sub-sampling by desired time/period and creating an accumulated intuitive dictionary of the pertinent data.

Although not yet implemented the most interesting feature would be a "TOPEX" class, which would make very simple and intuitive to load, sub-sample and deal with these datasets.


ATENTION
--------

This is not an official package, so please do not complain with AVISO if you have any trouble with it.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

