PyAVISO
=======

.. image:: https://readthedocs.org/projects/aviso/badge/?version=latest
   :target: https://readthedocs.org/projects/aviso/?badge=latest
      :alt: Documentation Status

.. image:: https://img.shields.io/travis/castelao/aviso.svg
        :target: https://travis-ci.org/castelao/aviso

.. image:: https://img.shields.io/pypi/v/aviso.svg
        :target: https://pypi.python.org/pypi/aviso


A library to handle altimetric data produced by AVISO.

* Free software: BSD license
* Documentation: https://pyaviso.readthedocs.org.

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

Now try this, with proper username and password, to download some data into a NetCDF file:

AVISO_download -R -5/310/15/350 -D 1999-01-01/1999-06-01 -u user_at_AVISO_DAP -p password_at_AVISO_DAP --map='msla' --timestep=30 -o aviso.nc

Documentation
-------------

The full documentation is at http://pyaviso.readthedocs.org

ATENTION
--------

This is not an official package, so please do not complain with AVISO if you have any trouble with it.
