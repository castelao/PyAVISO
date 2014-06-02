PyAVISO
=======

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

Now try this to download some data into a NetCDF file:

AVISO_download -R -5/310/15/350 -D 1999-01-01/1999-06-01 -u user_at_AVISO_DAP -p password_at_AVISO_DAP --map='madt+msla' --timestep=30 -o aviso.nc

topex
-----
A very simple library to read the NetCDF altimetry data of TOPEX and JASON-1. The objective here is to automate the procedure of making the datasets ready to use. This include reading how all the required NetCDF files, sub-sampling by desired time/period and creating an accumulated intuitive dictionary of the pertinent data.

Although not yet implemented the most interesting feature would be a "TOPEX" class, which would make very simple and intuitive to load, sub-sample and deal with these datasets.


ATENTION
--------

This is not an official package, so please do not complain with AVISO if you have any trouble with it.
