PyAVISO
=======

aviso
-----

A small library to download data from AVISO. It uses the AVISO's DAP server, so you
  use just the necessary of your link, downloading only the selected subset. For that,
  you need to register a username and password at AVISO. It's free. The data is
  downloaded in blocks, so don't overload AVISO's server, and run safe in unstable
  networks.

Fast aviso howto
----------------
Just run in the shell

AVISO_download -R -5/310/15/350 -u user_at_AVISO_DAP -p password_at_AVISO_DAO --timestep=100

topex
-----
A very simple library to read the NetCDF altimetry data of TOPEX and JASON-1. The objective here is automate the procedure of making the datasets ready to use. This include reading how all the required NetCDF files, sub-sampling by desired time/period and creating an accumulated intuitive dictionary of the pertinent data.

Although not yet implemented the most interesting feature would be a "TOPEX" class, which would make very simple and intuitive to load, sub-sample and deal with these datasets.
