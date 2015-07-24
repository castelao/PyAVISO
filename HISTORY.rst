.. :changelog:

History
-------


0.7.3 - May-29-2013
-------------------

 * First prototype of a new class: Products. The command AVISO_download now can estimate the Okubo-Weiss parameter and relative vorticity after download from AVISO, and save it on the same NetCDF file

0.7.2
-----

* Migrating to netCDF4. It is based on HDF5, thus more robust system. Most important, it handles the append mode, important for the capability of estimate some products, to be included soon.

0.7.0 - May-21-2013
-------------------

* It now handles negative longitudes, as well as initial longitude bigger than the final one. And in that case it circles the Earth passing through Greenwich until the final longitude.

AVISO_download -R -15/-10/15/-60 -D 1999-01-01/1999-01-10 -u user_at_AVISO_DAP -p password_at_AVISO_DAP --timestep=5 -o aviso.nc

On the example ahead, it covers from 10W, pass through 0 and 180, until reaches 60W. While in the example bellow

AVISO_download -R -15/-60/15/-10 -D 1999-01-01/1999-01-10 -u user_at_AVISO_DAP -p password_at_AVISO_DAP --timestep=5 -o aviso.nc

it covers from 60W to 10W.

0.6.0 - May-14-2013
-------------------

AVISO_download now has the option -D, to handle dates, like
-D 1998-02-01/1999-01-01

0.4.1 - May-9-2013
------------------

* I'm refactoring the whole system. It will be organized with a fundamental class that downloads the data, and that's it. The other functionalities will be rebuild in different classes using this fundamental one.

AVISO_download now is based on the class AVISO_fetch.


Apr-9-2013
----------

* Initial prototype of bin/AVISO_download


Apr-4-2013
----------

I wrote pytopex somehere around 2006, during my PhD, while I was studing the North Brazil Current Rings, but I haven't working on that anymore for several years. topex.py is outdated, but still have some usefull functions if you need to work with the along track data.

aviso.py is a library to download data from AVISO. At this point it handles only absolute dynamic topography and sea surface anomaly, but could be easily adapted for other type of data. The download uses the AVISO's DAP server, and have a safe mode in case of slow and unreliable network. To avoid overload the server, it request the data in smaller blocks and if loose the connection in the process, it stand by and try again later from the point it stopped.
