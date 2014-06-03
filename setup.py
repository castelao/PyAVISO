# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
#from distutils.core import setup

import os
import sys
from distutils import log

here = os.path.abspath(os.path.dirname(__file__))
README = "." #open(os.path.join(here, 'README.rst')).read()
NEWS = "." #open(os.path.join(here, 'NEWS.txt')).read()


classifiers = """\
Development Status :: 4 - Beta
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: MIT License 
Operating System :: OS Independent
Programming Language :: Python :: 2.6
Programming Language :: Python :: 2.7
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries :: Python Modules
"""

version = '0.7.8'

install_requires=[
    "numpy >= 1.1",
    "pupynere >= 1.0.15",
    "netCDF4 >= 1.0.4",
    "Pydap >= 3.1.RC1",
    #"fluid >= 0.1.10",
],
# How to include pupynere or netCDF4 ?


setup(
    name = 'AVISO',
    version = version,
    description = "A library to handle altimetric data produced by AVISO",
    long_description=README + '\n\n' + NEWS,
    classifiers=filter(None, classifiers.split("\n")),
    keywords='altimetry, AVISO, TOPEX, JASON, oceanography, Sea Surface Height',
    author = 'Guilherme Castelao',
    author_email = 'guilherme@castelao.net',
    url = 'https://pyaviso.castelao.net',
    #download_url = 'http://cheeseshop.python.org/packages/source/t/topex/topex-0.1.tar.gz',
    license = 'MIT',
    platforms = ['any'],
    py_modules=['topex','aviso'],
    zip_safe=True,
    # Fix it !!
    #install_requires=[
    #    # -*- Extra requirements: -*-
    #    'dap.plugins.netcdf',
    #],
    install_requires=install_requires,
    entry_points="""
    # -*- Entry points: -*-
    """,
    scripts=["bin/AVISO_download"],
    )



