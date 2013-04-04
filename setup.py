from setuptools import setup, find_packages

classifiers = """\
Development Status :: 4 - Beta
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: MIT License 
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries :: Python Modules
"""

version = '0.3.0'

setup(name = 'pyaviso',
      version = version,
      description = "A very simple library to interpret and load TOPEX/JASON altimetry data",
      long_description = """\
A very simple library to read the NetCDF altimetry data of TOPEX and JASON-1. The objective here is automate the procedure of making the datasets ready to use. This include reading how all the required NetCDF files, sub-sampling by desired time/period and creating an accumulated intuitive dictionary of the pertinent data.

Although not yet implemented the most interesting feature would be a "TOPEX" class, which would make very simple and intuitive to load, sub-sample and deal with these datasets.""",
      classifiers=filter(None, classifiers.split("\n")),
      keywords='altimetry, TOPEX, JASON-1, oceanography, Sea Surface Height',
      author = 'Guilherme Castelao',
      author_email = 'guilherme@castelao.net',
      url = 'https://pypi.python.org/pypi/pyaviso',
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
      entry_points="""
      # -*- Entry points: -*-
      """,
     )



