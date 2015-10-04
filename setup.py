#!/usr/bin/env python

'''
Created on Fri Aug 14 15:34:38 2015

@author: Dmitry A. Duev
'''

# compile fortran code using f2py
from numpy.distutils.core import Extension

import os

# fortran module to be compiled with f2py:
vintflib = Extension(name    = 'vintflib',
                     sources = ['src/pypride/vintflib.f'])

if __name__ == '__main__':
    # package name
    NAME = 'pypride'
    # executable command line scripts
    BIN = 'bin/'
    ls_bin = os.listdir(BIN)
    ls_bin = [os.path.join(BIN,l) for l in ls_bin if l[-3:]=='.py']
    
    from numpy.distutils.core import setup
    setup(name = NAME,
          description   = 'Python Tools for Planetary Interferometry and Doppler Experiments',
          author        = 'Dmitry A. Duev',
          author_email  = 'duev@jive.eu',
          url           = 'https://github.com/dmitryduev/pypride',
          platforms     = ['Linux', 'MacOS X'],
          licence       = 'GNU GPL v2',
          version       = '1.0.11',
          packages      = [NAME],
          package_dir   = {NAME: 'src/pypride'},
#          py_modules    = ['mod1', 'pkg.mod2'],
          ext_package   = NAME,
          ext_modules   = [vintflib],
          scripts       = ls_bin,
          package_data  = {NAME : ['cats/*', 'jpl_eph/*', 'inp.cfg']},
          install_requires = [ 'numpy',
                               'scipy',
                               'matplotlib',
                               'astropy',
                               'paramiko',
                               'numba',
                               'sklearn'],
          classifiers   = ['Development Status :: 5 - Production/Stable',
                           'Environment :: Console',
                           'Framework :: IDLE',
                           'Framework :: IPython',
                           'Intended Audience :: Science/Research',
                           'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
                           'Natural Language :: English',
                           'Operating System :: MacOS :: MacOS X',
                           'Operating System :: POSIX :: Linux',
                           'Programming Language :: Fortran',
                           'Programming Language :: Python',
                           'Programming Language :: Python :: 2.7',
                           'Programming Language :: Python :: 2 :: Only',
                           'Topic :: Scientific/Engineering :: Astronomy',
                           'Topic :: Scientific/Engineering :: Physics']
          )
