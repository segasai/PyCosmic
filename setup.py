#!/usr/bin/env python

from distutils.core import setup
setup(name='PyCosmic',
      version='0.3',
      description='Cosmic ray rejection algorithm on single frame CCD images',
      author='Bernd Husemann',
      author_email='bhuseman@eso.org',
      packages=['PyCosmic','PyCosmic.lib'],
      url='http://sourceforge.net/projects/pycosmic/',
      requires=['scipy','numpy','pyfits'],
      scripts=['PyCosmic/bin/PyCosmic'])
