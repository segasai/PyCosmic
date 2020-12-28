#!/usr/bin/env python3

from distutils.core import setup
setup(name='PyCosmic',
      version='0.6',
      description='Cosmic ray rejection algorithm on single frame CCD images',
      author='Bernd Husemann',
      author_email='husemann@gmx.de',
      packages=['PyCosmic'],
      url='http://sourceforge.net/projects/pycosmic/',
      requires=['scipy','numpy','astropy'],
      scripts=['bin/PyCosmic'])
