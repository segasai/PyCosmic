#!/usr/bin/env python

from distutils.core import setup
setup(name='PyCosmic',
      version='0.5',
      description='Cosmic ray rejection algorithm on single frame CCD images',
      author='Bernd Husemann',
      author_email='husemann@mpia.de',
      packages=['PyCosmic','PyCosmic.lib'],
      url='http://sourceforge.net/projects/pycosmic/',
      requires=['scipy','numpy','astropy'],
      scripts=['PyCosmic/bin/PyCosmic'])
