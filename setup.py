#!/usr/bin/env python

from distutils.core import setup
setup(name='PyCosmic',
      version='0.1',
      description='Cosmic ray rejection algorithm on single frame CCD images',
      author='Bernd Husemann',
      author_email='bhusemann@aip.de',
      packages=['PyCosmic'],
      scripts=['PyCosmic/bin/PyCosmic'])
