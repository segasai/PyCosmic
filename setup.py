#!/usr/bin/env python3

from distutils.core import setup
setup(name='PyCosmic',
      version='0.6',
      description='Cosmic ray rejection algorithm on single frame CCD images',
      author='Bernd Husemann',
      author_email='berndhusemann@gmx.de',
      classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3 :: Only',
            'Topic :: Scientific/Engineering :: Astronomy',
            'Topic :: Scientific/Engineering :: Image Processing'
      ],
      packages=['PyCosmic'],
      url='http://sourceforge.net/projects/pycosmic/',
      install_requires=['scipy>=1.3.0','numpy>=1.17','astropy>=3.0'],
      python_requires='>=3.5',
      scripts=['bin/PyCosmic'])
