#!/usr/bin/env python

from distutils.core import setup

setup(name='sdmpy',
      version='1.3',
      description='Python for ALMA/VLA Science Data Model',
      author='Paul Demorest',
      author_email='pdemores@nrao.edu',
      url='http://github.com/demorest/sdmpy',
      packages=['sdmpy'],
      package_data={'sdmpy':['xsd/*.xsd',]},
      scripts=['scripts/bdf_average.py']
     )
