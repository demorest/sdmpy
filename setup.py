#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name = 'sdmpy',
    version='1.36',
    description='Python for ALMA/VLA Science Data Model',
    author='Paul Demorest',
    author_email='pdemores@nrao.edu',
    url='http://github.com/demorest/sdmpy',
    packages = find_packages(),        # get all python scripts in realtime
    install_requires=['lxml', 'numpy'],
    package_data={'sdmpy':['xsd/*.xsd',]},
    scripts=['scripts/bdf_average.py','scripts/bdf_bin_split.py']
)
