#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Note: To use the 'upload' functionality of this file, you must:
#   $ pip install twine

from setuptools import setup, find_packages

setup(
    name='sdmpy',
    version='1.71.0',
    description='Python for ALMA/VLA Science Data Model',
    author='Paul Demorest',
    author_email='pdemores@nrao.edu',
    url='http://github.com/demorest/sdmpy',
    packages=find_packages(),        # get all python scripts in realtime
    install_requires=['lxml', 'numpy', 'future', 'progressbar2'],
    package_data={'sdmpy': ['xsd/*.xsd']},
    scripts=['scripts/bdf_average.py',
             'scripts/bdf_bin_split.py',
             'scripts/bdf_dedisperse_cut.py',
             'scripts/bdf_to_psrfits.py',
             'scripts/realfast_sdm_fix.py',
             ]
)
