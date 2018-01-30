#!/usr/bin/env python3

import os
from setuptools import setup, find_packages
import unittest

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# python3 setup.py sdist


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def my_test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='test_*.py')
    return test_suite

setup(name='Demultiplexer',
      version='0.4.0',
      description='Demultiplex Illumina lane .qseq files to sample .fastq files',
      author='Colin Farrell',
      author_email='colinpfarrell@gmail.com',
      license='MIT',
      include_package_data=True,
      package_data={'': ['tests/tes_qseq/*.txt', 'tests/tes_sample_files/*.txt', '*.txt']},
      packages=find_packages(),
      test_suite='setup.my_test_suite'
      )