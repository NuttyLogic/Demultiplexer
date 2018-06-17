#!/usr/bin/env python3

import os
from setuptools import setup, find_packages
import unittest

# python3 setup.py sdist


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def my_test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='test_*.py')
    return test_suite


setup(name='Demultiplexer',
      version='1.0.0',
      description='Illumina Lane Demultiplexer, .qseq lane files to .fastq sample files',
      author='Colin Farrell',
      author_email='colinpfarrell@gmail.com',
      license='MIT',
      include_package_data=True,
      package_data={'': ['tests/tes_qseq/*.txt', 'tests/test_sample_files/*.txt', '*.txt']},
      packages=find_packages(),
      test_suite='setup.my_test_suite'
      )
