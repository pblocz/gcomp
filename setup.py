#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import getpass
import codecs
from setuptools import setup
from setuptools.command.install import install

setup(
    name='gcomp',
    version="0.1",
    license='BSD New',
    author='Pablo Cabeza',
    author_email='lemniscata.lmn@gmail.com',
    description='auxiliary library for computational geometry',
    # packages=['gcomp'],

    install_requires=[
        'scipy',
        'numpy',
        'sympy',
        'matplotlib',
        'mayavi',
    ],

)
