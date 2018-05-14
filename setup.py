#!/usr/bin/env python

'''
setup.py file for mpscore2
'''
from setuptools import setup, find_packages
version = open('VERSION').readline().strip()

setup(
    name             = 'mpscore2',
    version          = version,
    description      = '''Proteomics post-search algorithm''',
    author           = 'Mark Ivanov & Lev Levitsky & Julia Bubis',
    author_email     = 'pyteomics@googlegroups.com',
    url              = 'https://bitbucket.org/markmipt/mpscore2',
    packages         = ['mpscore2', ],
    install_requires = [line.strip() for line in open('requirements.txt')],
    classifiers      = ['Intended Audience :: Science/Research',
                        'Programming Language :: Python :: 2.7',
                        'Programming Language :: Python :: 3',
                        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    license          = 'License :: OSI Approved :: Apache Software License',
    entry_points     = {'console_scripts': ['mpscore2 = mpscore2.search:run', ]}
    )
