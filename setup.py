#!/usr/bin/env python

from distutils.core import setup

setup(
    name='undr_rover',
    version='1.0.0',
    author='Roger Li',
    author_email='r.li3@student.unimelb.edu.au',
    packages=['undr_rover'],
    #entry_points={
    #    'console_scripts':['undr_rover = undr_rover.undr_rover:main']
    #},
    url='https://github.com/rli99/undr_rover',
    license='LICENSE.txt',
    description=(
        'UNDR ROVER: Unmapped primer directed read overlap variant caller.'),
    long_description=open('README.txt').read(),
    #install_requires=[
    #    "PyVCF",
    #    "biopython",
    #    "pyfaidx"
    #],
)
