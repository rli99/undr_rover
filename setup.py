#!/usr/bin/env python

from distutils.core import setup

setup(
    name='undr_rover',
    version='0.0.1',
    author='Roger Li',
    author_email='r.li3@student.unimelb.edu.au',
    packages=['undr_rover'],
    url='https://github.com/rli99/tbd',
    license='LICENSE.txt',
    description=(
        'UNDR ROVER: Unmapped primer directed read overlap variant caller.'),
    long_description=open('README.txt').read(),
    install_requires=[
        "PyVCF",
        "biopython",
        "pyfaidx"
    ],
)
