#!/usr/bin/env python

import os,sys
if sys.version_info[:2] < (2, 7):
    raise Exception('This version of gensim needs Python 2.7 or later. ')
elif sys.version_info[:2] > (3, 0):
    raise Exception('This version of gensim might have errors in Python 3 or later. ')

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="NanoMod", # Replace with your own username
    version="0.2.2",
    author="Qian Liu",
    author_email="liuqianhn@gmail.com",
    description="A computational tool to detect DNA modifications using Nanopore long-read sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/WGLab/NanoMod",
    #packages=setuptools.find_packages(),
    packages=['scripts'],
    package_dir={'scripts': 'bin/scripts'},
    classifiers=[
        "Programming Language :: Python :: 2.7",
        'Intended Audience :: Science/Research',
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: OS Independent",
    ],
)
