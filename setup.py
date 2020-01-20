#!/usr/bin/env python
from __future__ import absolute_import
from setuptools import setup, find_packages

if __name__ == '__main__':
    setup(packages=find_packages(),
          long_description=open('README.md').read(),
          long_description_content_type='text/markdown',
          name="manage_crystal",
          author="Daniele Ongari",
          author_email="daniele.ongari@epfl.ch",
          description="Tools for manipulating crystal formats",
          url="https://github.com/danieleongari/manage_crystal",
          license="Creative Commons",
          classifiers=["Programming Language :: Python"],
          version="0.1.0",
          install_requires=["numpy"],
          scripts=[
              "manage_crystal/manage_crystal", "manage_crystal/trajec_crystal"
          ],
          extras_require={
              "pre-commit": [
                  "pre-commit==1.21.0",
                  "yapf==0.29.0",
              ]
          })
