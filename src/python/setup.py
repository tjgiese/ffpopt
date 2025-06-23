#!/usr/bin/env python3

from setuptools import find_packages, setup
from glob import glob

install_requires = [
    'scipy',
    'numpy<2',
    'matplotlib',
    'ase',
    'geometric'
    ]

package_data = {
    "ffpopt": ["pkgdata/mace-off/*/*.model",
               "pkgdata/*/*.pb",
               "pkgdata/__init__.py"]
    }

scripts = glob("bin/*.py")

setup( name="ffpopt",
       version="0.1",
       description="force filed parameter optimzer",
       author="Timothy J. Giese",
       author_email="TimothyJGiese@gmail.com",
       platforms=["any"],
       license="MIT",
       url=None,
       python_requires='>3.5',
       install_requires=install_requires,
       include_package_data=True,
       package_data=package_data,
       scripts=scripts,
       packages=["ffpopt","ffpopt.pkgdata","ffpopt.ase","ffpopt.constants"],
       package_dir={"ffpopt": "./lib/ffpopt"} )

