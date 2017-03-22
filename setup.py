

import setuptools
from numpy.distutils.core import Extension, setup



classifiers = [
"Development Status :: 4 - Beta",
"Intended Audience :: Science/Research",
"License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
"Programming Language :: F",
"Programming Language :: Python :: 3.5",
"Topic :: Scientific/Engineering"]


spaic_module = Extension(name='spyclib.spaic',
            sources            = ['spyclib/spaic.f90'])



setup(
    name         = "spyclib",
    version      = "0.0.1",
    author       = "Justus Schwabedal",
    author_email = "JSchwabedal@gmail.com",
    description  = ("Solve the Schroedinger equation."),
    license      = "MIT",
    packages     = ['spyclib'],
    ext_modules	 = [spaic_module],
    classifiers	 = classifiers)
