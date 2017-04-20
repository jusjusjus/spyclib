

import setuptools
from numpy.distutils.core import Extension, setup



classifiers = [
"Development Status :: 4 - Beta",
"Intended Audience :: Science/Research",
"License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
"Programming Language :: F",
"Programming Language :: Python :: 3.5",
"Topic :: Scientific/Engineering"]


spaic_module = Extension(
        name    = 'spyclib.spaic',
        sources = [
            'spyclib/spaic.f90'
        ],
        # extra_f90_compile_args = ["-fopenmp"],
        # extra_link_args = ["-lgomp"]
)


spaic2_module = Extension(
        name    = 'spyclib.spaic2',
        sources = [
            'spyclib/spaic2_d2pot.f90',
            'spyclib/spaic2_betaSpect.f90'
        ],
        # extra_f90_compile_args = ["-fopenmp"],
        extra_link_args = ["-llapack", "-lblas"]
)



setup(
    name         = "spyclib",
    version      = "0.0.1",
    author       = "Justus Schwabedal",
    author_email = "JSchwabedal@gmail.com",
    description  = ("Solve the Schroedinger equation."),
    license      = "MIT",
    packages     = ['spyclib'],
    ext_modules	 = [spaic_module, spaic2_module],
    classifiers	 = classifiers)
