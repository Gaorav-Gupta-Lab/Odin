"""
Setup file for Odin
"""
import sys
from setuptools import setup, find_namespace_packages

if sys.version_info < (3, 5):
    sys.stdout.write("At least Python 3.5 is required.\n")
    sys.exit(1)

setup(
    name='Odin',
    version='1.0.0',
    packages=find_namespace_packages(include=['Valkyries.*', 'odin.*']),

    url='',
    license='',
    author='Dennis A. Simpson',
    author_email='dennis.simpson@rtpgenomics.com',
    long_description=open('README.md').read(),

    description='Package for detection of residual disease in AML patients',
    install_requires=['pathos', 'sortedcontainers', 'natsort', 'python-magic', 'pysam', 'pyfaidx', 'pyensembl', 'numpy',
                      'pybedtools', 'pybedtools', 'scipy', 'cython', 'setuptools', ' PyVCF', 'python-Levenshtein',
                      'more-itertools', 'pandas']
    )
