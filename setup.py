import os
from setuptools import setup, find_packages, Extension

def read(fname):
	return open(os.path.join(os.path.dirname(__file__), fname)).read()

preprocess = Extension('waistcoat.preprocess', 
		sources=['waistcoat/preprocess.c',])

setup(
    name = "waistcoat",
    version = "0.01",
    packages = find_packages(),
		ext_modules = [preprocess,],

    install_requires = ['pysam>=0.7'],

		test_suite = 'test',

    # metadata for upload to PyPI
    author = "Haydn King",
    author_email = "hjk734@gmail.com",
    description = "A script to map RNA-seq data to a genome",
		long_description = read("README.md"),
    license = "GPLv2",

    # could also include long_description, download_url, classifiers, etc.
)
