import os
from setuptools import setup, find_packages

def read(fname):
	return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "waistcoat",
    version = "0.01",
    packages = find_packages(),

    install_requires = ['docutils>=0.3'],

		test_suite = 'tests'

    # metadata for upload to PyPI
    author = "Haydn King",
    author_email = "hjk734@gmail.com",
    description = "A script to map RNA-seq data to a genome",
		long_description = read("README")
    license = "GPLv2",

    # could also include long_description, download_url, classifiers, etc.
)
