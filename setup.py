#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as f:
    README = f.read()

setup(
    name='tb-rnap-compensation',
    version='0.1.0',
    description='Infer if a non-resistant mutation is associated with a resistance mutation using Fishers Exact Test',
    author='Viktoria Brunner',
    author_email='viktoria.brunner@dtc.ox.ac.uk',
    url='https://github.com/fowler-lab/tb-rnap-compensation',
    scripts=['bin/calculate-fisher-tests.py','bin/results-evaluation.py'],
    long_description = README,
    install_requires=[
        'pandas',
        'fisher',
        'tqdm',
        #'xlsxwriter'
        ],
    packages = ['tb_rnap_compensation'],
    package_data={'': ['tables/*']},
    python_requires='>=3.8',
    zip_safe=False
    )
