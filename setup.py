from setuptools import setup, find_namespace_packages

import pkg_resources
#include a pyproject.toml file?

setup(
    name='LR2TF-py',
    version='0.0.1',
    install_requires=[
    ],
    packages=find_namespace_packages(where="scr"),
    package_dir={"": "scr"},
    package_data={"LR2TF_py.data": ["*.csv"]},
    include_package_data=True
    
)
