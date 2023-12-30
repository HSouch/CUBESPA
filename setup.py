from distutils.core import setup
from setuptools import find_packages
import cubespa

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setup(
    name='cubespa',
    packages=find_packages(),
    version=cubespa.__version__,
    license='BSD-3-Clause',
    description = 'Cube SPectral Analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Harrison Souchereau",
    author_email='harrison.souchereau@yale.edu',
    url='https://github.com/HSouch/cubespa',
    keywords='galaxies datacubes spectra',
    scripts=[
        
    ],
    install_requires=[
        # 'numpy',
        # 'photutils',
        # 'astropy',
        # 'matplotlib',
        # "reproject"
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.9',
    ],
)
