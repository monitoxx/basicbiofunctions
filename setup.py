from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'My first python package'
LONG_DESCRIPTION = ('My first python package has the intention of providing ease of use to junior bioinformatics'
                    'for handling, processing and plotting large amounts of sequences.')

setup(
    name="basicbiofunctions",
    version=VERSION,
    author="Luis Javier Jimenez Bernal",
    author_email="<ljimenezbe@unal.edu.co>",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=[re, random, pandas],

    keywords=['python', 'first package'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)