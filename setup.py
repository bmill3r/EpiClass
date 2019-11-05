'''
M
E
T
H
U ser
S earch for
E pigenetic
L ocus
A ssessment of
H eterogeneity

#########################\\>==+,,++=|\\################################
#######################\,......___,__.-.+\\############################
######################\+,__..___..._+=>>>==\\##########################
#####################\=,____,,++,,_...,=|\\<+>\########################
####################|=,,+++++,,,,+,,,,,,+=<\\\<\#######################
###################|=,,,,,,,,_____,,_,,,,,,>|\\=+<#####################
##################>|+,__...____......._,,,,+|\\|<\#####################
#################\=|+,_..__________....__,,+<\\\\\#####################
################\\>|+,__________.....______+<\\\\\#####################
###############<,_,\|+=====+,___,,,,,,____,<\\\\\\#####################
###############<+_+\|><\\\\\\\\\<\\\||\\||\\\\=,__\####################
###############\=,>\+,,++=<|\=<\+=<<===>\<+||+,..,>####################
################\++<=,,++++\=,_==,++++++=,,>++=__+\####################
#################\,__,,,++=,____,<>++++__,=+_++,+\#####################
##################\_.._..++,_.._,,=,...._,+_.__+#######################
###################\\+___,_________,...__,,_._|########################
#####################\,__,,,,,,,,,___.._,,=\|\#########################
######################|_____,,,_____..._,,<############################
#######################=_____________.___,|############################
########################=,,________,,____+\############################
########################\>,++,,,+,,,,___+\\<\##########################
########################\\>,,+,==___,,=\\\\\+>\\#######################
#####################\\\\\\|++,_,,_+>\\\\\\\>==,_,=|\\#################
################\\\<+<|\\\\\\\>+,++\#\\\\\\\|<<>++,,,,,,+=>|\\\########
##########\\<=+=|\<+<>\###\\\\\\\\\\\\\\\\\\\\|<>>>>>==+++++,++=<\#####

Optimizing and predicting performance of DNA methylation biomarkers using sequence methylation density information.

2019  Brendan F. Miller
bmille79 <at> jh <dot> edu

-----------
  PUBLIC DOMAIN NOTICE
 
  This software is "United States Government Work" under the terms of the United
  States Copyright Act. It was written as part of the authors' official duties
  for the United States Government and thus cannot be copyrighted. This software
  is freely available to the public for use without a copyright
  notice. Restrictions cannot be placed on its present or future use.
 
  Although all reasonable efforts have been taken to ensure the accuracy and
  reliability of the software and associated data, the National Human Genome
  Research Institute (NHGRI), National Institutes of Health (NIH) and the
  U.S. Government do not and cannot warrant the performance or results that may
  be obtained by using this software or data. NHGRI, NIH and the U.S. Government
  disclaim all warranties as to performance, merchantability or fitness for any
  particular purpose.
 
  Please cite the authors in any work or product based on this material.
-------

'''

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="methuselah",
    version="2.0.0-beta",
    author="Brendan F. Miller",
    author_email="bmille79@jh.edu",
    description="Optimizing and predicting performance of DNA methylation biomarkers using sequence methylation density information.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bmill3r/methuselah",
    packages=['methuselah'],
    package_dir={'methuselah': 'methuselah'},
    include_package_data=True,
    license='Public Domain',
    zip_safe=False,

    install_requires=[
           'numpy==1.16.*',
           'pandas==0.25.*',
           'scikit-learn==0.21.*',
           'scipy==1.3.*',
           'matplotlib==3.1.*',
           'tables==3.5.*'],

    entry_points={
        'console_scripts': [
            'methuselah = methuselah.__main__:main'
        ]
    },

    test_suite='nose.collector',
    tests_require=[
        'nose',
        'pytest',
        'pytest-cov',
        'coverage'],

    classifiers=[
        'Programming Language :: Python :: 3 :: Only',
        "Operating System :: OS Independent",
        'Development Status :: 4 - Beta',
        'License :: Public Domain',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Environment :: Console'
    ]
)