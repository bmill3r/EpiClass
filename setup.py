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

Optimizing and predicting performance of DNA methylation biomarkers using methylation density information.

Copyright (C) 2019  Brendan F. Miller
bmille79 <at> jh <dot> edu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.

'''

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="methuselah",
    version="2.0.0-beta",
    author="Brendan F. Miller",
    author_email="bmille79@jh.edu",
    description="Optimizing and predicting performance of DNA methylation biomarkers using methylation density information.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bmill3r/methuselah",
    packages=['methuselah'],
    package_dir={'methuselah': 'methuselah'},
    include_package_data=True,
    license='Public Domain',
    zip_safe=False,

    install_requires=[
           'numpy=1.16.*',
           'pandas=0.25.*',
           'scikit-learn=0.21.*',
           'scipy=1.3.*',
           'matplotlib=3.1.*',
           'pytables=3.5.*'],

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
