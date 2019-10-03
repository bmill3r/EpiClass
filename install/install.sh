#!/bin/bash

git clone https://github.com/bmill3r/methuselah

conda env create --file methuselah/environment.yml
conda activate methuselah

cd methuselah; python setup.py build; python setup.py install

methuselah -V