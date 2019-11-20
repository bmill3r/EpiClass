#!/bin/bash

git clone https://github.com/bmill3r/EpiClass

cd EpiClass

conda env create --file epiclass_env.yml
conda activate epiclass

python setup.py build; python setup.py install

epiclass -V