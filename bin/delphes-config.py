#!/usr/bin/env python
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--cflags', action='store_true', help='Help text')
parser.add_argument('--incdir', action='store_true', help='Help text')
parser.add_argument('--libs', action='store_true', help='Help text')
parser.add_argument('--ldflags', action='store_true', help='Help text')
args = parser.parse_args()

# prefix = os.path.realpath(os.path.dirname(__file__))
conda_prefix = os.getenv('CONDA_PREFIX')
if conda_prefix is None:
    raise EnvironmentError('CONDA_PREFIX is not set. Please activate the conda environment.')

include_dir = os.path.join(conda_prefix, "include")
lib_dir = os.path.join(conda_prefix, "lib")

if args.cflags:
    print(f'-std=c++17 -I{include_dir}')

if args.incdir:
    print(f'{include_dir}')

if args.libs:
    print(f'-L{lib_dir} -lDelphes -Wl,-rpath,{lib_dir} -pthread -lm -ldl -rdynamic')

if args.ldflags:
    print()
