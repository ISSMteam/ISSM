#!/bin/bash
set -eu

#Some cleanup
rm -rf install

#Download development version
git clone https://github.com/SciCompKL/MeDiPack.git install
