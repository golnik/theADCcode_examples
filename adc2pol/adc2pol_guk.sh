#!/bin/bash

#!!! IMPORTANT:
# note that guk v6.x performs atomic simulations in C1 group!
# make sure you specify symgroup="C1" and sym_nmbs="1" if you use guk v6.x!

#source ../header_guk6.sh
source ../header_guk7.sh
source ../header_theADCcode.sh

title="Ar ADC(2) excited states calculations"

geometry=$(cat <<-END
zmat angstroms
  Ar
end
END
)

basisset=$(cat <<-END
basis 6-31G
END
)

# note that ADC2POL can only compute singlets, so the spin keyword is inactive
spin="1"  #spins

active="1 to 13"   # active orbital space

symgroup="C1"                 # symmetry
sym_nmbs="1"     # irreps

calc_dir=$PWD/out_guk  # calc directory

mkdir -p $calc_dir
cd $calc_dir

###########################
# calculations start here #
###########################

if [ $symgroup == 'C1' ]; then
  nosym="adapt off"
fi

# gamess-environment:
export ed3=dfile
export ed6=vfile

# SCF-Calculation:

echo "#begin<scf>"
$gamess_bin << EOF
title
$title
charge 0
$geometry
$nosym
integral high
form high
$basisset
SCFTYPE MP2
Threshold 8
vectors extguess
enter 1
EOF
echo "#end<scf>"

# SCF-Transformation

echo "#begin<trans>"
$gamess_bin << EOF
restart
title
$title
charge 0
$geometry
$nosym
integral high
$basisset
bypass hf
runtype transform
active
$active
end
threshold 8
vectors 1
enter 1
EOF
echo "#end<trans>"

# theADCcode

echo "#begin<theADCcode>"
$theADCcode_bin << EOF
&frontend
guk
&propagator
adc2pol; spin $spin; sym $sym_nmbs
SYMGRP $symgroup
&self-energy
infinite
&diagonalizer
lanczos
iter 20
&eigen
ps 5.0
thresh 0.5
EOF
echo "#end<theADCcode>"

exit 0