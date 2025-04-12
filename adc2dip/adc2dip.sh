#!/bin/bash

#source ../header_guk6.sh
source ../header_guk7.sh
source ../header_theADCcode.sh

title="CO DIP-ADC(2)"

geometry=$(cat <<-END
geometry angstroms
  0.0 0.0 0.0     6.0   C
  0.0 0.0 1.128   8.0   O
end
END
)

basisset=$(cat <<-END
basis cc-pVDZ
END
)

# note that ADC2DIP can compute singlets and triplets
spin="1 3" #spins

active="3 to 30"   # active orbital space

symgroup="C2v"         # symmetry
sym_nmbs="1 2 3 4"     # irreps

calc_dir=$PWD/out  # calc directory

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
adc2dip; spin $spin; sym $sym_nmbs
SYMGRP $symgroup
&self-energy
fplus
&diagonalizer
full
&eigen
ps 5.0
thresh 0.5
EOF
echo "#end<theADCcode>"

exit 0