#!/bin/bash

#source ../header_guk6.sh
source ../header_guk7.sh
source ../header_theADCcode.sh

title="CO ndADC(3)-IP"

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

# note that ADC-IP can only compute doublets, so the spin keyword is inactive
spin="2" #spin

active="1 to 26"   # active orbital space

symgroup="C1"         # symmetry
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
ndadc3ip; spin $spin; sym $sym_nmbs
SYMGRP $symgroup
&self-energy
infinite
&diagonalizer
full
&eigen
ps 1.0
thresh 0.5
EOF
echo "#end<theADCcode>"

exit 0