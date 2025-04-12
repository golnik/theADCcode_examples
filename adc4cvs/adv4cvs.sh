#!/bin/bash

#!!! IMPORTANT:
# note that at the moment adc4 works with guk v6.x only!

source ../header_guk6.sh
source ../header_theADCcode.sh

title="NH3 CVS-ADC(4)-IP"

geometry=$(cat <<-END
zmat angstrom
H  
N 1 hn  
H 2 hn 1 hnh  
H 2 hn 1 hnh 3 113.7187      
constants
hn 1.012 
hnh 106.67 
end
END
)

basisset=$(cat <<-END
basis DZ
END
)

# note that ADC-IP can only compute doublets, so the spin keyword is inactive
spin="2" #spin

active="1 to 16"   # active orbital space

symgroup="Cs"      # symmetry
sym_nmbs="1"       # irreps

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
adc4cvs; spin $spin; sym $sym_nmbs
SYMGRP $symgroup
&self-energy
infinite
&adc4
MCORE 1
&diagonalizer
davi
nroots 55
maxdavsp 900
iter 30 
memory 1779900918
&eigen
ps 1.0
thresh 0.1
EOF
echo "#end<theADCcode>"

exit 0