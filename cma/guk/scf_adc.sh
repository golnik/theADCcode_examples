#!/bin/bash

#source ../header_guk6.sh
source ../../header_guk7.sh
source ../../header_theADCcode.sh

title="propiolic ndADC(3)-IP"

geometry=$(cat <<-END
geometry angstroms
          0.0000000000          0.4869330000          0.0000000000 6.0 C
         -0.2294990000         -0.9388900000          0.0000000000 6.0 C
         -0.4780750000         -2.1108070000          0.0000000000 6.0 C
          1.3173440000          0.7815830000          0.0000000000 8.0 O
         -0.8733170000          1.3157160000          0.0000000000 8.0 O
         -0.6978930000         -3.1501090000          0.0000000000 1.0 H
          1.3911280000          1.7483030000          0.0000000000 1.0 H
end
END
)

basisset=$(cat <<-END
basis DZ
END
)

# note that ADC-IP can only compute doublets, so the spin keyword is inactive
spin="2" #spin

active="6 to 54"   # active orbital space

symgroup="C1"      # symmetry
sym_nmbs="1"       # irreps

calc_dir=$PWD/out_ADC  # calc directory

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

#
#############################################################################
#
# theADCcode

echo "#begin<theADCcode>"
$theADCcode_bin << EOF
&frontend
guk
&propagator
ndadc3ip; sym $sym_nmbs; spin $spin
ADC_matrix $calc_dir/ADC_File
f_matrix $calc_dir/fmatrix_nd
SYMGRP $symgroup
&self-energy
infinite
MAXSTA 67392
sig_matrix $calc_dir/sig_matrix
rho_matrix $calc_dir/rho_matrix
&savefiles
f_matrix
evecs
ADC_matrix
se_matrix
&diagonalizer
lanczos
iter 100
roots 1..10
evecs $calc_dir/evecs
&eigen
ps 0.0
thresh 0.1
EOF
echo "#end<theADCcode>"
