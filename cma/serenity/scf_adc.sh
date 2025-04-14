#!/bin/bash

#!!! IMPORTANT:
# note that a modified version of serenity code is used for these calculations
# see here: https://github.com/golnik/serenity, branch data_saver

source ../../header_serenity.sh
source ../../header_theADCcode.sh

# note that ADC-IP can only compute doublets, so the spin keyword is inactive
spin="2" #spin

active="6 54"   # active orbital space

# note that serenity can only perform calculations in C1 group!
symgroup="C1"      # symmetry
sym_nmbs="1"       # irreps

calc_dir=$PWD/out_ADC  # calc directory

mkdir -p $calc_dir
cd $calc_dir

cat << EOF > molecule.xyz
7
propiolic acid
C          0.0000000000          0.4869330000          0.0000000000
C         -0.2294990000         -0.9388900000          0.0000000000
C         -0.4780750000         -2.1108070000          0.0000000000
O          1.3173440000          0.7815830000          0.0000000000
O         -0.8733170000          1.3157160000          0.0000000000
H         -0.6978930000         -3.1501090000          0.0000000000
H          1.3911280000          1.7483030000          0.0000000000
EOF

basisset="DZ"

###########################
# calculations start here #
###########################

# run serenity

echo "#begin<serenity>"

cat << EOF > input.inp
+system
  name molecule
  geometry ./molecule.xyz
  method hf
  +basis
    label $basisset
    makeSphericalBasis false
    integralThreshold 1E-10
  -basis
-system
+task scf
  system molecule
-task
+task SAVEDATA
  system molecule
  activeOrbs {$active}
-task
EOF

$serenity_bin input.inp

echo "#end<serenity>"

# theADCcode

echo "#begin<theADCcode>"
$theADCcode_bin << EOF
&frontend
hdf5
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
