#!/bin/bash

#!!! IMPORTANT:
# note that a modified version of serenity code is used for these calculations
# see here: https://github.com/golnik/serenity, branch data_saver

source ../header_serenity.sh
source ../header_theADCcode.sh

spin="1" #spin

active="1 13"   # active orbital space

# note that serenity can only perform calculations in C1 group!
symgroup="C1"         # symmetry
sym_nmbs="1"          # irreps

calc_dir=$PWD/out_serenity  # calc directory

mkdir -p $calc_dir
cd $calc_dir

cat << EOF > molecule.xyz
1
carbon mono oxide
Ar 0.0 0.0 0.0
EOF

basisset="6-31G"

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
