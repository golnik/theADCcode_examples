#!/bin/bash

#!!! IMPORTANT:
# note that a modified version of serenity code is used for these calculations
# see here: https://github.com/golnik/serenity, branch data_saver

source ../header_serenity.sh
source ../header_theADCcode.sh

# note that ADC-IP can only compute doublets, so the spin keyword is inactive
spin="2" #spin

active="1 19"   # active orbital space

# note that serenity can only perform calculations in C1 group!
symgroup="C1"         # symmetry
sym_nmbs="1"          # irreps

calc_dir=$PWD/out  # calc directory

mkdir -p $calc_dir
cd $calc_dir

cat << EOF > molecule.xyz
3
Water molecule
O          0.00000        0.00000        0.11779
H          0.00000        0.75545       -0.47116
H          0.00000       -0.75545       -0.47116
EOF

basisset="6-311G"

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