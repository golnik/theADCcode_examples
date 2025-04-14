#!/bin/bash

source ../../header_theADCcode.sh

# note that ADC-IP can only compute doublets, so the spin keyword is inactive
spin="2" #spin

symgroup="C1"      # symmetry
sym_nmbs="1"       # irreps

cdir=$PWD                         # current directory
ADCdir="$cdir/out_ADC"                # folder with ADC calculations
propdir="$cdir/out_prop"              # folder with propagated WF
cmadir="$cdir/out_cma"		 	     # cma output folder

mkdir -p $cmadir
cp $ADCdir/data.hdf5 $cmadir
cd $cmadir

###########################
# calculations start here #
###########################

echo "#begin<theADCcode>"
$theADCcode_bin << EOF
&frontend
hdf5
&propagator
adcload; sym $sym_nmbs; spin $spin
ADC_matrix $ADCdir/ADC_File
f_matrix $ADCdir/fmatrix_nd
SYMGRP $symgroup
&self-energy
load
sig_matrix $ADCdir/sig_matrix
rho_matrix $ADCdir/rho_matrix
&cma
psi
prop_vecs $propdir/prop_vecs
z_axis y
z_min -7.0
z_max 7.0
z_npts 100
Q_file Q.dat
EOF
echo "#end<theADCcode>"

#plot Q density
python $cdir/../plotQ.py $cmadir/Q.dat.1 $cdir/Q.png