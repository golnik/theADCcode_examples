#!/bin/bash

source ../../header_FTprop.sh

ADCdir="$PWD/out_ADC"	
propdir="$PWD/out_prop"

prop_inp="propagator.inp"

mkdir -p $propdir
cd $propdir

cat << EOF > $prop_inp
; input file for FileTools propagation
[propagation]
abserr   = 1e-06           	; absolute error
dt       = 0.02            	; dt in fs
tstart   = 0              	; in dt intervals
tend     = 500            	; in dt intervals
tprint   = 1               	; print evolution and save files every tprint steps

[state]
vecfile     = $ADCdir/fmatrix_nd.1
get_columns = 12           ; get the specified columns from vector file
coeffs      = 1.0          ; initial state will be created by summing up the indicated columns with these coefficiets

[system]
H0 			= $ADCdir/ADC_File_dia.1,$ADCdir/ADC_File_off.1
basis       = $ADCdir/evecs.1
project     = false
states      = 1..4

[output]
prop_vecs   = $propdir/prop_vecs.1   ; propagated vectors will be saved in this file in FileTool's format
popout      = $propdir/popul.dat     ; evolution of populations will be saved to this file
EOF

time $FTprop_bin $prop_inp

exit 0


