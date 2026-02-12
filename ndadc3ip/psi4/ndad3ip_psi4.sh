#!/bin/bash

source ../../header_theADCcode.sh

# note that ADC-IP can only compute doublets, so the spin keyword is inactive
spin="2" #spin

# active orbital space
i_act="1"
f_act="26"

# note that serenity can only perform calculations in C1 group!
symgroup="C1"         # symmetry
sym_nmbs="1"          # irreps

root_dir=$PWD
calc_dir=$root_dir/out_psi4  # calc directory

mkdir -p $calc_dir
cd $calc_dir

geometry=$(cat <<-END
  C  0.0  0.0  0.0
  O  0.0  0.0  1.128
END
)

basisset="6-311G"

###########################
# calculations start here #
###########################

# run psi4

echo "#begin<psi4>"

cat << EOF > script.py
import psi4

import sys
sys.path.append('$root_dir')
from hdf5_exporter import export_adc_hdf5

mol = psi4.geometry("""
$geometry
symmetry $symgroup
""")


psi4.set_options({
    "basis": "$basisset",
    "reference": "rhf",
    "scf_type": "pk",
    "e_convergence": 1e-10,
    "molden_write": False,
})

E, wfn = psi4.energy("scf", return_wfn=True)

export_adc_hdf5(
    wfn=wfn,
    mol=mol,
    filename="data.hdf5",
    iActOrb=$i_act,
    fActOrb=$f_act,
    eri_threshold=1.0e-10,
)
EOF

python script.py

echo "#end<psi4>"

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
