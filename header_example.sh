#!/bin/bash

# gamess-uk
#module add guk/6.2
#or
module add guk/7.0

gamess_bin=gamess

# serenity


# theADCcode
module add gcc/12.3.0
module add mkl/2024.0
module add boost/1.82.0/gcc8.5.0
module add hdf5/intel2024.0.2

theADCcode_bin=$HOME/programs/ADC/theADCcode/build/theADCcode