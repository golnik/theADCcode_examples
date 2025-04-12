# theADCcode examples repository

This repository contains a set of examples demonstrating the usage of theADCcode program [see https://github.com/golnik/theADCcode] together with comparison of the results with other software (when available).

## Before start

Each directory contains a bash script that runs a frontend program to execute the HF calculations and then starts theADCcode with appropriate keywords. 
Since the way frontends programs and theADCcode are installed on a host system could vary, we separate the instructions for running these programs into corresponding header files. 
Before one can use theADCcode examples, an end user must create their own header files with the following names:

* __header_guk6.sh__: for running gamess-UK version 6.x
```
# This is an example. Content of this file could vary.
# But it must specify the command gamess_bin that runs guk v6.x on your system.
module add guk/6.2
gamess_bin=gamess
```

* __header_guk7.sh__: for running gamess-UK version 7.x
```
# This is an example. Content of this file could vary.
# But it must specify the command gamess_bin that runs guk v7.x on your system.
module add guk/7.0
gamess_bin=gamess
```

* __header_theADCcode.sh__: for running theADCcode
```
# This is an example. Content of this file could vary.
# But it must specify the command theADCcode_bin that runs theADCcode on your system.
module add gcc/12.3.0
module add mkl/2024.0
module add boost/1.82.0/gcc8.5.0
module add hdf5/intel2024.0.2
theADCcode_bin=$HOME/programs/ADC/theADCcode/build/theADCcode
```
