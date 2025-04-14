# theADCcode examples repository

This repository contains a set of examples demonstrating the usage of theADCcode program (see https://github.com/golnik/theADCcode) together with comparison of the results with other software (when available).

## Before Start

Each directory contains a bash script that runs a frontend program to execute the HF calculations and then starts theADCcode with appropriate keywords. 
Since the way frontends programs and theADCcode are installed on a host system could vary, we separate the instructions for running these programs into corresponding header files. 
Before one can use theADCcode examples, an end user should modify the corresponding header files to make the required codes working.

Header files correspond to the following packages:

* __header_theADCcode.sh__: for running theADCcode
* __header_guk6.sh__: for running gamess-UK version 6.x
* __header_guk7.sh__: for running gamess-UK version 7.x
* __header_serenity.sh__: for running serenity code. Note that a modified version of the serenity is required (https://github.com/golnik/serenity, branch data_saver).
* __header_FTprop.sh__: for running FileToolsPropagator (this is required to run cma dynamics example, see [here](https://bitbucket.org/golnik/theadccode/wiki/FileTool's%20propagation) for more details)
