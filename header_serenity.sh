# This is an example. Content of this file could vary.
# But it must specify the command serenity_bin that runs serenity on your system.

module add gcc/12.3.0
module add hdf5/gcc12.3.0
module add boost/1.82.0/gcc8.5.0
module add eigen/3.4.0

root=$HOME/software/serenity/
serenity_bin=$root/build/bin/serenity

export SERENITY_HOME=$root
export SERENITY_RESOURCES=$SERENITY_HOME/data/
export SERENITY_BIN=$SERENITY_HOME/bin/