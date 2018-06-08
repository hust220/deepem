source env.sh

####################### PLEASE EDIT THE THREE VARIABLES BLOW!!! ######################
MPI=true # true: use openmpi for paralleling, false: use openmp for paralleling

COMPILATION_THREADS=8 # set the number of threads for compilation

DEBUG=false # true: debug mode, false: release mode

INFO=0 # 0: no information printed, 1: information of root node, 2: information of all nodes

RAND=true
######################################################################################

make -j${COMPILATION_THREADS} DEBUG=${DEBUG} INFO=${INFO} MPI=${MPI} RAND=${RAND}



