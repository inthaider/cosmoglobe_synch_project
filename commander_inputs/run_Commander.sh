#!/usr/bin/bash

killall -9 commander3


export OPTYPE=sample                        # {sample,optimize}
export SYNCH_BETA_LNL=chisq                 # log-likelihood type {prior,chisq,ridge,marginal}. All evaluated at smoothing scale, local sampling
export MD_VERSION=v02                       # version of md*.dat file used for run
export NGIBBS_ITER=100                       # Length of each Markov chain
export OUTPUT_DIR_COMM=/mn/stornext/d16/cmbco/AST9240/2022/jibran/commander_chains/${OPTYPE}11_beta${SYNCH_BETA_LNL}_md${MD_VERSION}_bd005017_cp020406_c1_k${NGIBBS_ITER}
rm -rf $OUTPUT_DIR_COMM; mkdir -p $OUTPUT_DIR_COMM

PARAM_TEMP_FILE='/mn/stornext/d16/cmbco/AST9240/2022/jibran/commander_inputs/tmpenv_param_jibran_LWA1.txt'
PARAM_FILE=/mn/stornext/d16/cmbco/AST9240/2022/jibran/commander_inputs/params/param_jibran_LWA1_md${MD_VERSION}.txt

MYVARS='$OPTYPE:$SYNCH_BETA_LNL:$MD_VERSION:$NGIBBS_ITER:$OUTPUT_DIR_COMM'
envsubst "$MYVARS" < $PARAM_TEMP_FILE > $PARAM_FILE
# printf "%s" "$(<$PARAM_FILE)"


NUM_PROC=64
HEALPIX_ROOT='/mn/stornext/u3/jibran/Commander/build_owl2528_oneapi/install/healpix'
COMM3_BIN='/mn/stornext/u3/jibran/Commander/build_owl2528_oneapi/install/bin/commander3'
# CHAINS_DIR='/mn/stornext/d16/cmbco/AST9240/2022/jibran/commander_chains/'
CHAINS_DIR=$OUTPUT_DIR_COMM
OUT_FILE='slurm.txt'
COMMANDER_PARAMS_DEFAULT_ROOT='/mn/stornext/u3/jibran/Commander/commander3/parameter_files/defaults'

export OMP_NUM_THREADS=1
export HEALPIX=$HEALPIX_ROOT
export COMMANDER_PARAMS_DEFAULT=$COMMANDER_PARAMS_DEFAULT_ROOT
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'/mn/stornext/u3/jibran/Commander/build_owl2528_oneapi/lib':'/mn/stornext/u3//jibran/Commander/build_owl2528_oneapi/healpix/libexport'



mpirun -env I_MPI_FABRICS=shm -env I_MPI_PIN_DOMAIN=numa -np $NUM_PROC $COMM3_BIN $PARAM_FILE 2>&1 | tee $CHAINS_DIR/$OUT_FILE



# mpirun -env I_MPI_FABRICS shm -np $NUM_PROC $COMM3_BIN $PARAM_FILE 2>&1 | tee $CHAINS_DIR/$OUT_FILE

# build=owl2528
# pfile=param_WMAP_only.txt
# dir=chains_WMAP_220825
# mkdir -p $dir
# cp $pfile $dir/param_orig.txt
# mpiexec -env I_MPI_FABRICS shm -n $n /mn/stornext/u3/jibran/Commander/build_$build"_oneapi/install/bin/commander3" $pfile --OUTPUT_DIRECTORY=$dir 2>&1| tee $dir/slurm.txt

# PARAM_FILE=/mn/stornext/u3/jibran/Commander/commander3/parameter_files/param_BP10_final_v2.txt
# CHAINS_DIR=/mn/stornext/u3/jibran/Commander_Outputs/BP10_rep/chains