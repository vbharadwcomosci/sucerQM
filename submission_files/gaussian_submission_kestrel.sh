#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --partition=standard
#SBATCH --job-name=opt
#SBATCH --output=std.out
#SBATCH --error=std.err
#SBATCH --account=rg2
#SBATCH --exclusive
#----------------------------------------------------------#
module load gaussian
module list
cd $SLURM_SUBMIT_DIR
FULL_NAME=$(find -type f -name "*.com")
TEMP_BASENAME="${FULL_NAME:2}"
INPUT_BASENAME="${TEMP_BASENAME%.*}"
echo $INPUT_BASENAME
GAUSSIAN_EXEC=g16
if [ -e /dev/nvme0n1 ]; then
echo 'This node has a local storage and will use as the scratch path'
SCRATCH=$TMPDIR
else
echo 'This node does not have a local storage drive and will use /scratch as the scratch path'
SCRATCH=/scratch/$USER/$SLURM_JOB_ID
fi
mkdir -p $SCRATCH
export GAUSS_SCRDIR=$SCRATCH
# Run gaussian NREL script
g16_nrel

#Setup Linda paprameters
if [ $SLURM_JOB_NUM_NODES -gt 1 ]; then 
export GAUSS_LFLAGS='-vv -opt 'Tsnet.Node.lindarsharg: ssh'' 
export GAUSS_EXEDIR=$g16root/g16/linda-exe:$GAUSS_EXEDIR
fi
####################
# Run Gaussian job #
####################
srun $GAUSSIAN_EXEC < $INPUT_BASENAME.com >& $INPUT_BASENAME.log
echo Optimization is done


if [ -d $SCRATCH ]
then
   rm -rf $SCRATCH
fi
