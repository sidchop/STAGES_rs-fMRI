#!/bin/bash
#SBATCH --job-name=STAGES_fmriprep
#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
# SBATCH --mail-user=kristina.sabaroedin@monash.edu
# SBATCH --mail-type=BEGIN
# SBATCH --mail-type=FAIL
# SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=8000
#SBATCH --qos=normal
#SBATCH -A kg98

subject=$s


datadir=/home/scho0011/kg98/Sid/STAGES/STAGES_fmriprep/data
output_dir=$datadir/derivatives/fmriprep
fslicense=/home/scho0011/kg98/Sid/license.txt
work_dir=$datadir/fmriprep_work/$subject
bids_dir=$datadir/rawdata/

module load fmriprep/1.4.1
unset PYTHONPATH

# Run fmriprep
fmriprep \
--participant-label ${subject} \
--output-spaces func fsnative MNI152NLin2009cAsym:res-2 MNI152Lin:res-2 MNI152NLin6Asym:res-2 \
--mem_mb 80000 \
-n-cpus 8 \
--fs-license-file ${fslicense} \
--use-aroma \
--cifti-output \
--ignore fieldmaps \
--use-syn-sdc \
-w ${work_dir} \
${bids_dir} \
${output_dir} \
participant
