#!/bin/bash
#SBATCH --job-name=STAGES_mriqc
#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1-00:00:00
# SBATCH --mail-user=kristina.sabaroedin@monash.edu
# SBATCH --mail-type=BEGIN
# SBATCH --mail-type=FAIL
# SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=8000
#SBATCH --qos=normal
#SBATCH -A kg98

module load mriqc/0.9.7
subject=$s
mriqc /home/scho0011/kg98/Sid/STAGES/STAGES_fmriprep/data/rawdata /home/scho0011/kg98/Sid/STAGES/STAGES_fmriprep/data/derivatives/mriqc participant --nprocs 12 --participant_label ${subject}


