#!/bin/bash
#SBATCH --job-name=ESTEBAN
#SBATCH --time=71:59:00
#SBATCH --output=test_%A_%a.out
#SBATCH --error=test_%A_%a.err
#SBATCH --array=1-12
module load R/4.0.0-foss-2020a  # NOTE

key=(ETH KEN LSO MOZ MWI NGA SWZ TZA UGA ZAF ZMB ZWE) # country key
iso=${key[$SLURM_ARRAY_TASK_ID-1]}                    # as iso3
# R --slave --vanilla --args < ../R/inference.R $iso    # run
# sensitivity analyses
R --slave --vanilla --args < ../R/inference.R $iso EPTBlo
R --slave --vanilla --args < ../R/inference.R $iso EPTBhi
