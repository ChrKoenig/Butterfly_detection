#PBS -N MSOM_3
#PBS -l ncpus=4,mem=60GB,walltime=1000:00:00
#PBS -q maoekq
#PBS -m bea
#PBS -o ./Data/models_fit/MSOM/MSOM_3_output.txt
#PBS -e ./Data/models_fit/MSOM/MSOM_3_error.txt
cd $PBS_O_WORKDIR
module use /srv/ag-zurell/privatemodules/
module load jags
Rscript --vanilla ~/Butterfly_project/Butterfly_detection/2_fit_MSOM3.R
