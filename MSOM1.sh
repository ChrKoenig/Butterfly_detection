#PBS -N MSOM_1_par
#PBS -l ncpus=4,mem=100GB,walltime=1000:00:00
#PBS -q maoekq
#PBS -m bea
#PBS -o ./Data/models_fit/MSOM/MSOM_1_rjags_output.txt
#PBS -e ./Data/models_fit/MSOM/MSOM_1_rjags_error.txt
cd $PBS_O_WORKDIR
module use /srv/ag-zurell/privatemodules/
module load jags
Rscript --vanilla ~/Butterfly_project/Butterfly_detection/2_fit_MSOM1.R
