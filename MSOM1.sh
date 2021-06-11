#PBS -N MSOM_1
#PBS -l ncpus=4,mem=100GB,walltime=1000:00:00
#PBS -q maoekq
#PBS -m bea
#PBS -o ./MSOM_1_output.txt
#PBS -e ./MSOM_1_error.txt
cd $PBS_O_WORKDIR
module use /srv/ag-zurell/privatemodules/
module load jags
Rscript --vanilla ~/Butterfly_project/Butterfly_detection/2_fit_MSOM1.R
