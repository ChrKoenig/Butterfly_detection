#PBS -N MSOM_2
#PBS -l ncpus=4,mem=100GB,walltime=1000:00:00
#PBS -q maoekq
#PBS -m bea
#PBS -o ./MSOM_2_output.txt
#PBS -e ./MSOM_2_error.txt
cd $PBS_O_WORKDIR
module use /srv/ag-zurell/privatemodules/
module load jags
Rscript --vanilla ~/Butterfly_project/Butterfly_detection/2_fit_MSOM2.R
