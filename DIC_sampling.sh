#PBS -N DIC_sample
#PBS -l ncpus=4,mem=100GB,walltime=1000:00:00
#PBS -q maoekq
#PBS -m bea
#PBS -o ./dic_output.txt
#PBS -e ./dic_error.txt
cd $PBS_O_WORKDIR
module use /srv/ag-zurell/privatemodules/
module load jags
Rscript --vanilla ~/Butterfly_project/Butterfly_detection/3_dic_sampling.R
