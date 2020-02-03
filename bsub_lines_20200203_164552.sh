bsub -q short-serial -n 1 -R "rusage[mem=20000]" -M 20000 -W 23:59 -J HRS_modelling_Job[12401-12600]%1 -oo Outputs/Log_files/R-%J-%I.o -eo Outputs/Log_files/R-%J-%I.e < multi_array.job
bsub -q short-serial -n 1 -R "rusage[mem=20000]" -M 20000 -W 23:59 -J HRS_modelling_Job[1-200]%1 -oo Outputs/Log_files/R-%J-%I.o -eo Outputs/Log_files/R-%J-%I.e < multi_array.job
bsub -q short-serial -n 1 -R "rusage[mem=20000]" -M 20000 -W 23:59 -J HRS_modelling_Job[8815-9000]%1 -oo Outputs/Log_files/R-%J-%I.o -eo Outputs/Log_files/R-%J-%I.e < multi_array.job
