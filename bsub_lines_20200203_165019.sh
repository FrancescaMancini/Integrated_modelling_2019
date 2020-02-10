bsub -q short-serial -n 1 -R "rusage[mem=20000]" -M 20000 -W 23:59 -J Integrated_modelling_Job[12401-12600]%1 -oo Outputs/Log_files_int/R-%J-%I.o -eo Outputs/Log_files_int/R-%J-%I.e < multi_array_int.job
bsub -q short-serial -n 1 -R "rusage[mem=20000]" -M 20000 -W 23:59 -J Integrated_modelling_Job[1-200]%1 -oo Outputs/Log_files_int/R-%J-%I.o -eo Outputs/Log_files_int/R-%J-%I.e < multi_array_int.job
bsub -q short-serial -n 1 -R "rusage[mem=20000]" -M 20000 -W 23:59 -J Integrated_modelling_Job[8829-9000]%1 -oo Outputs/Log_files_int/R-%J-%I.o -eo Outputs/Log_files_int/R-%J-%I.e < multi_array_int.job
bsub -q short-serial -n 1 -R "rusage[mem=20000]" -M 20000 -W 23:59 -J Integrated_modelling_Job[10477-10600]%1 -oo Outputs/Log_files_int/R-%J-%I.o -eo Outputs/Log_files_int/R-%J-%I.e < multi_array_int.job
bsub -q short-serial -n 1 -R "rusage[mem=20000]" -M 20000 -W 23:59 -J Integrated_modelling_Job[7756-7800]%1 -oo Outputs/Log_files_int/R-%J-%I.o -eo Outputs/Log_files_int/R-%J-%I.e < multi_array_int.job