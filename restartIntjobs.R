source("restart_failures.R")
source("job_status.R")

r_path <- restart_failures(dir = getwd(), originalJob = "bsub_lines_int.sh", remove_logs = TRUE, datafiles_path = "Outputs/Integrated/", 
                           consolefiles_path = "Outputs/Console_int", logfiles_path = "Outputs/Log_files_int", 
                           roster_file = "Data/parameters.csv", job_file = "multi_array_int.job")
