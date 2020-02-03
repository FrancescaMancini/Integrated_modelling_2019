source("restart_failures.R")
source("job_status.R")

r_path <- restart_failures(dir = getwd(), originalJob = "bsub_lines.sh", remove_logs = TRUE, datafiles_path = "Outputs/HRS/", 
                           consolefiles_path = "Outputs/Console", logfiles_path = "Outputs/Log_files", 
                           roster_file = "Data/parameters.csv", job_file = "multi_array.job", first.it = 200000, target.it = 400000)
                           