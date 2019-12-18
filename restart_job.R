# A function to restart a job by looking at the data files and roster 
# to work out the starting points
# dir - top directory of the job
# originalJob - the original job file

restart_job <- function(dir = getwd(), params_dir  = getwd(),  
                        data_dir = getwd(), output_dir = getwd(),
                        console_files_dir = getwd(),
                        oe_files_dir = getwd(), originalJob, 
                        MbMemory = NULL, walltime = NULL,
			                  q = NULL, remove_logs = FALSE){
  
  # Get the roster
  roster <- read.csv(file = paste0(data_dir, 'parameters.csv'), 
                     stringsAsFactors = FALSE)
  
  # Get parameters from the original
  originalJob_eg <- readLines(originalJob)[1]
  if(is.null(q)) q <- gsub('\\-q\\s', '', stringr::str_extract(originalJob_eg, '\\-q\\s[[:alpha:][:punct:]]+'))
  if(is.null(MbMemory)) MbMemory <- as.numeric(gsub('mem=', '', stringr::str_extract(originalJob_eg, 'mem=[[:digit:]]+')))
  if(is.null(walltime)) walltime <- gsub('\\-W\\s', '', stringr::str_extract(originalJob_eg, '\\-W\\s[[:digit:][:punct:]]+'))
  
  data_files <- list.files(path = output_dir, pattern = '.rdata$')
  consolefiles <- list.files(path = console_files_dir, pattern = paste0('^console'),
                             full.names = TRUE)
  oe_files <- list.files(path = oe_files_dir, pattern = "\\.[oe]$",
                         full.names = TRUE)
  
  # Work out which lines in the roster have output
  expected_data <- paste0(paste(roster$species, roster$end, roster$chain, sep = '_'), '.rdata')
  roster$data <- expected_data %in% data_files
  
  # Get the bit of roster to re-run
  roster_sec <- roster[!roster$data, ]
  
  doc <- NULL
  for(sp in unique(roster_sec$species)){
    
    roster_sp <- na.omit(roster_sec[roster_sec$species == sp, ])
    
    for(chain in unique(roster_sp$chain)){
      
      roster_ch <- na.omit(roster_sp[roster_sp$chain == chain, ])
      
      job_line <- paste0('-J JAGS[', min(as.numeric(row.names(roster_ch))),
                         '-', max(as.numeric(row.names(roster_ch))), ']%1')
      submit_line <- paste('bsub', # submit
                           paste('-q', q), #queue
                           '-n 1', #N cores,
                           '-R', paste0('"rusage[mem=', MbMemory, ']"'),
                           '-M', as.integer(MbMemory*1000),
                           paste('-W', walltime), #wall time
                           job_line, #job array
                           paste0('-oo ', oe_files_dir, 'R-%J-%I.o'), #log
                           paste0('-eo ', oe_files_dir, 'R-%J-%I.e'), #log
                           '< multi_array.job') #call to R
      
      doc <- c(doc, submit_line)
      
      # remove the log files that may exist for these runs
      if(remove_logs){
        
        # delete failed console files
        con_del_index <- as.numeric(gsub('^console', '', gsub('.Rout$', '', basename(consolefiles)))) %in% min(as.numeric(row.names(roster_ch))):max(as.numeric(row.names(roster_ch)))
        unlink(x = consolefiles[con_del_index], force = TRUE)
        
        # delete failed .o and .e files
        oe_del_index <- as.numeric(gsub('\\.[oe]', '', stringr::str_extract(basename(oe_files), '[:digit:]+\\.[oe]$'))) %in% min(as.numeric(row.names(roster_ch))):max(as.numeric(row.names(roster_ch)))
        unlink(x = oe_files[oe_del_index], force = TRUE)
        
      }
    }
  }
  
  # write the sh script
  f_name <- paste0("bsub_lines_", format(Sys.time(),"%Y%m%d_%H%M%S"), '.sh')
  write.table(doc, file = file(f_name, "wb"),
              append = TRUE, col.names = FALSE,
              row.names = FALSE, quote = FALSE)
  closeAllConnections()
  
  return(f_name)
}