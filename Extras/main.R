source("ActionPotential.R")
source("vars.R")

main <- function() {
  args = commandArgs(trailingOnly=TRUE)
  cores = 35
  if (length(args) > 0) {
    if (args[1] == "--cores") {
      cores = strtoi(args[2]) 
    }
  }
  #test()
  AP_ex = ActionPotential(par_0, trace_data_ex, time_data_ex)
  render(write=F, trace=F, group=F, group_trace=T, stat=F, run=T, cores=cores)
}

render <- function(write=FALSE, trace=FALSE, group=FALSE, group_trace=FALSE, stat=FALSE, run=FALSE, cores=1) {
  if (write) {
    rmarkdown::render("AP_writeup.Rmd", output_format ="html_document", 
                      output_file = paste(notebook_folder, "Writeup ", Sys.Date(), sep=""),
                      params=list(run=run) )
  }
  
  if (trace) {
    rmarkdown::render("AP_read_traces.Rmd", output_format ="html_document", 
                      output_file = paste(notebook_folder, "Trace ", Sys.Date(), sep=""),
                      params=list(run=run,cores=cores) )
  }
  
  if (group) {
    rmarkdown::render("AP_group_analysis.Rmd", output_format ="html_document", 
                      output_file = paste(notebook_folder, "Group ", Sys.Date(), sep=""),
                      params=list(run=run, cores=cores) )
  }
  
  if (group_trace) {
    rmarkdown::render("AP_group_trace.Rmd", output_format ="html_document", 
                      output_file = paste(notebook_folder, "Group_Data ", Sys.Date(), sep=""),
                      params=list(run=run, cores=cores) )
  }
  
  if (stat) {
    rmarkdown::render("AP_stats.Rmd", output_format ="html_document", 
                      output_file = paste(notebook_folder, "Stat ", Sys.Date(), sep=""),
                      params=list(run=run, cores=cores) )
  }
}

test <- function() {
  source("ActionPotential.R")
  source("vars.R")
  #file = "Atratus_WT.csv"
  #trace_file = "trace_par"
  group_file = "group_par"
  all_tables = vector("list", 25)
  foreach(s = 1:25, .errorhandling = 'remove') %do% {
    all_tables[[s]] = readRDS(file = paste(output_folder, "Group/", group_file, s, "data_table.Rds", sep=""))
  }
  
  table_means = foreach (d = 1:length(all_tables), .combine = rbind, .errorhandling = 'remove') %do% {
    dt = all_tables[[d]]
    group_data1 = dt[dt$group==1,]
    group_data2 = dt[dt$group==2,]
    
    mean_diff = colMeans(group_data2) - colMeans(group_data1)
  }
  print("Generated table means")
  print(colnames(table_means))
  
  for (col in 1:ncol(table_means)) {
    hist(table_means[,col], main=paste("Histogram for", colnames(table_means)[col]), xlab="Parameter Average")
  }
}

main()