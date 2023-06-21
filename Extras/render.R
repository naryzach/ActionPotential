library(stringr)
source("vars.R")

run=TRUE
cores=15
name="triplewise"
rmarkdown::render(paste("AP_", name, ".Rmd", sep=""), output_format ="html_document", 
                  output_file = paste(notebook_folder, str_to_title(name), " ", Sys.Date(), sep=""),
                  params=list(run=run,cores=cores) )