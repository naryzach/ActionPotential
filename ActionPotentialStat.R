library(SobolSequence)
library(doParallel)
library(foreach)

# Class for ActionPotential statistics
ActionPotentialStat <- setRefClass("ActionPotentialStat", contains="ActionPotential",
  fields = list(
    opt_params="vector",
    # Files
    stat_file="character",
    og_stat_file="character",
    stat_ref_data="character",
    opt_ref_data="character"
  )
)

# Statistics methods for ActionPotential
ActionPotentialStat$methods(
  initialize = function(...) {
    stat_file <<- "stat_par"
    og_stat_file <<- "stat_par_og"
    stat_ref_data <<- "stat_AP_ref"
    opt_ref_data <<- "opt_par_ref"
    callSuper(...)
  },
  
  # Perform the statistical analysis
  run_statistical_analysis = cmpfun(function(cores=(detectCores() - 1), opt_params=all_params) {
    num_fits = 5000
    range = 2.5
    center = ifelse(range < 1, 1, range + 0.1)
    
    # Generate trace data to use for comparison
    generate_trace()
    trace_data_gen = Vs[(stabil_time/dt):(stabil_time/dt+length(trace_data))]
    time_data_gen = c(0, dt)
    
    # Save reference data to refer to later
    AP_ref = ActionPotential(all_params, trace_data_gen, time_data_gen)
    saveRDS(AP_ref, file = paste(output_folder, stat_ref_data, ".Rds", sep=""))
    saveRDS(opt_params, file = paste(output_folder, opt_ref_data, ".Rds", sep=""))
    
    print(paste("Fits: ", num_fits, ". Range: ", range, ". Center: ", center, ". Cores: ", cores, sep = ""))
    
    # Generate quasi-random parameter sets
    init_par = as.data.frame(sobolSequence.points(length(opt_params), count=num_fits)*(2*range) + (center-range))
    colnames(init_par) = names(opt_params)
    for(r in 1:nrow(init_par)) {
      init_par[r,] = init_par[r,] * opt_params
    }
    # Add in fixed parameters
    non_opt_params = all_params[! names(all_params) %in% names(opt_params)]
    if(length(non_opt_params) > 0) {
      for(i in 1:length(non_opt_params)) {
        init_par = cbind(init_par, rep(non_opt_params[i], nrow(init_par)))
        colnames(init_par)[ncol(init_par)] = names(non_opt_params[i])
      }
    }
    saveRDS(init_par, file = paste(output_folder, og_stat_file, ".Rds", sep=""))
    
    if (file.exists(paste(stat_file, ".Rds", sep=""))) {
      file.remove(paste(stat_file, ".Rds", sep=""))
    }
    
    # Run computation for identifiability analysis
    clust = makeCluster(cores)
    clusterExport(clust, c("init_par", "opt_params", "trace_data_gen", "time_data_gen"), envir=environment())
    registerDoParallel(clust)
    data_one <- foreach(i = 1:num_fits, .combine = rbind, .errorhandling = 'remove') %dopar% { 
      source("ActionPotential.R")
      cur_pars = unlist(init_par[i,])
      cur_opt_pars = cur_pars[names(cur_pars) %in% names(opt_params)]
      AP_stat = ActionPotential(cur_pars, trace_data_gen, time_data_gen)
      out_par = AP_stat$optimize(cur_opt_pars)
      
      data.frame(t(sapply(out_par$par,c)), values=out_par$value, convergence=out_par$convergence)
    }
    stopCluster(clust)
    
    data_one = (data_one[order(data_one$values),])
    min_val = data_one[,"values"][1]
    new_max_val = min_val*1.2
    max_range = length(which(data_one[,"values"] <= new_max_val)) # take all values 20% away from the min value
    data_one = data_one[1:ifelse(max_range<10,10,max_range),]
    data_one = data_one[,!(names(data_one) %in% c("values","convergence"))]
    mult = 10
    for (i in 1:nrow(data_one)) {
      for (j in 1:mult) {
        data_one = rbind(data_one, data_one[i,])
      }
    }
    num_fits_two = nrow(data_one)
    
    # Run computation for identifiability analysis a second time
    clust = makeCluster(cores)
    clusterExport(clust, c("data_one", "opt_params", "trace_data_gen", "time_data_gen"), envir=environment())
    registerDoParallel(clust)
    data <- foreach(i = 1:num_fits_two, .combine = rbind, .errorhandling = 'remove') %dopar% { 
      source("ActionPotential.R")
      cur_pars_tmp = unlist(data_one[i,])
      cur_non_opt_pars = cur_pars_tmp[! names(cur_pars_tmp) %in% names(opt_params)]
      cur_opt_pars = cur_pars_tmp[names(cur_pars_tmp) %in% names(opt_params)]
      cur_opt_pars = cur_opt_pars * runif(length(cur_opt_pars), min=0.95, max=1.05)
      cur_pars = c(cur_opt_pars, cur_non_opt_pars)
      AP_stat = ActionPotential(cur_pars, trace_data_gen, time_data_gen)
      out_par = AP_stat$optimize(cur_opt_pars)
      
      data.frame(t(sapply(out_par$par,c)), values=out_par$value, convergence=out_par$convergence)
    }
    stopCluster(clust)
    
    # Save output data
    saveRDS(data, file = paste(output_folder, stat_file, ".Rds", sep=""))
    saveRDS(data, file = paste(output_folder, stat_file, num_fits, "-", range, ".Rds", sep=""))
  }),
  
  # Display the results of the statistical analysis
  show_statistical_analysis = cmpfun(function() {
    tmp_data <- subset(readRDS(file = paste(output_folder, stat_file, ".Rds", sep="")), convergence != 10)
    data <- (tmp_data[order(tmp_data$values),])#[1:(nrow(tmp_data)*0.2),] # Take top 20%
    
    # 3D plot setup; need library(rgl)
    #plot3d(data$M_4,data$M_5,data$M_6,col=colorRampPalette(c("green","blue","red"))(nrow(data)))
    
    # Read reference data used to generate the statistics
    AP_ref = readRDS(file = paste(output_folder, stat_ref_data, ".Rds", sep=""))
    print("Reference Data")
    print(AP_ref$all_params)
    print("Optimized paramters")
    print(readRDS(file = paste(output_folder, opt_ref_data, ".Rds", sep="")))
    AP_ref$display_action_potential()
    
    # Show the plots generated by optimized stats data
    display_stat_APs(data, AP_ref)
    
    # Calcualte coefficients of variation
    print(paste("n=", nrow(data), " "))
    for (i in 1:ncol(data)) {
      print( paste("CV for", colnames(data)[i], sd(data[,i])/mean(data[,i])) )
    }
    
    # Plot pairwise comparisons of parameters
    par(mfrow=c(1,1))
    for(i in 1:ncol(data)){
      for(j in 1:ncol(data)) {
        if(i<j){
          plot(data[,i],data[,j],
               main=paste(colnames(data)[i], "vs", colnames(data)[j]),
               xlab=colnames(data)[i],
               ylab=colnames(data)[j],
               col=colorRampPalette(c("blue","red"))(nrow(data)) )
          print(summary(lm(data[,i]~data[,j])))
        }
      }
    }
    
    pairs(data)
    return(data)
  }),
  
  # Show action potentials generated for each optimized parameter set
  display_stat_APs = function(data, AP_ref) {
    AP_ref$update_model(unlist(data[1,]))
    with(AP_ref$generate_trace(PLOT=TRUE), {
      plot(AP_ref$t, AP_ref$Vs, 
           main="All Action Potentials",
           ylab="Voltage (mV)", 
           xlab="Time (mS)",
           type='l',
           ylim=c(-100,100))
      abline(h=c(30,-80), col="black", lty=3)
    })
    for(i in 2:nrow(data)) {
      AP_ref$update_model(unlist(data[i,]))
      with(AP_ref$generate_trace(PLOT=TRUE), {
        lines(AP_ref$t, AP_ref$Vs, col=rgb(nrow(data)-i,0,i, max=nrow(data)))
      })
    }
  }
)