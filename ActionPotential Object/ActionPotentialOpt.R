# Optimization methods for ActionPotential
ActionPotential$methods(
  # Optimize action potential against data
  optimize = cmpfun(function(opt_params_in=all_params) {
    peak_dur = pk_dur_est
    setup_dur = 5 # Used to "calm down" beginning of AP
    Ind = (t>(tot_wait-setup_dur) & t<(tot_wait)) | 
      (t>(tot_wait+peak_dur) & t<(tot_wait+length(trace_data)))
    Ind_pk = t>(tot_wait) & t<(tot_wait+peak_dur)
    
    # Save original paramter set
    orig_param = all_params
    
    # Parition the paramter set
    opt_params = all_params[names(all_params) %in% 
                              names(opt_params_in)] #safeguard
    non_opt_params = all_params[! names(all_params) %in% 
                                  names(opt_params_in)]
    
    # Value to return when a optim condition is unwanted
    obj_reject = Inf
    
    # Penalty check for unresonably sized parameters
    check_penalty <- function(params) {
      max_up = 5 #number represented as a percent
      max_down = 0.1
      par_norm = params/orig_param[names(orig_param) %in% 
                                     names(params)]
      return(sum((par_norm < max_down) | (par_norm > max_up)))
    }
    
    # Objective function
    objective <- function(params) {
      if(check_penalty(params)) {
        return(obj_reject)
      }
      update_model(c(params,non_opt_params))
      # weighted action potential on peak
      score = 5*sum((Vs[Ind_pk]-AP_val[Ind_pk])^2) + 
        sum((Vs[Ind]-AP_val[Ind])^2)
      return(ifelse(is.na(score), obj_reject, score))
    }
    out_par = optim_standard(objective, opt_params)
    
    # Add back in paramters that were not optimized
    out_par$par = c(out_par$par,non_opt_params)
    update_subunits(out_par$par)
    return(out_par)
  }),
  
  # Multistart implementation
  nudge_optimize = cmpfun(function(opt_params=all_params, 
                                   num_fits=1, range=0.1) {
    center = ifelse(range < 1, 1, range + 0.5)
    #could change to "wander" with better fits
    time_data_loc = c(0,dt)
    
    init_par = as.data.frame(sobolSequence.points(
      length(opt_params), count=num_fits)*(2*range) + 
        (center-range))
    colnames(init_par) = names(opt_params)
    for(r in 1:nrow(init_par)) {
      init_par[r,] = init_par[r,] * opt_params
    }
    # Add in fixed parameters
    non_opt_params = all_params[! names(all_params) %in% 
                                  names(opt_params)]
    if(length(non_opt_params) > 0) {
      for(i in 1:length(non_opt_params)) {
        init_par = cbind(init_par, rep(non_opt_params[i], 
                                       nrow(init_par)))
        colnames(init_par)[ncol(init_par)] = 
          names(non_opt_params[i])
      }
    }
    
    # Optimize a number of paramter sets
    tmp_data = foreach(i = 1:num_fits, .combine = rbind, 
                       .errorhandling = 'remove') %do% {
      cur_pars = unlist(init_par[i,])
      cur_opt_pars = cur_pars[names(cur_pars) %in% 
                                names(opt_params)]
      AP_nudge = ActionPotential(cur_pars, trace_data, 
                                 time_data_loc)
      out_par = AP_nudge$optimize(cur_opt_pars)
      
      data.frame(t(sapply(out_par$par,c)), values=out_par$value, 
                 convergence=out_par$convergence)
    }
    
    # Remove "unruly" parameter sets and return the best
    tmp_data = subset(tmp_data, convergence != 10)
    data = (tmp_data[order(tmp_data$values),])
    data = data[,!(names(tmp_data)%in% 
                     c("values","convergence"))]
    update_model(unlist(data[1,]))
    return(optimize(opt_params))
  }),
  
  # Optimization routine that splits the trace into multiple 
  #   sections to optimize independently
  optimize_split = cmpfun(function() {
    peak_dur = pk_dur_est
    setup_dur = 5 # Used to "calm down" beginning of AP
    Ind_K = t>(tot_wait-setup_dur) & t<(tot_wait) | 
      t>(tot_wait+(peak_approx*0.9)) & 
      t<(tot_wait+length(trace_data))
    Ind_Na = t>(tot_wait) & t<(tot_wait+(peak_approx*1.1))
    Ind_all = t>(tot_wait-setup_dur) &
      t<(tot_wait+length(trace_data))
    
    # Gather initial parameters
    orig_param = all_params
    opt_params_K = all_params[names(all_params) %in% 
      c('N_1','N_2','N_6','N_7','g_K')]
    opt_params_Na = all_params[names(all_params) %in% 
      c('M_1','M_2','M_6','M_7','H_3','H_4','H_5','H_6','g_Na')]
    opt_params_other = all_params[names(all_params) %in% 
      c('g_Leak')]
    
    # Penalty check for unresonably sized parameters
    check_penalty <- function(params) {
      max_dev = 0.90 #number represented as a percent
      par_norm = params/orig_param[names(orig_param) %in% 
                                     names(params)]
      return(sum((par_norm < (1-max_dev)) | 
                   (par_norm > (1+max_dev))))
    }
    
    # Objective functions
    objective_all <- function(opt_params) {
      if(check_penalty(opt_params)) {
        return(obj_reject)
      }
      non_opt_params = all_params[! names(all_params) %in% 
                                    names(opt_params)]
      update_model(c(opt_params, non_opt_params))
      score = sum((Vs[Ind_all]-AP_val[Ind_all])^2)
      return(ifelse(is.na(score), obj_reject, score))
    }
    objective_K <- function(K_params) {
      if(check_penalty(K_params)) {
        return(obj_reject)
      }
      update_model(c(K_params, opt_params_Na, opt_params_other))
      score = sum((Vs[Ind_K]-AP_val[Ind_K])^2)
      return(ifelse(is.na(score), obj_reject, score))
    }
    objective_Na <- function(Na_params) {
      if(check_penalty(Na_params)) {
        return(obj_reject)
      }
      update_model(c(Na_params, opt_params_K, opt_params_other))
      score = sum((Vs[Ind_Na]-AP_val[Ind_Na])^2)
      return(ifelse(is.na(score), obj_reject, score))
    }
    
    # Rejection return value
    orig_score = sum((Vs[Ind_all]-AP_val[Ind_all])^2)
    obj_reject = Inf #abs(orig_score + 10) * 1000000
    
    # Optimize the Na and K paramters:
    # Start with a baseline with K
    out_par_K = optim_standard(objective_K, opt_params_K)
    opt_params_K = out_par_K$par
    update_model(c(opt_params_K, opt_params_Na, 
                   opt_params_other))
    
    # Cycle though Na and K
    for(i in 1:15) {
      out_par_Na = optim_standard(objective_Na, opt_params_Na)
      opt_params_Na = out_par_Na$par
      update_model(c(opt_params_K, opt_params_Na, 
                     opt_params_other))
      
      out_par_K = optim_standard(objective_K, opt_params_K)
      opt_params_K = out_par_K$par
      update_model(c(opt_params_K, opt_params_Na, 
                     opt_params_other))
    }
    
    # Finally find a fit over entire trace to show feasibility
    all_params_Na = all_params[names(all_params) %in% 
                        c('M_1','M_2','M_6','M_7','g_Na')]
    all_params_not_Na = all_params[! names(all_params) %in% 
                        names(all_params_Na)]
    out_par = nudge_optimize(all_params_not_Na, 15, 0.15)
    out_par = optimize()

    return(out_par)
  }),
  
  # Find the foot of the action potential
  find_foot = cmpfun(function() {
    # Objective function for the foot finder
    objective = function(foot_params) {
      with(as.list(foot_params), {
        if(t_0 < 0 | t_0 > peak_approx | stim_dim < 1 | 
           stim_d > peak_approx) return(obj_reject)
        configure_stim(stim_d, stim_dim, 
                ((stim_dim+1)*stim_A)/stim_d, stabil_time + t_0)
        
        # Create an array of integration values to optimize 
        #   against
        Ind_stim = t>(tot_wait) & t<(tot_wait+stim_d)
        len_stim = length(AP_val[Ind_stim])
        stim_int = rep(0, len_stim)
        for(i in 1:len_stim) {
          stim_int[i] = integrate(.self$stim_function, 
                              lower=0, upper=(i*dt))$value + RMP
        }
        
        # Compare the trace to the integral values
        return(stim_h*sum(((stim_int - 
                AP_val[Ind_stim])/AP_val[Ind_stim])^2)/len_stim)
      })
    }
    
    # Get a baseline to use to reject senseless parameters
    obj_reject = (objective(c(t_0=0, stim_d=stim_d, 
                              stim_dim=stim_dim)) + 10) * 1000
    
    # Scan start times to coax convergence
    num_fits = 10
    tmp_data = foreach(i = 1:num_fits, .combine = rbind, 
                       .errorhandling = 'remove') %do% {
      out_par = tryCatch(
        {
          optim(fn=objective, par=c(t_0=(i/num_fits), 
                              stim_d=stim_d, stim_dim=stim_dim), 
                method="Nelder-Mead", control=list(reltol=10^-8, 
                                                   maxit=1000))
        },
        error=function(e) {
          stop(paste("Error with optim while finding the foot.\n",e))
        }
      )
      
      data.frame(t(sapply(out_par$par,c)), values=out_par$value, 
                 dur=stim_d, dim=stim_dim, ht=stim_h, 
                 wait=tot_wait)
    }
    tmp_data = (tmp_data[order(tmp_data$values),])
    configure_stim(unlist(tmp_data[1,]["dur"]), 
                   unlist(tmp_data[1,]["dim"]), 
                   unlist(tmp_data[1,]["ht"]), 
                   unlist(tmp_data[1,]["wait"]))
    data = tmp_data[,!(names(tmp_data) %in% 
                    c("values", "dur", "dim", "ht", "wait"))]
    return(unlist(data[1,]))
  }),
  
  optim_standard = cmpfun(function(fcn, input) {
    opt_set = tryCatch(
      {
        optim(fn=fcn, par=input, method="Nelder-Mead", 
              control=list(parscale=input, reltol=10^-12, 
                           maxit=100000))
      },
      error=function(e) {
        stop(paste(e, AP_name))
      }
    )
    return(opt_set)
  })
)