# Libraries needed for ActionPotential
library(methods)
library(deSolve)
library(compiler)
source("Subunit.R")

# Action Potential object declaration
ActionPotential <- setRefClass("ActionPotential",
  fields = list(
    # Initial Data
    all_params="vector",trace_data="vector",dt="numeric",
    AP_name="character",
    
    # Sub-unit Declarations
    n="Subunit", m="Subunit",h="Subunit",
    K="Channel", Na="Channel", Leak="Channel",
    
    # Stimulus Values
    stim_d="numeric", stim_h="numeric", stim_dim="numeric",
    stim_A="numeric",
    
    # Action Potential Values
    time="numeric", stabil_time="numeric", tot_wait="numeric",
    t="vector", peak_approx="numeric", pk_dur_est="numeric",
    AP_max="numeric", AP_min="numeric", AP_height="numeric",
    RMP="numeric", V_init="numeric", Vs="vector", AP_val="vector"
  )
)

# Core methods for ActionPotential
ActionPotential$methods(
  # Set up action potential class
  initialize = function(in_all_params, in_trace_data, 
                        in_time_data, name="Nameless") {
    # Set ActionPotential name
    AP_name <<- name
    
    # Default values
    time <<- 30
    stabil_time <<- 10
    
    # Grab data from input
    all_params <<- in_all_params
    dt <<- round((in_time_data[2]-in_time_data[1]),3)
    
    # Setup Subunits
    n <<- new("Subunit")
    m <<- new("Subunit")
    h <<- new("Subunit")
    K <<- K_orig
    Na <<- Na_orig
    Leak <<- Leak_orig
    update_subunits(all_params)
    
    # Stimulus variables
    stim_d <<- 0.64
    stim_h <<- 54.874
    stim_dim <<- 2
    stim_A <<- integrate(.self$stim_function, lower=0, 
                         upper=stim_d)$value
    tot_wait <<- stabil_time
    
    # Approximate information about experimental traces
    peak_approx <<- 2 # Approximate time of peak after stimulus
    pk_dur_est <<- 5 # Approximate duration of the peak
    
    # Prepare action potential
    trace_data <<- in_trace_data
    AP_max <<- max(trace_data)
    AP_min <<- min(trace_data)
    AP_height <<- AP_max-AP_min
    RMP <<- min(trace_data[1:(peak_approx/dt)]) + 
      (max(trace_data[((pk_dur_est+1)/dt):((pk_dur_est+1.2)/dt)])
      -min(trace_data[((pk_dur_est+1)/dt):((pk_dur_est+1.2)/dt)]) 
      )/2
    V_init <<- RMP
    l = length(trace_data)
    AP_val <<- c(rep(RMP,(stabil_time)/dt+1), 
                 trace_data, rep(RMP,(time-(stabil_time))/dt-l))
    
    # Set up timing
    t <<- seq(0,time,dt) # time in ms
    Vs <<- rep(V_init, length(t))
    
    # Preemptively calulate Vs
    generate_trace()
  },
  
  # Calculate the trace values
  generate_trace = cmpfun(function(PLOT=FALSE) {
    # Initialize empty model data
    Vs_loc = rep(V_init, length(t))
    
    # Include tracking if plotting
    if (PLOT) {
      Is = rep(0, length(t))
      Istims = rep(0, length(t))
      INas = rep(0, length(t))
      IKs = rep(0, length(t))
      ILeaks = rep(0, length(t))
      Istims = rep(0, length(t))
      ns = rep(infty(n, V_init), length(t))
      ms = rep(infty(m, V_init), length(t))
      hs = rep(infty(h, V_init), length(t))
    }
    
    # Keep track of iterations
    i = 2
    
    # Initialize sub-unit conditions
    n_t = infty(n, V_init)
    m_t = infty(m, V_init)
    h_t = infty(h, V_init)
    V_m = V_init
    for (timestep in tail(t,-1)) {
      # Inject pulse after stabilization delay
      if (timestep > tot_wait && timestep < (tot_wait + stim_d)){
        I_stim = stim_function(timestep-tot_wait)
      } else {
        I_stim = 0
      }
      
      # Compute time step
      n_t = n_t + dt*dn(n, n_t, V_m)
      m_t = m_t + dt*dm(m, m_t, V_m)
      h_t = h_t + dt*dh(h, h_t, V_m)
      I_K = K@g*n_t^4*(V_m-K@E)
      I_Na = Na@g*m_t^3*h_t*(V_m-Na@E)
      I_Leak = Leak@g*(V_m-RMP)
      I_m = I_K + I_Na + I_Leak
      V_m = V_m + dt*dV(I_m, I_stim)
      
      # Log time step
      Vs_loc[i] = V_m
      if (PLOT) {
        Is[i] = I_m
        Istims[i] = I_stim
        INas[i] = I_Na
        IKs[i] = I_K
        ILeaks[i] = I_Leak
        Istims[i] = I_stim
        ns[i] = n_t
        ms[i] = m_t
        hs[i] = h_t
      }
      i = i+1
    }
    
    # Return plotting information if requested
    if (PLOT) {
      return(list(Vs=Vs_loc, Is=Is, INas=INas, IKs=IKs, 
                  ILeaks=ILeaks, Istims=Istims, 
                  ns=ns, ms=ms, hs=hs))
    }
    
    # Update Vs
    Vs <<- Vs_loc
  }),
  
  # Define derivative functions
  dn = cmpfun(function(n, n_prob, V) 
    return(n@alpha(V)*(1-n_prob)-n@beta(V)*n_prob)),
  dm = cmpfun(function(m, m_prob, V) 
    return(m@alpha(V)*(1-m_prob)-m@beta(V)*m_prob)),
  dh = cmpfun(function(h, h_prob, V) 
    return(h@alpha(V)*(1-h_prob)-h@beta(V)*h_prob)),
  dV = cmpfun(function(I_m, I_stim) 
    return(I_stim-I_m)),
  
  # Function to generate trace using ode method
  gen_trace_ode = cmpfun(function(params=all_params) {
    trace_soln = ode(func=.self$AP_ode, parms=params,
                     y=c(V=V_init, n=infty(n, V_init), 
                         m=infty(m, V_init), h=infty(h, V_init)),
                     times=t, method="ode23")
    Vs <<- unlist(trace_soln[,2])
  }),
  
  # Function of time, t, with in_vect = c(V,n,m,h), and params
  AP_ode = cmpfun(function(t, in_vect, params) {
    with(as.list(c(in_vect,params)), {
      n_alpha = ifelse(V<N_2, 0, (N_1*V-N_1*N_2))
      n_beta = exp((V+N_7)/N_6)
      m_alpha = ifelse(V<M_2, 0, (M_1*V-M_1*M_2))
      m_beta = exp((V+M_7)/M_6)
      h_alpha = exp((V+H_6)/H_3)
      h_beta = 1/(1+exp((V+H_4)/H_5))
      
      dn = n_alpha*(1-n)-n_beta*n
      dm = m_alpha*(1-m)-m_beta*m
      dh = h_alpha*(1-h)-h_beta*h
      
      I_K = K@g*n^4*(V-K@E)
      I_Na = Na@g*m^3*h*(V-Na@E)
      I_Leak = Leak@g*(V-RMP)
      I_m = I_K + I_Na + I_Leak
      I_stim = ifelse( t > tot_wait && t < (tot_wait + stim_d), 
                       stim_function(t-tot_wait), 0) 
      dV = I_stim-I_m
      
      list(c(dV, dn, dm, dh))
    })
  }),
  
  # Function for stimulus shape
  stim_function = cmpfun(function(x) {
    return(stim_h*(x/stim_d)^stim_dim)
  }),
  
  # Update the alpha/beta functions with new parameters
  update_subunits = cmpfun(function(params) {
    all_params <<- params
    with(as.list(all_params), {
      n@alpha <<- function(V) ifelse(V<N_2, 0, (N_1*V-N_1*N_2)) #N_1*(V+N_2)/(exp((V+N_2)/N_3)-1)
      n@beta <<- function(V) exp((V+N_7)/N_6)
      m@alpha <<- function(V) ifelse(V<M_2, 0, (M_1*V-M_1*M_2)) #M_1*(V+M_2)/(exp((V+M_2)/M_3)-1)
      m@beta <<- function(V) exp((V+M_7)/M_6)
      h@alpha <<- function(V) exp((V+H_6)/H_3)
      h@beta <<- function(V) 1/(1+exp((V+H_4)/H_5))
      K@g <<- g_K; Na@g <<- g_Na; Leak@g <<- g_Leak
    })
  }),
  
  # Update the model with new parameters
  update_model = cmpfun(function(params) {
    update_subunits(params)
    generate_trace() # can use gen_trace_ode()
  }),
  
  # Update the values for the stimulus
  configure_stim = function(stim_d, stim_dim, stim_h, tot_wait) {
    stim_d <<- stim_d
    stim_dim <<- stim_dim
    stim_h <<- stim_h
    tot_wait <<- tot_wait
  },
  
  # Return the trace_data
  get_trace_data = function() {
    return(Vs)
  },
  
  trace_data_as_ex_data = function() {
    generate_trace()
    return(Vs[(stabil_time/dt):
                (stabil_time/dt+length(trace_data))])
  }
)

# Source Supporting method files of ActionPotential
source("ActionPotentialShow.R")
source("ActionPotentialOpt.R")
source("ActionPotentialStat.R")