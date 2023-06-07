# Plotting methods for ActionPotential
ActionPotential$methods(
  # Method for displaying the action potential
  display_action_potential = function() {
    update_model(all_params)
    with (generate_trace(PLOT=TRUE), {
      plot(t, Vs, 
           main="Action Potential",
           xlab="Time (ms)",
           ylab="Voltage (mV)",
           type='l',
           ylim=c(-100,100))
      lines(t, Is, col="gray")
      lines(t, Istims, col="blue")
      abline(h=c(30,-80), col="black", lty=3)
      
      lines(x = t, y=AP_val, col="purple")
      legend(0, 100, 
             legend=c("AP", expression(I[m]), "Pulse", "Exp data"),
             col=c('black', 'gray', 'blue', 'purple'),
             lty=1)
    })
  },
  
  # Plot the decomposed currents for the model
  display_currents = function() {
    with(generate_trace(PLOT=TRUE), {
      y_min = min(c(Is, IKs, INas, ILeaks, Istims))
      y_max = max(c(Is, IKs, INas, ILeaks, Istims))
      plot(t, Vs, 
           main="Current Decomposition",
           xlab="Time (ms)",
           ylab="Current (mA)",
           type='l',
           ylim=c(y_min,y_max))
      lines(t, Is, col="gray")
      lines(t, IKs, col="blue")
      lines(t, INas, col="springgreen4")
      lines(t, ILeaks, col="pink")
      lines(t, Istims, col="aquamarine3")
      
      legend(0, y_max, 
             legend=c("AP", expression(I[m]), expression(I[K]), expression(I[Na]), expression(I[Leak]), "Pulse"),
             col=c('black', 'gray', 'blue', 'springgreen4', 'pink', 'aquamarine3'), 
             lty=1)
    })
  },
  
  # Plot the decomposed sub-unit conductance
  display_conductance = function() {
    with(generate_trace(PLOT=TRUE), {
      plot(t, ns, 
           main="Conductance Decomposition",
           xlab="Time (ms)",
           ylab="Particle Proability Active",
           type='l',
           col="blue",
           ylim=c(0,1))
      lines(t, ms, col="springgreen4")
      lines(t, hs, col="red")
      
      legend(0, 1, 
             legend=c("n", "m", "h"),
             col=c('blue', 'springgreen4', 'red'), 
             lty=1)
      abline(v=tot_wait, col="black", lty=3)
    })
  },
  
  # Test alpha and beta
  display_alpha_beta = function() {
    V = seq(-150,100,1)
    plot_max = 2
    plot(V, n@alpha(V), 
         main="Alpha: open; Beta: closed",
         xlab="Voltage (mV)",
         ylab="Rate constant (1/mS)", 
         type='l', 
         col="blue",
         lty=1,
         ylim=c(0,plot_max))
    lines(V, n@beta(V), col="blue", lty=2)
    lines(V, m@alpha(V), col="springgreen4", lty=1)
    lines(V, m@beta(V), col="springgreen4", lty=2)
    lines(V, h@alpha(V), col="red", lty=1)
    lines(V, h@beta(V), col="red", lty=2)
    legend(50, plot_max, 
           legend=c(expression(alpha[n]), expression(beta[n]), 
                    expression(alpha[m]), expression(beta[m]), 
                    expression(alpha[h]), expression(beta[h])),
           col=c('blue', 'blue', 'springgreen4', "springgreen4", "red", "red"), 
           lty=c(1,2,1,2,1,2))
  },
  
  # Testing Subunit functions
  display_subunits = function() {
    #V_init = -80
    ns = c(infty(n,V_init))
    ms = c(infty(m,V_init))
    hs = c(infty(h,V_init))
    for (timestep in tail(t,-1)) {
      n_t = tail(ns, 1)
      m_t = tail(ms, 1)
      h_t = tail(hs, 1)
      V = 0
      if (timestep>(time/2))
        V = -80
      ns = append(ns, n_t + dt*dn(n, n_t, V))
      ms = append(ms, m_t + dt*dm(m, m_t, V))
      hs = append(hs, h_t + dt*dh(h, h_t, V))
    }
    plot(t, ns, 
         main="Subunit Conductance",
         xlab="Time (ms)",
         ylab="Particle Proability Active", 
         ylim=c(0,1),
         type='l', 
         col="blue",
         lty=1)
    lines(t, ms, col="springgreen4", lty=1)
    lines(t, hs, col="red", lty=1)
    lines(t, ns^4, col="blue", lty=2)
    lines(t, ms^3*hs, col="springgreen4", lty=2)
    legend(time-5, 1, 
           legend=c('n', 'm', 'h', "K", "Na"),
           col=c('blue', 'springgreen4', 'red', "blue", "springgreen4"), 
           lty=c(1,1,1,2,2))
  },
  
  # Generate plots for Taus
  display_taus = function() {
    V = seq(-100,100,1)
    plot_max = 10
    plot(V, tau(n,V), 
         main="Time Constants", 
         xlab="Voltage (mV)",
         ylab="Tau (ms)", 
         type='l',
         ylim=c(0,plot_max), 
         col="blue")
    lines(V, tau(m,V), col="springgreen4")
    lines(V, tau(h,V), col="red")
    legend(50, plot_max, 
           legend=c('n', 'm', 'h'),
           col=c('blue', 'springgreen4', 'red'), 
           lty=1)
  },
  
  # Generate plots for Infinities
  display_infinities = function() {
    V = seq(-150,100,1)
    plot(V, infty(n,V), 
         main="Steady-Sate Solutions", 
         xlab="Voltage (mV)",
         ylab="Particle Probabiliy Active", 
         type='l',
         col="blue",
         ylim=c(0,1))
    lines(V, infty(m,V), col="springgreen4")
    lines(V, infty(h,V), col="red")
    legend(50, 0.9, 
           legend=c('n', 'm', 'h'),
           col=c('blue', 'springgreen4', 'red'), 
           lty=1)
  },
  
  display_all_plots = function() {
    display_action_potential()
    display_currents()
    display_conductance()
    display_subunits()
    display_alpha_beta()
    display_taus()
    display_infinities()
  },
  
  # Plot the AP simulation across a range of values for each parameter
  plot_param_range = function() {
    for(i in 1:(length(all_params))) {
      update_model(all_params)
      par_tmp = all_params[i]
      max_seq = 2*par_tmp
      plot(t, Vs, 
        main=paste(names(par_tmp), ":", par_tmp),
        xlab="Time (ms)",
        ylab="Voltage (mV)", 
        type='l',
        ylim=c(-100,100),
        col="yellow")
      #lines(t, Istims, col="green")
      
      legend(20, 75, legend=c(0, max_seq), 
             fill = c(rgb(max_seq,0,0, max=max_seq),rgb(0,0,max_seq, max=max_seq)))
      
      par = seq(0,max_seq,length.out=100)
      for(val in par) {
        all_params[i] <<- val
        update_model(all_params)
        lines(t, Vs, col=rgb(max_seq-val,0,val, max=max_seq))
      }
      print(par_tmp)
      all_params[i] <<- par_tmp
    }
  }
)