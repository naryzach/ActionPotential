# The ActionPotential(parameters, trace_data, time_data, name="Name") 
#   implementation creates an action potetnial object with given
#   data about the trace.
AP = ActionPotential(par_0, trace_data_ex, time_data_ex,"Main_Doc")

# find_foot() finds the start of the action potential and 
#   determines the shape of the foot
AP$find_foot()

# optimize(paramters) is a member function that optimizes the 
#   given paramters against the trace data supplied
AP$optimize(opt_par_0)

# display_action_potential() shows all of the plots for 
#   corresponding to the current paramters
AP$display_all_plots()