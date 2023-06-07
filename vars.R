par_0 <- c(N_1=0.004,N_2=-119,N_6=-23,N_7=64.5,
          M_1=0.523,M_2=-69.7,M_6=-17.6,M_7=30.8,
          H_3=-21.1,H_4=68.3,H_5=-3.67,H_6=161.5,
          g_K=36, g_Na=120, g_Leak=0.3)

with(as.list(par_0), {
  opt_par_0 <<- c(N_6=N_6,N_7=N_7,N_1=N_1,N_2=N_2,
                  M_6=M_6,M_7=M_7,M_1=M_1,M_2=M_2,
                  H_4=H_4,H_5=H_5,H_6=H_6,H_3=H_3,
                  g_Na=g_Na, g_K=g_K, g_Leak=g_Leak)
})

data_folder = paste(getwd(), "/Cleaned_Data/", sep="")
output_folder = paste(getwd(), "/Parameter_Files/", sep="")
notebook_folder = paste(getwd(), "/Notebook/", sep="")

# Example Values
AP_data_ex = read.csv(paste(data_folder, "AP_trace.csv", sep=""))
trace_data_ex = AP_data_ex$WT_avg_10ms
time_data_ex = AP_data_ex$Time_ms