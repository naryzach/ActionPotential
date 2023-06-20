from numpy import pi, sin, exp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

def signal(Vs):
	print(t/dt)
	return Vs[int(t/dt)]

# Action Potential Information
N_1 = 0.004; N_2 = -119; N_6 = -23; N_7 = 64.5
M_1 = 0.523; M_2 = -69.7; M_6 = -17.6; M_7 = 30.8
H_3 = -21.1; H_4 = 68.3; H_5 = -3.67; H_6 = 161.5;
g_K = 36; g_Na = 120; g_Leak = 0.3
E_K = -83.6; E_Na = 69.3; E_Leak = -77

stim_d = 0.64; stim_h = 55; stim_dim = 2

V_init = RMP = -71
tot_wait = 10
dt = 0.005
time = 30.0

t = np.arange(0.0, time, dt)

def alpha(V, P_1, P_2):
	if (V<P_2):
		return 0
	else:
		return P_1*V-P_1*P_2
def beta(V, P_6, P_7):
	return exp((V+P_7)/P_6)

def alpha_h(V, H_3, H_6):
	return exp((V+H_6)/H_3)
def beta_h(V, H_4, H_5):
	return 1/(1+exp((V+H_4)/H_5))

def AP(N_1, N_2, N_6, N_7,
		M_1, M_2, M_6, M_7,
		H_3, H_4, H_5, H_6,
		g_K, g_Na, g_Leak,
		stim_d, stim_h, stim_dim, V_init,
		mode):
	# Initialize empty model data
	data = [0 for i in range(len(t))]

	# Initialize sub-unit conditions
	V = RMP = V_init
	n = alpha(V, N_1, N_2) / (alpha(V, N_1, N_2) + beta(V, N_6, N_7))
	m = alpha(V, M_1, M_2) / (alpha(V, M_1, M_2) + beta(V, M_6, M_7))
	h = alpha_h(V, H_3, H_6) / (alpha_h(V, H_3, H_6) + beta_h(V, H_4, H_5))

	# Add in first data_point
	if mode == 'V':
		data[0] = V_init
	elif mode == 'I' or mode == 'N' or mode == 'K' or mode == 'L':
		data[0] = 0
	elif mode == 'n':
		data[0] = n
	elif mode == 'm':
		data[0] = m
	elif mode == 'h':
		data[0] = h
	for i, timestep in enumerate(t[1:]):
		# Inject pulse after stabilization delay
		if (timestep > tot_wait) and (timestep < (tot_wait + stim_d)):
			I_stim = stim_h*(( timestep-tot_wait )/stim_d)**stim_dim
		else:
			I_stim = 0

		# Compute alphas/betas
		n_alpha = alpha(V, N_1, N_2)
		n_beta = beta(V, N_6, N_7)

		m_alpha = alpha(V, M_1, M_2)
		m_beta = beta(V, M_6, M_7)

		h_alpha = alpha_h(V, H_3, H_6)
		h_beta = beta_h(V, H_4, H_5)

		# Compute timestep
		n = n + dt * (n_alpha*(1-n)-n_beta*n)
		m = m + dt * (m_alpha*(1-m)-m_beta*m)
		h = h + dt * (h_alpha*(1-h)-h_beta*h)
		I_K = g_K * (n**4) * (V-E_K)
		I_Na = g_Na * (m**3) * h * (V-E_Na)
		I_Leak = g_Leak * (V-RMP)
		I_m = I_K + I_Na + I_Leak
		V = V + dt * (I_stim-I_m)

		# Log time step
		if mode == 'V':
			data[i+1] = V
		elif mode == 'I':
			data[i+1] = I_m
		elif mode == 'N':
			data[i+1] = I_Na
		elif mode == 'K':
			data[i+1] = I_K
		elif mode == 'L':
			data[i+1] = I_Leak
		elif mode == 'n':
			data[i+1] = n
		elif mode == 'm':
			data[i+1] = m
		elif mode == 'h':
			data[i+1] = h
	return data

axis_color = 'lightgoldenrodyellow'

fig = plt.figure()
ax = fig.add_subplot(111)

# Adjust the subplots region to leave some space for the sliders and buttons
fig.subplots_adjust(left=0.25, bottom=0.6)

# Draw the initial plot
# The 'line' variable is used for modifying the line later
Vs_orig = AP(N_1, N_2, N_6, N_7,
		M_1, M_2, M_6, M_7,
		H_3, H_4, H_5, H_6,
		g_K, g_Na, g_Leak,
		stim_d, stim_h, stim_dim, V_init, 'V')
[line] = ax.plot(t, Vs_orig, linewidth=2, color='red')
ax.set_xlim([0, time])
ax.set_ylim([-100, 100])
line.set_color("purple")

# Add two sliders for tweaking the parameters

# Define an axes area and draw all sliders
stim_d_slider_ax  = fig.add_axes([0.25, 0.5, 0.65, 0.03], facecolor=axis_color)
stim_d_slider = Slider(stim_d_slider_ax, 'stim dur', 0.1, 10.0, valinit=stim_d)

stim_h_slider_ax = fig.add_axes([0.25, 0.475, 0.65, 0.03], facecolor=axis_color)
stim_h_slider = Slider(stim_h_slider_ax, 'stim ht', 0.1, 100, valinit=stim_h)

stim_dim_slider_ax = fig.add_axes([0.25, 0.45, 0.65, 0.03], facecolor=axis_color)
stim_dim_slider = Slider(stim_dim_slider_ax, 'stim dim', 1.0, 10.0, valinit=stim_dim)

V_init_slider_ax = fig.add_axes([0.25, 0.425, 0.65, 0.03], facecolor=axis_color)
V_init_slider = Slider(V_init_slider_ax, 'V_init/RMP', -120, -0.1, valinit=V_init)

# Potassium
N_1_slider_ax  = fig.add_axes([0.25, 0.4, 0.65, 0.03], facecolor=axis_color)
N_1_slider = Slider(N_1_slider_ax, 'N_1', 0.0001, 0.01, valinit=N_1)

N_2_slider_ax = fig.add_axes([0.25, 0.375, 0.65, 0.03], facecolor=axis_color)
N_2_slider = Slider(N_2_slider_ax, 'N_2', -200, -0.1, valinit=N_2)

N_6_slider_ax = fig.add_axes([0.25, 0.35, 0.65, 0.03], facecolor=axis_color)
N_6_slider = Slider(N_6_slider_ax, 'N_6', -50, -0.1, valinit=N_6)

N_7_slider_ax = fig.add_axes([0.25, 0.325, 0.65, 0.03], facecolor=axis_color)
N_7_slider = Slider(N_7_slider_ax, 'N_7', 0.1, 150, valinit=N_7)

# Sodium activation
M_1_slider_ax  = fig.add_axes([0.25, 0.3, 0.65, 0.03], facecolor=axis_color)
M_1_slider = Slider(M_1_slider_ax, 'M_1', 0.01, 1.0, valinit=M_1)

M_2_slider_ax = fig.add_axes([0.25, 0.275, 0.65, 0.03], facecolor=axis_color)
M_2_slider = Slider(M_2_slider_ax, 'M_2', -150, -0.1, valinit=M_2)

M_6_slider_ax = fig.add_axes([0.25, 0.25, 0.65, 0.03], facecolor=axis_color)
M_6_slider = Slider(M_6_slider_ax, 'M_6', -30, -0.1, valinit=M_6)

M_7_slider_ax = fig.add_axes([0.25, 0.225, 0.65, 0.03], facecolor=axis_color)
M_7_slider = Slider(M_7_slider_ax, 'M_7', 0.1, 100, valinit=M_7)

# Sodium deactivation
H_3_slider_ax  = fig.add_axes([0.25, 0.2, 0.65, 0.03], facecolor=axis_color)
H_3_slider = Slider(H_3_slider_ax, 'H_3', -50, -0.1, valinit=H_3)

H_4_slider_ax = fig.add_axes([0.25, 0.175, 0.65, 0.03], facecolor=axis_color)
H_4_slider = Slider(H_4_slider_ax, 'H_4', 0.1, 150, valinit=H_4)

H_5_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axis_color)
H_5_slider = Slider(H_5_slider_ax, 'H_5', -10, -0.01, valinit=H_5)

H_6_slider_ax = fig.add_axes([0.25, 0.125, 0.65, 0.03], facecolor=axis_color)
H_6_slider = Slider(H_6_slider_ax, 'H_6', 0.1, 300, valinit=H_6)

# Conductances
g_K_slider_ax  = fig.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axis_color)
g_K_slider = Slider(g_K_slider_ax, 'g_K', 0.1, 100, valinit=g_K)

g_Na_slider_ax = fig.add_axes([0.25, 0.075, 0.65, 0.03], facecolor=axis_color)
g_Na_slider = Slider(g_Na_slider_ax, 'g_Na', 0.1, 500, valinit=g_Na)

g_Leak_slider_ax = fig.add_axes([0.25, 0.05, 0.65, 0.03], facecolor=axis_color)
g_Leak_slider = Slider(g_Leak_slider_ax, 'g_Leak', 0.001, 1.0, valinit=g_Leak)

# Add a set of radio buttons for changing data view
mode_ax = fig.add_axes([0.025, 0.7, 0.15, 0.20], facecolor=axis_color)
mode_radios = RadioButtons(mode_ax, ('Voltage', 'Current', 'Na_Current',
		'K_Current', 'Leak_Current', 'n', 'm', 'h'), active=0)

def mode_decode(mode_label):
	if mode_label == "Voltage":
		return 'V'
	elif mode_label == "Current":
		return 'I'
	elif mode_label == "Na_Current":
		return 'N'
	elif mode_label == "K_Current":
		return 'K'
	elif mode_label == "Leak_Current":
		return 'L'
	elif mode_label == "n":
		return 'n'
	elif mode_label == "m":
		return 'm'
	elif mode_label == "h":
		return 'h'

# Define an action for modifying the line when any slider's value changes
def sliders_on_changed(val):
	mode = mode_decode(mode_radios.value_selected)
	Vs = AP(N_1_slider.val, N_2_slider.val, N_6_slider.val, N_7_slider.val,
			M_1_slider.val, M_2_slider.val, M_6_slider.val, M_7_slider.val,
			H_3_slider.val, H_4_slider.val, H_5_slider.val, H_6_slider.val,
			g_K_slider.val, g_Na_slider.val, g_Leak_slider.val,
			stim_d_slider.val, stim_h_slider.val, stim_dim_slider.val, V_init_slider.val,
			mode)
	line.set_ydata(Vs)
	fig.canvas.draw_idle()
N_1_slider.on_changed(sliders_on_changed)
N_2_slider.on_changed(sliders_on_changed)
N_6_slider.on_changed(sliders_on_changed)
N_7_slider.on_changed(sliders_on_changed)

M_1_slider.on_changed(sliders_on_changed)
M_2_slider.on_changed(sliders_on_changed)
M_6_slider.on_changed(sliders_on_changed)
M_7_slider.on_changed(sliders_on_changed)

H_3_slider.on_changed(sliders_on_changed)
H_4_slider.on_changed(sliders_on_changed)
H_5_slider.on_changed(sliders_on_changed)
H_6_slider.on_changed(sliders_on_changed)

g_K_slider.on_changed(sliders_on_changed)
g_Na_slider.on_changed(sliders_on_changed)
g_Leak_slider.on_changed(sliders_on_changed)

stim_d_slider.on_changed(sliders_on_changed)
stim_h_slider.on_changed(sliders_on_changed)
stim_dim_slider.on_changed(sliders_on_changed)
V_init_slider.on_changed(sliders_on_changed)

def mode_on_clicked(label):
	mode = mode_decode(label)
	if mode == 'V':
		ax.set_ylim([-100, 100])
	elif mode == 'I' or mode == 'N' or mode == 'K' or mode == 'L':
		ax.set_ylim([-200, 200])
	elif mode == 'n' or mode == 'm' or mode == 'h':
		ax.set_ylim([0, 1])
	Vs = AP(N_1_slider.val, N_2_slider.val, N_6_slider.val, N_7_slider.val,
			M_1_slider.val, M_2_slider.val, M_6_slider.val, M_7_slider.val,
			H_3_slider.val, H_4_slider.val, H_5_slider.val, H_6_slider.val,
			g_K_slider.val, g_Na_slider.val, g_Leak_slider.val,
			stim_d_slider.val, stim_h_slider.val, stim_dim_slider.val, V_init_slider.val,
			mode)
	line.set_ydata(Vs)
	fig.canvas.draw_idle()
mode_radios.on_clicked(mode_on_clicked)

# Add a button for resetting the parameters
reset_button_ax = fig.add_axes([0.8, 0.015, 0.1, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
	N_1_slider.reset()
	N_2_slider.reset()
	N_6_slider.reset()
	N_7_slider.reset()

	M_1_slider.reset()
	M_2_slider.reset()
	M_6_slider.reset()
	M_7_slider.reset()

	H_3_slider.reset()
	H_4_slider.reset()
	H_5_slider.reset()
	H_6_slider.reset()

	g_K_slider.reset()
	g_Na_slider.reset()
	g_Leak_slider.reset()

	stim_d_slider.reset()
	stim_h_slider.reset()
	stim_dim_slider.reset()
	V_init_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)

plt.show()
