#python extinction_diagram.py 5.26341231209 2604310.66576 14.8993552947 0.114036183311 13.5563713727 0.661686204463 13.6316983972 0.530479078769 12.8038903036 2.45264187586 5.61771223327 2600405.84776 13.4310792758 0.604426684665 5.60684697857 2491094.08031 J0017+1307_uvb.spec

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams["font.family"] = "Times New Roman"
from matplotlib.backends.backend_pdf import PdfPages
import sys
import matplotlib.ticker as ticker


##################################LOG_TO_REAL_ERROR##################################


def log_to_real_err(log_array, log_error):
	linear_array = 10**(log_array)
	linear_err = np.abs(linear_array * (np.log(10) * log_error))
	return linear_err


##################################LOG_TO_REAL_ERROR##################################





##################################REAL_TO_LOG_ERROR##################################


def real_to_log_err(linear_array, linear_error):
	log_array = np.log10(linear_array)
	log_err = np.abs(linear_error/(linear_array*(np.log(10))))
	return log_err


##################################REAL_TO_LOG_ERROR##################################





##################################EXCITATION_TEMPERATURE_FROM_H2(J2->0)#################################


def excitation_temperature(linn0, linn2, E_k, g0, g2):	
	T_ex = (-E_k) / (np.log( (linn2/linn0) * (g0/g2) ))
	return (T_ex)


##################################EXCITATION_TEMPERATURE_FROM_H2(J2->0)#################################




############################################SLOPE############################################

def slope(x1, y1, x2, y2):
    m = (y2-y1)/(x2-x1)
    return m

############################################SLOPE############################################


############################################LINE############################################

def fit_func(x, a, b):
    return a*x + b

############################################LINE############################################









#Loading H2 rotational level column densities and errors

#'''
log_NJ0 = float(sys.argv[1])
log_NJ0_err = float(sys.argv[2])

log_NJ1 = float(sys.argv[3])
log_NJ1_err = float(sys.argv[4])

log_NJ2 = float(sys.argv[5])
log_NJ2_err = float(sys.argv[6])

log_NJ3 = float(sys.argv[7])
log_NJ3_err = float(sys.argv[8])

log_NJ4 = float(sys.argv[9])
log_NJ4_err = float(sys.argv[10])

log_NJ5 = float(sys.argv[11])
log_NJ5_err = float(sys.argv[12])

log_NJ6 = float(sys.argv[13])
log_NJ6_err = float(sys.argv[14])

log_NJ7 = float(sys.argv[15])
log_NJ7_err = float(sys.argv[16])



#Defining a name for the observed object 
str1 = str(sys.argv[17])
name_DLA = str1[:10]




#Removing unnecessary warnings created when column density is below detection limit
if (log_NJ0_err>10):
	log_NJ0_err = 1e-5

if (log_NJ1_err>10):
	log_NJ1_err = 1e-5

if (log_NJ2_err>10):
	log_NJ2_err = 1e-5

if (log_NJ3_err>10):
	log_NJ3_err = 1e-5

if (log_NJ4_err>10):
	log_NJ4_err = 1e-5

if (log_NJ5_err>10):
	log_NJ5_err = 1e-5

if (log_NJ6_err>10):
	log_NJ6_err = 1e-5

if (log_NJ7_err>10):
	log_NJ7_err = 1e-5




#Changing column density from logscale to linear scale

lin_NJ0 = 10**(log_NJ0)
lin_NJ1 = 10**(log_NJ1)
lin_NJ2 = 10**(log_NJ2)
lin_NJ3 = 10**(log_NJ3)
lin_NJ4 = 10**(log_NJ4)
lin_NJ5 = 10**(log_NJ5)
lin_NJ6 = 10**(log_NJ6)
lin_NJ7 = 10**(log_NJ7)

lin_NJ0_err = log_to_real_err(log_NJ0, log_NJ0_err)
lin_NJ1_err = log_to_real_err(log_NJ1, log_NJ1_err)
lin_NJ2_err = log_to_real_err(log_NJ2, log_NJ2_err)
lin_NJ3_err = log_to_real_err(log_NJ3, log_NJ3_err)
lin_NJ4_err = log_to_real_err(log_NJ4, log_NJ4_err)
lin_NJ5_err = log_to_real_err(log_NJ5, log_NJ5_err)
lin_NJ6_err = log_to_real_err(log_NJ6, log_NJ6_err)
lin_NJ7_err = log_to_real_err(log_NJ7, log_NJ7_err)




#Creating an array of column densities for different rotational levels of the object

H2_array = [lin_NJ0, lin_NJ1, lin_NJ2, lin_NJ3, lin_NJ4, lin_NJ5, lin_NJ6, lin_NJ7]

H2_array_err = [lin_NJ0_err, lin_NJ1_err, lin_NJ2_err, lin_NJ3_err, lin_NJ4_err, lin_NJ5_err, lin_NJ6_err, lin_NJ7_err]


#Computing the excitaiton paramters

H2_j = [0,1,2,3,4,5,6,7]

B = 85.3

k = 8.6173303e-5

g_j = np.zeros([len(H2_j)])

E_J0 = np.zeros([len(H2_j)])

y_axis_lin = np.zeros([len(H2_j)])
y_axis_err_lin = np.zeros([len(H2_j)])

for i in range(len(H2_j)):
	if (i%2==0):
		g_j[i] = (2*(H2_j[i]) + 1)*((2*0)+1)
	else:
		g_j[i] = (2*(H2_j[i]) + 1)*((2*1)+1)
		

E_K=np.array([0.,170.48,509.87,1015.08,1681.6,2503.67,3474.45,4586.31])



#######################################Gry_et_al_2002#######################################


#gry_2002 = np.array([np.nan, np.nan, 2.5e18, 1.1e17, 6.0e15, 4.5e14, np.nan, np.nan])
gry_2002_J, gry_2002_g, gry_2002_HD102065, gry_2002_HD108927, gry_2002_HD96675 = np.loadtxt('h2_excitation_diagram_files/Gry2002_2.dat', comments='#', unpack=True)

gry_2002_E_K=np.array([0.,170.48,509.87,1015.08,1681.6,2503.67])

gry_2002_HD102065_2 = np.zeros([len(gry_2002_g)])
gry_2002_HD108927_2 = np.zeros([len(gry_2002_g)])
gry_2002_HD96675_2 = np.zeros([len(gry_2002_g)])



for i in range(len(gry_2002_g)):
	gry_2002_HD102065_2[i] = (gry_2002_HD102065[i]/gry_2002_g[i])/(gry_2002_HD102065[0]/gry_2002_g[0])
	gry_2002_HD108927_2[i] = (gry_2002_HD108927[i]/gry_2002_g[i])/(gry_2002_HD108927[0]/gry_2002_g[0])
	gry_2002_HD96675_2[i] = (gry_2002_HD96675[i]/gry_2002_g[i])/(gry_2002_HD96675[0]/gry_2002_g[0])


#######################################Gry_et_al_2002#######################################


#######################################Boisse_et_al_2005#######################################


#boisse_2005 = np.array([3.2e20, 3.2e20, 1.8e19, 6.2e18, 7.1e17, 3.3e17, 4.0e15, 2.5e15])
boisse_2005_nu, boisse_2005_J, boisse_2005_g, boisse_2005_E, boisse_2005_num, boisse_2005_logN, boisse_2005_lower, boisse_2005_upper, boisse_2005_translucent_cloud, boisse_2005_C_shock, boisse_2005_hot_PDR, boisse_2005_total = np.loadtxt('h2_excitation_diagram_files/Boisse2005.dat', comments='#', unpack=True)

boisse_2005_logN_2 = np.zeros([len(boisse_2005_logN)])
boisse_2005_lower_2 = np.zeros([len(boisse_2005_logN)])
boisse_2005_upper_2 = np.zeros([len(boisse_2005_logN)])

for i in range(len(gry_2002_g)):
	boisse_2005_logN_2[i] = (boisse_2005_logN[i]/boisse_2005_g[i])/(boisse_2005_logN[0]/boisse_2005_g[0])
	boisse_2005_lower_2[i] = (boisse_2005_lower[i]/boisse_2005_g[i])/(boisse_2005_lower[0]/boisse_2005_g[0])
	boisse_2005_upper_2[i] = (boisse_2005_upper[i]/boisse_2005_g[i])/(boisse_2005_upper[0]/boisse_2005_g[0])


#######################################Boisse_et_al_2005#######################################


#######################################Noterdaeme_et_al_2015a#######################################


noterdaeme_2015 = np.array([19.84, 19.81, 17.96, 17.76, 15.88, 15.17, np.nan, np.nan])
noterdaeme_2015_err = np.array([0.09, 0.04, 0.14, 0.40, 0.26, 0.16, np.nan, np.nan])

noterdaeme_2015_err = log_to_real_err(noterdaeme_2015, noterdaeme_2015_err)
noterdaeme_2015 = 10**noterdaeme_2015

noterdaeme_2015_lin = np.zeros([len(g_j)])
noterdaeme_2015_err_lin = np.zeros([len(g_j)])

for i in range(len(g_j)):
	noterdaeme_2015_lin[i] = (noterdaeme_2015[i]/g_j[i])/(noterdaeme_2015[0]/g_j[0])
	noterdaeme_2015_err_lin[i] = np.abs(noterdaeme_2015_lin[i]) * np.sqrt( ((noterdaeme_2015_err[i]/g_j[i])/(noterdaeme_2015[i]/g_j[i]))**2 + ((noterdaeme_2015_err[0]/g_j[0])/(noterdaeme_2015[0]/g_j[0]))**2 )


noterdaeme_2015_err = np.zeros([len(noterdaeme_2015_err_lin)])
for i in range(len(noterdaeme_2015_err_lin)):
	noterdaeme_2015_err[i] = real_to_log_err(noterdaeme_2015_lin[i], noterdaeme_2015_err_lin[i])


#######################################Noterdaeme_et_al_2015a#######################################


#######################################Balashev_et_al_2017#######################################


balashev_2017 = np.array([20.71, 21.06, 19.59, 18.74, 17.23, 16.23, 14.95, 14.89])
balashev_2017_err = np.array([0.02, 0.02, 0.02, 0.02, 0.13, 0.04, 0.02, 0.02])

balashev_2017_err = log_to_real_err(balashev_2017, balashev_2017_err)
balashev_2017 = 10**balashev_2017

balashev_2017_lin = np.zeros([len(g_j)])
balashev_2017_err_lin = np.zeros([len(g_j)])

for i in range(len(g_j)):
	balashev_2017_lin[i] = (balashev_2017[i]/g_j[i])/(balashev_2017[0]/g_j[0])
	balashev_2017_err_lin[i] = np.abs(balashev_2017_lin[i]) * np.sqrt( ((balashev_2017_err[i]/g_j[i])/(balashev_2017[i]/g_j[i]))**2 + ((balashev_2017_err[0]/g_j[0])/(balashev_2017[0]/g_j[0]))**2 )




balashev_2017_err = np.zeros([len(balashev_2017_err_lin)])
for i in range(len(balashev_2017_err_lin)):
	balashev_2017_err[i] = real_to_log_err(balashev_2017_lin[i], balashev_2017_err_lin[i])


#######################################Balashev_et_al_2017#######################################


for i in range(len(H2_j)):
	y_axis_lin[i] = (H2_array[i]/g_j[i])/(H2_array[0]/g_j[0])
	y_axis_err_lin[i] = np.abs(y_axis_lin[i]) * np.sqrt( ((H2_array_err[i]/g_j[i])/(H2_array[i]/g_j[i]))**2 + ((H2_array_err[0]/g_j[0])/(H2_array[0]/g_j[0]))**2 )


x_axis = E_K

y_axis_err = np.zeros([len(y_axis_err_lin)])

for i in range(len(y_axis_err_lin)):
	y_axis_err[i] = real_to_log_err(y_axis_lin[i], y_axis_err_lin[i])





#Estimating the gas temperature from the slope fitted to lower rotational levels

#######################################J=0 to J=2#######################################


count = 2
y_axis = np.log10(y_axis_lin)
#y_axis = y_axis_lin

xnew = np.array(x_axis[0:(count+1)])
ynew = np.array(y_axis[0:(count+1)])

popt, pcov  = curve_fit(fit_func, xnew, ynew)
anew = popt[0]
bnew = popt[1]
perr = np.sqrt(np.diag(pcov))
anew_err = perr[0]
bnew_err = perr[1]

a1 = anew - anew_err
a2 = anew + anew_err
b1 = bnew - bnew_err
b2 = bnew + bnew_err

slope = np.log(10**(anew))
slope_err = np.log(10**(anew_err))


T_02_fit = (-1/slope)
T_02_fit_err = np.abs(T_02_fit * (-slope_err/slope))


#######################################J=0 to J=2#######################################



#Estimating the excitation caused from radiation

#######################################J=3 to J=4#######################################


#count = 2
xnew2 = np.array(x_axis[2:5])
ynew2 = np.array(y_axis[2:5])


popt2, pcov2  = curve_fit(fit_func, xnew2, ynew2)
anew2 = popt2[0]
bnew2 = popt2[1]
perr2 = np.sqrt(np.diag(pcov2))
anew_err2 = perr2[0]
bnew_err2 = perr2[1]

a1_2 = anew2 - anew_err2
a2_2 = anew2 + anew_err2
b1_2 = bnew2 - bnew_err2
b2_2 = bnew2 + bnew_err2

slope2 = np.log(10**(anew2))
slope_err2 = np.log(10**(anew_err2))

T_34_fit = (-1/slope2)
T_34_fit_err = np.abs(T_34_fit * (-slope_err2/slope2))


#######################################J=3 to J=4#######################################







#Plotting

size_of_font = 25
fig = plt.figure(figsize=(9, 9), dpi=300,)
ax = fig.add_subplot(111)

ax.plot(gry_2002_E_K, np.log10(gry_2002_HD102065_2), 'k.', marker="s", markersize=size_of_font/2.0, markerfacecolor='white', label='diffuse MW ISM', zorder=1, rasterized=True)
ax.plot(boisse_2005_E, np.log10(boisse_2005_logN_2), 'k.', marker="+", markersize=size_of_font/2.0, label='HD 34078', zorder=2, rasterized=True)
ax.errorbar(E_K, np.log10(noterdaeme_2015_lin), yerr=noterdaeme_2015_err, color='green', alpha=0.8, fmt='o', ecolor='g', elinewidth=1.0, capsize=1.0, capthick=1, label="J2140-0321", zorder=3, rasterized=True)
ax.errorbar(E_K, np.log10(balashev_2017_lin), yerr=balashev_2017_err, color='blue', alpha=0.8, fmt='o', ecolor='b', elinewidth=1.0, capsize=1.0, capthick=1, label="J0843+0221", zorder=4, rasterized=True)
ax.errorbar(x_axis[0:7], np.log10(y_axis_lin[0:7]), yerr=y_axis_err[0:7], color='red', alpha=1.0, fmt='o', ecolor='r', elinewidth=1.0, capsize=1.0, capthick=1, label=str(name_DLA), zorder=5, rasterized=True)
ax.plot(x_axis[0:5], fit_func(x_axis, *popt)[0:5], 'r--', zorder=6, rasterized=True)
#ax.text(1500, -8.5, r'$\rm T_{01}\, \simeq\, T_{02}\, $ =%2d$\rm \pm$%2dK' %(T_02_fit, T_02_fit_err), fontsize=size_of_font, bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
ax.text(1500, -8.5, r'$T_{01}\, \simeq\, T_{02}\, $ =%2d$\rm \pm$%2dK' %(T_02_fit, T_02_fit_err), fontsize=size_of_font, rasterized=True)

#plt.ylim(10**(-10), 10**(1))
#plt.xlim(-200,5600)
#plt.yscale('log')
ax.set_ylim(-10, 0.9)
ax.set_xlim(-200,5600)
ax.tick_params(axis='both', which='major', direction='in', length=size_of_font/2, labelsize=size_of_font)
ax.tick_params(axis='both', which='minor', direction='in', length=size_of_font/4)

ax.legend(fontsize=size_of_font)
ax.set_ylabel(r'log$N_i/g_i$ - log$N_0/g_0$', fontsize=size_of_font)
ax.set_xlabel(r'Energy, $\rm cm^{-1}$', fontsize=size_of_font)

from matplotlib.ticker import MaxNLocator

ax.xaxis.set_minor_locator(MaxNLocator(29))
ax.yaxis.set_minor_locator(MaxNLocator(22))
ax.xaxis.set_major_locator(MaxNLocator(6))
ax.yaxis.set_major_locator(MaxNLocator(11))

#plt.subplots_adjust(right = 0.0, left  = 0.0, bottom=0.0, top=0.0, wspace=0.0, hspace=0.0)
plt.subplots_adjust(wspace=0.0, hspace=0.0)


plt.title(name_DLA, fontsize=size_of_font, rasterized=True)

name_fig = str(sys.argv[17])
plt.tight_layout()


#Changing the tick params
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
ax.xaxis.set_major_locator(ticker.MultipleLocator(1000))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(200))
ax.tick_params(axis = 'both', which = 'major', direction='in', length=15, width=2, colors='k')
ax.tick_params(axis = 'both', which = 'minor', direction='in', length=7, width=1, colors='k')


#Saving the plot
plt.savefig(name_fig)
print ('Extinction Diagram Figure Saved....')
#plt.show()
