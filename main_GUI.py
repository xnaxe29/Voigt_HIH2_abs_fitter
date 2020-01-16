from __future__ import print_function, division

#numpy
import numpy
import numpy as np
import numpy.polynomial.chebyshev as cheb
import numpy.polynomial.polynomial as poly
import numpy.ma as ma
from numpy.polynomial import polynomial as P

#scipy
import scipy
import scipy as sp
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy.optimize import fmin
from scipy import integrate
from scipy.integrate import quad
from scipy import signal

#matplotlib
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons, RectangleSelector
#from matplotlib.widgets import TextBox
import matplotlib.patches as mpatches

#pyastronomy
from PyAstronomy import pyasl
from PyAstronomy.pyaC import pyaErrors as PE
from PyAstronomy import funcFit as fuf

#others
import sys
from pathlib import Path
import os
import os.path
import pyfits
import csv as csv
import itertools
from tabulate import tabulate

#Custom made functions used by this code
from relevant_functions import *

#speed of light in cm/s
c = 2.99792e10



initial_guesses = np.genfromtxt('initial_guess.dat', dtype=[('mystring','S30')], comments='#')

file_name1 = str(initial_guesses[0][0])
z = float(initial_guesses[1][0])
z_HI = float(initial_guesses[1][0])
z_H2 = float(initial_guesses[2][0])
z_qso = float(initial_guesses[3][0])
fwhm_vel = float(initial_guesses[4][0])
resolution = (c/1e5)/fwhm_vel
region1_start = float(initial_guesses[5][0])
region1_end = float(initial_guesses[6][0])
region2_start = float(initial_guesses[7][0])
region2_end = float(initial_guesses[8][0])
region3_start = float(initial_guesses[9][0])
region3_end = float(initial_guesses[10][0])

logN_HI = float(initial_guesses[11][0])
logN_H2J0 = float(initial_guesses[12][0])
logN_H2J1 = float(initial_guesses[13][0])
logN_H2J2 = float(initial_guesses[14][0])
logN_H2J3 = float(initial_guesses[15][0])
logN_H2J4 = float(initial_guesses[16][0])
logN_H2J5 = float(initial_guesses[17][0])
logN_H2J6 = float(initial_guesses[18][0])
logN_H2J7 = float(initial_guesses[19][0])
logN_HDJ0 = float(initial_guesses[20][0])
logN_HDJ1 = float(initial_guesses[21][0])
logN_HDJ2 = float(initial_guesses[22][0])

b_val = float(initial_guesses[23][0])
b_HI = float(initial_guesses[23][0])
b_H2 = float(initial_guesses[24][0])

fix_HI_continuum = int(initial_guesses[25][0])
fix_H2_continuum = int(initial_guesses[26][0])
fix_H2_Redshift = int(initial_guesses[27][0])
fix_H2J0_col_den = int(initial_guesses[28][0])
fix_H2J2_col_den = int(initial_guesses[29][0])
fix_HD_col_den = int(initial_guesses[30][0])

window_HI = float(initial_guesses[31][0])
window_H2 = float(initial_guesses[32][0])
window_HD = float(initial_guesses[33][0])

size_of_font = int(initial_guesses[34][0])

#Atomic Database Filename
atomic_database_file = str(initial_guesses[35][0])



font = {'family' : 'times new roman',
        #'weight' : 'bold',
        'size'   : size_of_font}

matplotlib.rc('font', **font)




#
#Choose the ions that you would like to fit. The default is HI. See the commented line just below for the correct format to request multiple species fitting. Please not that the interactive GUI is just created for HI, H2 and HD ions.
#
#atomic_string_name = [HI]
atomic_string_name = ['HI', 'H2', 'HD']
#
atom_name = atomic_string_name
#






#Loading the file and defining objects for continuum and fitting lines
wave, flux, flux_err = np.loadtxt(file_name1, unpack=True)
cont2 = []
fit = []

#Defining output file names
rest = file_name1.split('.')[0]
file2 = 'plotting_files/' + str(rest) + '_selected_data_H2.txt'
file3 = 'plotting_files/' + str(rest) + '_continuum_points_H2.txt'
file4 = 'plotting_files/' + str(rest) + '_selected_data2_H2.spec'
file5 = 'plotting_files/' + str(rest) + '_selected_data2_fitted_HI.spec'
file6 = 'plotting_files/' + str(rest) + '_selected_data2_fitted_H2.spec'
file7 = 'plotting_files/' + str(rest) + '_fitted_combined.spec'
file8 = 'plotting_files/' + str(rest) + '_selected_combined.spec'
file9 = 'plotting_files/' + str(rest) + '_fit_params.txt'
file10 = 'plotting_files/' + str(rest) + '_HI_H2_excel.txt'





#Loading the quasar template
qso_template_wave_rest, qso_template_flux, qso_template_flux_err = np.loadtxt('quasar_template_selsing_2015.dat', comments='#', unpack=True)
#Flux calibration parameter (set to 1 by default assuming that the quasar flux would be on the same levels as the template)
flux_red = 1.0
























##########################################CODE_STARTS_HERE#########################################


######################KEY_PRESS_EVENTS##################################


def on_key(event):
	if ((event.key).lower() == 'a'):
		global ix, iy
    		ix, iy = event.xdata, event.ydata
    		global coordsx, coordsy
    		coordsx = np.append(coordsx, (ix))
    		coordsy = np.append(coordsy, (iy))

	elif ((event.key).lower() == 'r'):
    		global ix2, iy2
    		ix2, iy2 = event.xdata, event.ydata
		ixnew = find_nearest(coordsx, ix2)
		ixnew2 = np.where(coordsx==ixnew)	
    		coordsx = np.delete(coordsx, ixnew2)
    		coordsy = np.delete(coordsy, ixnew2)

	elif ((event.key).lower() == 'm'):
    		global ix3, iy3
    		ix3, iy3 = event.xdata, event.ydata
		ixnew3 = find_nearest(coordsx, ix3)
		ixnew4 = np.where(coordsx==ixnew3)
    		coordsx[ixnew4] = ix3
		coordsy[ixnew4] = iy3




	'''
	points = zip(coordsx, coordsy)
	points = sorted(points, key=lambda point: point[0])
	x1, y1 = zip(*points)
	new_length = len(wave)
	l1 = np.searchsorted(wave, min(x1))
	l2 = np.searchsorted(wave, max(x1))
	new_x = []
	new_y = []
	global contdata
	if (len(x1)>3):
		new_x = np.linspace(min(x1), max(x1), (l2-l1))
	    	new_y = sp.interpolate.splrep(x1, y1)
    	    	contdata[l1:l2] = sp.interpolate.splev(new_x, new_y, der=0)

	'''


	points = zip(coordsx, coordsy)
	points = sorted(points, key=lambda point: point[0])
	x1, y1 = zip(*points)
	new_length = len(wave)
	l1 = np.searchsorted(wave, min(x1))
	l2 = np.searchsorted(wave, max(x1))
	new_x = []
	new_y = []
	global contdata
	global contzz
	contzz = np.zeros([(len(wave))])
	if (len(x1)>3):
		new_x = np.linspace(min(x1), max(x1), (l2-l1))
		#new_y = sp.interpolate.interp1d(x1, y1, kind='cubic')(new_x)
		new_y = sp.interpolate.splrep(x1, y1)
		#cont[l1:l2] = new_y
		contzz = sp.interpolate.splev(wave, new_y, der=0)


	contdata = contzz
   
	line5.set_data(x1, y1)
	line6.set_data(wave, contdata)
	fig.canvas.draw()
    
   
	#Disconnect if you press 'x'
	if event.key == 'x':
		fig.canvas.mpl_disconnect(cid)
   		plt.close(1)

	return


def press(event):
	print('press', event.key)
	sys.stdout.flush()


######################KEY_PRESS_EVENTS##################################





############################CHARACTERISTIC_DEFINITIONS#####################################

def parameter_declarations2(atom_name, z_HI, z_H2, logN_HI, logN_H2J0, logN_H2J1, logN_H2J2, logN_H2J3, logN_H2J4, logN_H2J5, logN_H2J6, logN_H2J7, logN_HDJ0, logN_HDJ1, logN_HDJ2):
	name, wave, osc, tau = atomic_list_extractor(atom_name)
	logN = np.zeros([len(wave)])
	z = np.zeros([len(wave)])
	b = np.zeros([len(wave)])	
	for i in range(len(wave)):
		if ('HI' in str(name[i])):
			logN[i] = logN_HI
			z[i] = z_HI
			b[i] = b_HI
		elif ('H2J0' in str(name[i])):
			logN[i] = logN_H2J0
			z[i] = z_H2
			b[i] = b_H2
		elif ('H2J1' in str(name[i])):
			logN[i] = logN_H2J1
			z[i] = z_H2
			b[i] = b_H2
		elif ('H2J2' in str(name[i])):
			logN[i] = logN_H2J2
			z[i] = z_H2
			b[i] = b_H2
		elif ('H2J3' in str(name[i])):
			logN[i] = logN_H2J3
			z[i] = z_H2
			b[i] = b_H2
		elif ('H2J4' in str(name[i])):
			logN[i] = logN_H2J4
			z[i] = z_H2
			b[i] = b_H2
		elif ('H2J5' in str(name[i])):
			logN[i] = logN_H2J5
			z[i] = z_H2
			b[i] = b_H2
		elif ('H2J6' in str(name[i])):
			logN[i] = logN_H2J6
			z[i] = z_H2
			b[i] = b_H2
		elif ('H2J7' in str(name[i])):
			logN[i] = logN_H2J7
			z[i] = z_H2
			b[i] = b_H2
		elif ('HDJ0' in str(name[i])):
			logN[i] = logN_HDJ0
			z[i] = z_H2
			b[i] = b_H2
		elif ('HDJ1' in str(name[i])):
			logN[i] = logN_HDJ1
			z[i] = z_H2
			b[i] = b_H2
		elif ('HDJ2' in str(name[i])):
			logN[i] = logN_HDJ2
			z[i] = z_H2
			b[i] = b_H2
	return (name, wave, osc, tau, logN, z, b)
	


############################CHARACTERISTIC_DEFINITIONS#####################################




####################################################################################
##############################FITTING_FUNCTIONS#####################################
####################################################################################

def profile_func_complete(wave, res, cheb_order_HI, cheb_order_H2, cheb_order_H2_2, *pars):
	pars_new = []
	for i in pars:
		pars_new = np.append(pars_new, i)

	cheb_order = cheb_order_HI+cheb_order_H2+cheb_order_H2_2
	cheb_coeff = pars_new[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):]
	cheb_coeff_HI = cheb_coeff[0:cheb_order_HI]
	cheb_coeff_H2 = cheb_coeff[cheb_order_HI:cheb_order_HI+cheb_order_H2]
	cheb_coeff_H2_2 = cheb_coeff[cheb_order_HI+cheb_order_H2:cheb_order_HI+cheb_order_H2+cheb_order_H2_2]
	cont_new_HI = chebyshev_disp(xdata2, cheb_coeff_HI/1e5)
	cont_new_H2 = chebyshev_disp(xdata3, cheb_coeff_H2/1e5)
	cont_new_H2_2 = chebyshev_disp(xdata4, cheb_coeff_H2_2/1e5)
	
	logN_HI, logN_H2J0, logN_H2J1, logN_H2J2, logN_H2J3, logN_H2J4, logN_H2J5, logN_H2J6, logN_H2J7, logN_HDJ0, logN_HDJ1, logN_HDJ2, z_HI, z_H2 = (pars_new[:(len(pars_new)-cheb_order)])
	name, wave_atoms, osc, tau, logN, z_ind, b_ind = parameter_declarations2(atom_name, z_HI/1e5, z_H2/1e5, logN_HI, logN_H2J0, logN_H2J1, logN_H2J2, logN_H2J3, logN_H2J4, logN_H2J5, logN_H2J6, logN_H2J7, logN_HDJ0, logN_HDJ1, logN_HDJ2)
	profile_HI = group_voigt2_H2(wave_atoms, osc, tau, xdata2, logN, b_ind, z_ind, resolution)*cont_new_HI
	profile_H2 = group_voigt2_H2(wave_atoms, osc, tau, xdata3, logN, b_ind, z_ind, resolution)*cont_new_H2
	profile_H2_2 = group_voigt2_H2(wave_atoms, osc, tau, xdata4, logN, b_ind, z_ind, resolution)*cont_new_H2_2
	
	return (profile_HI, profile_H2, profile_H2_2)





def profile_func_complete_new(wave, res, cheb_order_HI, cheb_order_H2, cheb_order_H2_2, *pars):
	pars_new = []
	for i in pars:
		pars_new = np.append(pars_new, i)

	cheb_order = cheb_order_HI+cheb_order_H2+cheb_order_H2_2
	cheb_coeff = pars_new[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):]
	cheb_coeff_HI = cheb_coeff[0:cheb_order_HI]
	cheb_coeff_H2 = cheb_coeff[cheb_order_HI:cheb_order_HI+cheb_order_H2]
	cheb_coeff_H2_2 = cheb_coeff[cheb_order_HI+cheb_order_H2:cheb_order_HI+cheb_order_H2+cheb_order_H2_2]
	cont_new_HI = chebyshev_disp(xdata2, cheb_coeff_HI/1e5)
	cont_new_H2 = chebyshev_disp(xdata3, cheb_coeff_H2/1e5)
	cont_new_H2_2 = chebyshev_disp(xdata4, cheb_coeff_H2_2/1e5)
	
	logN_HI, logN_H2J0, logN_H2J1, logN_H2J2, logN_H2J3, logN_H2J4, logN_H2J5, logN_H2J6, logN_H2J7, logN_HDJ0, logN_HDJ1, logN_HDJ2, z_HI, z_H2, b_val_new = (pars_new[:(len(pars_new)-cheb_order)])
	name, wave_atoms, osc, tau, logN, z_ind, b_ind = parameter_declarations2(atom_name, z_HI/1e5, z_H2/1e5, logN_HI, logN_H2J0, logN_H2J1, logN_H2J2, logN_H2J3, logN_H2J4, logN_H2J5, logN_H2J6, logN_H2J7, logN_HDJ0, logN_HDJ1, logN_HDJ2)
	b_ind[:] = b_val_new
	profile_HI = group_voigt2_H2(wave_atoms, osc, tau, xdata2, logN, b_ind, z_ind, resolution)*cont_new_HI
	profile_H2 = group_voigt2_H2(wave_atoms, osc, tau, xdata3, logN, b_ind, z_ind, resolution)*cont_new_H2
	profile_H2_2 = group_voigt2_H2(wave_atoms, osc, tau, xdata4, logN, b_ind, z_ind, resolution)*cont_new_H2_2
	
	return (profile_HI, profile_H2, profile_H2_2)






def fitting_func_complete(wave, *pars_complete):
	pars_new = []
	for i in pars_complete:
		pars_new = np.append(pars_new, i)

	profile_HI, profile_H2, profile_H2_2 = profile_func_complete(wave, res, cheb_order_HI, cheb_order_H2, cheb_order_H2_2, *pars_new)
	profile_HI_selected = profile_HI[numpy.logical_not(numpy.isnan(ynew2))]
	profile_H2_selected = profile_H2[numpy.logical_not(numpy.isnan(ynew3))]
	profile_H2_2_selected = profile_H2_2[numpy.logical_not(numpy.isnan(ynew4))]
	total_profile = []
	total_profile = np.append(total_profile, profile_HI_selected)
	total_profile = np.append(total_profile, profile_H2_selected)
	total_profile = np.append(total_profile, profile_H2_2_selected)

	return (total_profile)





def fitting_func_complete_new(wave, *pars_complete):
	pars_new = []
	for i in pars_complete:
		pars_new = np.append(pars_new, i)

	profile_HI, profile_H2, profile_H2_2 = profile_func_complete_new(wave, res, cheb_order_HI, cheb_order_H2, cheb_order_H2_2, *pars_new)
	profile_HI_selected = profile_HI[numpy.logical_not(numpy.isnan(ynew2))]
	profile_H2_selected = profile_H2[numpy.logical_not(numpy.isnan(ynew3))]
	profile_H2_2_selected = profile_H2_2[numpy.logical_not(numpy.isnan(ynew4))]
	total_profile = []
	total_profile = np.append(total_profile, profile_HI_selected)
	total_profile = np.append(total_profile, profile_H2_selected)
	total_profile = np.append(total_profile, profile_H2_2_selected)

	return (total_profile)





def viewing_func_complete3(wave, *pars_complete):
	pars_new = []
	for i in pars_complete:
		pars_new = np.append(pars_new, i)

	profile_HI, profile_H2, profile_H2_2 = profile_func_complete(wave, res, cheb_order_HI, cheb_order_H2, cheb_order_H2_2, *pars_new)
	profile_HI_selected = profile_HI
	profile_H2_selected = profile_H2
	profile_H2_2_selected = profile_H2_2
	cheb_order = cheb_order_HI+cheb_order_H2+cheb_order_H2_2
	cheb_coeff = pars_new[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):]
	cheb_coeff_HI = cheb_coeff[0:cheb_order_HI]
	cheb_coeff_H2 = cheb_coeff[cheb_order_HI:cheb_order_HI+cheb_order_H2]
	cheb_coeff_H2_2 = cheb_coeff[cheb_order_HI+cheb_order_H2:cheb_order_HI+cheb_order_H2+cheb_order_H2_2]
	cont_new_HI = chebyshev_disp(xdata2, cheb_coeff_HI/1e5)
	cont_new_H2 = chebyshev_disp(xdata3, cheb_coeff_H2/1e5)
	cont_new_H2_2 = chebyshev_disp(xdata4, cheb_coeff_H2_2/1e5)
	total_profile = []
	total_profile = np.append(total_profile, profile_HI_selected)
	total_profile = np.append(total_profile, profile_H2_selected)
	total_profile = np.append(total_profile, profile_H2_2_selected)
	total_cont = []
	total_cont = np.append(total_cont, cont_new_HI)
	total_cont = np.append(total_cont, cont_new_H2)
	total_cont = np.append(total_cont, cont_new_H2_2)

	return (total_profile, total_cont)





def fit_curvefit_complete(p0, datax, datay, yerror, function, method_str, **kwargs):
	
	bounds_lower = np.zeros([len(p0)])
	bounds_upper = np.zeros([len(p0)])
	bounds_lower[:] = -np.inf
	bounds_upper[:] = np.inf
	cheb_order = cheb_order_HI+cheb_order_H2+cheb_order_H2_2
	

	if (fix_HI_continuum==1 and fix_H2_continuum==1):
		bounds_lower[-(cheb_order):-(cheb_order_H2+cheb_order_H2_2)] = p0[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):-(cheb_order_H2+cheb_order_H2_2)]-0.1
		bounds_upper[-(cheb_order):-(cheb_order_H2+cheb_order_H2_2)] = p0[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):-(cheb_order_H2+cheb_order_H2_2)]+0.1
		bounds_lower[-(cheb_order_H2+cheb_order_H2_2):] = p0[-(cheb_order_H2+cheb_order_H2_2):]-0.1
		bounds_upper[-(cheb_order_H2+cheb_order_H2_2):] = p0[-(cheb_order_H2+cheb_order_H2_2):]+0.1
		print ('HI and H2 continuum bounded')

	elif (fix_HI_continuum==1 and fix_H2_continuum==0):
		bounds_lower[-(cheb_order):-(cheb_order_H2+cheb_order_H2_2)] = p0[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):-(cheb_order_H2+cheb_order_H2_2)]-0.1
		bounds_upper[-(cheb_order):-(cheb_order_H2+cheb_order_H2_2)] = p0[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):-(cheb_order_H2+cheb_order_H2_2)]+0.1
		print ('HI continuum bounded')

	elif (fix_HI_continuum==0 and fix_H2_continuum==1):
		bounds_lower[-(cheb_order_H2+cheb_order_H2_2):] = p0[-(cheb_order_H2+cheb_order_H2_2):]-0.1
		bounds_upper[-(cheb_order_H2+cheb_order_H2_2):] = p0[-(cheb_order_H2+cheb_order_H2_2):]+0.1
		print ('H2 continuum bounded')


	if (fix_H2_Redshift==1):
		bounds_lower[-(cheb_order+1)] = p0[-(cheb_order+1)] - 1e-8
		bounds_upper[-(cheb_order+1)] = p0[-(cheb_order+1)] + 1e-8
		print ('H2 Redshift bounded')


	#Column density bounds (logscale)
	bounds_lower[:12] = 5.	
	bounds_upper[:12] = 25.

	if (fix_H2J0_col_den==1):
		bounds_lower[1] = p0[1] - 1e-3
		bounds_upper[1] = p0[1] + 1e-3
		bounds_lower[2] = p0[2] - 1e-3
		bounds_upper[2] = p0[2] + 1e-3
		print ('H2J0/H2J1 Column Density bounded')


	if (fix_H2J2_col_den==1):
		bounds_lower[3:9] = p0[3:9] - 1e-3
		bounds_upper[3:9] = p0[3:9] + 1e-3
		print ('H2J2+ Column Density bounded')


	if (fix_HD_col_den==1):
		bounds_lower[9:12] = p0[9:12] - 1e-3
		bounds_upper[9:12] = p0[9:12] + 1e-3
		print ('HD Column Density bounded')


	pfit, pcov = curve_fit(function,datax,datay,p0=p0, bounds=((bounds_lower), (bounds_upper)), sigma=yerror, method=method_str, **kwargs)


	error = [] 
	for i in range(len(pfit)):
		try:
			error.append(np.absolute(pcov[i][i])**0.5)
		except:
			error.append(0.00)
	pfit_curvefit = pfit
	perr_curvefit = np.array(error)

	return pfit_curvefit, perr_curvefit


####################################################################################



def fit_curvefit_complete_new(p0, datax, datay, yerror, function, method_str, **kwargs):
	
	bounds_lower = np.zeros([len(p0)])
	bounds_upper = np.zeros([len(p0)])
	bounds_lower[:] = -np.inf
	bounds_upper[:] = np.inf
	cheb_order = cheb_order_HI+cheb_order_H2+cheb_order_H2_2
	

	if (fix_HI_continuum==1 and fix_H2_continuum==1):
		bounds_lower[-(cheb_order):-(cheb_order_H2+cheb_order_H2_2)] = p0[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):-(cheb_order_H2+cheb_order_H2_2)]-0.1
		bounds_upper[-(cheb_order):-(cheb_order_H2+cheb_order_H2_2)] = p0[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):-(cheb_order_H2+cheb_order_H2_2)]+0.1
		bounds_lower[-(cheb_order_H2+cheb_order_H2_2):] = p0[-(cheb_order_H2+cheb_order_H2_2):]-0.1
		bounds_upper[-(cheb_order_H2+cheb_order_H2_2):] = p0[-(cheb_order_H2+cheb_order_H2_2):]+0.1
		print ('HI and H2 continuum bounded')

	elif (fix_HI_continuum==1 and fix_H2_continuum==0):
		bounds_lower[-(cheb_order):-(cheb_order_H2+cheb_order_H2_2)] = p0[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):-(cheb_order_H2+cheb_order_H2_2)]-0.1
		bounds_upper[-(cheb_order):-(cheb_order_H2+cheb_order_H2_2)] = p0[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):-(cheb_order_H2+cheb_order_H2_2)]+0.1
		print ('HI continuum bounded')

	elif (fix_HI_continuum==0 and fix_H2_continuum==1):
		bounds_lower[-(cheb_order_H2+cheb_order_H2_2):] = p0[-(cheb_order_H2+cheb_order_H2_2):]-0.1
		bounds_upper[-(cheb_order_H2+cheb_order_H2_2):] = p0[-(cheb_order_H2+cheb_order_H2_2):]+0.1
		print ('H2 continuum bounded')


	if (fix_H2_Redshift==1):
		bounds_lower[-(cheb_order+1)] = p0[-(cheb_order+1)] - 1e-8
		bounds_upper[-(cheb_order+1)] = p0[-(cheb_order+1)] + 1e-8
		print ('H2 Redshift bounded')



	

	#Column density bounds (logscale)
	bounds_lower[:12] = 5.	
	bounds_upper[:12] = 25.

	if (fix_H2J0_col_den==1):
		bounds_lower[1] = p0[1] - 1e-3
		bounds_upper[1] = p0[1] + 1e-3
		bounds_lower[2] = p0[2] - 1e-3
		bounds_upper[2] = p0[2] + 1e-3
		print ('H2J0/H2J1 Column Density bounded')


	if (fix_H2J2_col_den==1):
		bounds_lower[3:9] = p0[3:9] - 1e-3
		bounds_upper[3:9] = p0[3:9] + 1e-3
		print ('H2J2+ Column Density bounded')


	if (fix_HD_col_den==1):
		bounds_lower[9:12] = p0[9:12] - 1e-3
		bounds_upper[9:12] = p0[9:12] + 1e-3
		print ('HD Column Density bounded')



	bounds_lower[14] = p0[14]-1e-5	
	bounds_upper[14] = p0[14]+1e-5



	

	pfit, pcov = curve_fit(function,datax,datay,p0=p0, bounds=((bounds_lower), (bounds_upper)), sigma=yerror, method=method_str, **kwargs)


	error = [] 
	for i in range(len(pfit)):
		try:
			error.append(np.absolute(pcov[i][i])**0.5)
		except:
			error.append(0.00)
	pfit_curvefit = pfit
	perr_curvefit = np.array(error)

	return pfit_curvefit, perr_curvefit


####################################################################################
##############################FITTING_FUNCTIONS#####################################
####################################################################################








#Declaring multiple data variable (for ease of use) 

xdata = wave
ydata = flux
ydata_err = flux_err
contdata = np.zeros([len(wave)])
name, wave_atoms, osc, tau, logN, z_ind, b_ind = parameter_declarations2(atom_name, z_HI, z_H2, logN_HI, logN_H2J0, logN_H2J1, logN_H2J2, logN_H2J3, logN_H2J4, logN_H2J5, logN_H2J6, logN_H2J7, logN_HDJ0, logN_HDJ1, logN_HDJ2)
profile = group_voigt2_H2(wave_atoms, osc, tau, xdata, logN, b_ind, z_ind, resolution)*contdata




#Plotting

fig = plt.figure(1)
ax = fig.add_subplot(111)
coordsx = np.array([])
coordsy = np.array([])
cid = fig.canvas.mpl_connect('key_press_event', on_key)
line1, = ax.plot(xdata, ydata, 'b-', alpha=0.25, zorder=7)
line3, = ax.plot(xdata, profile, color='c', linewidth=2, alpha=0.75, zorder=8)
point, = ax.plot([], [], marker="o", color="crimson", zorder=2)
line2, = ax.plot([], [], 'r.', zorder=3)
line4, = ax.plot([], [], 'k.', zorder=1, alpha=0.1, linewidth=1)
line5, = ax.plot([],[], 'yo', zorder=5)
line6, = ax.plot([],[], 'g--', zorder=6)
line7, = ax.plot(qso_template_wave_rest*(1+z_qso), qso_template_flux, 'k-')



#Check buttons indicating the visibility of lines (model, dataset, fit and continuum) in the GUI

rax = plt.axes([0.00, 0.74, 0.1, 0.15])
check = CheckButtons(rax, ('model', 'data', 'fit', 'cont'), (True, True, True, True))
def func(label):
	if label == 'model':
		line7.set_visible(not line7.get_visible())
	elif label == 'data':
		line1.set_visible(not line1.get_visible())
	elif label == 'fit':
		line3.set_visible(not line3.get_visible())
	elif label == 'cont':
		line5.set_visible(not line5.get_visible())	
		line6.set_visible(not line6.get_visible())


	fig.canvas.draw()
check.on_clicked(func)



#Radio buttons indicating the whether you would like to vary the continuum while fitting

rax2 = plt.axes([0.905, 0.82, 0.07, 0.095])
radio = RadioButtons(rax2, ('Vary_HI_Cont', 'Fix_HI_Cont'))

def hzfunc(label):
	s0 = 0
	s1 = 1
	hzdict = {'Vary_HI_Cont': s0, 'Fix_HI_Cont': s1}
	global fix_HI_continuum
	fix_HI_continuum = hzdict[label]

	if (fix_HI_continuum==1):
		print ('HI continuum fixed')
	elif (fix_HI_continuum==0):
		print ('HI continuum free')

	
radio.on_clicked(hzfunc)



rax3 = plt.axes([0.905, 0.72, 0.07, 0.095])
radio2 = RadioButtons(rax3, ('Vary_H2_Cont', 'Fix_H2_Cont'))

def hzfunc2(label2):
	s0 = 0
	s1 = 1
	hzdict2 = {'Vary_H2_Cont': s0, 'Fix_H2_Cont': s1}
	global fix_H2_continuum
	fix_H2_continuum = hzdict2[label2]
	
	if (fix_H2_continuum==1):
		print ('H2 continuum fixed')
	elif (fix_H2_continuum==0):
		print ('H2 continuum free')

	
radio2.on_clicked(hzfunc2)






#Radio button indicating the whether you would like to vary the H2 redshift while fitting

rax4 = plt.axes([0.905, 0.62, 0.085, 0.095])
radio3 = RadioButtons(rax4, ('Vary_H2_Redshift', 'Fix_H2_Redshift'))

def hzfunc3(label3):
	s0 = 0
	s1 = 1
	hzdict3 = {'Vary_H2_Redshift': s0, 'Fix_H2_Redshift': s1}
	global fix_H2_Redshift
	fix_H2_Redshift = hzdict3[label3]
	
	if (fix_H2_Redshift==1):
		print ('H2-z fixed')
	elif (fix_H2_Redshift==0):
		print ('H2-z free')

	
radio3.on_clicked(hzfunc3)





#Radio buttons indicating the whether you would like to vary the H2 and HD column densities while fitting

rax5 = plt.axes([0.905, 0.52, 0.07, 0.095])
radio5 = RadioButtons(rax5, ('Vary_H2J0/H2J1', 'Fix_H2J0/H2J1'))

def hzfunc5(label5):
	s0 = 0
	s1 = 1
	hzdict5 = {'Vary_H2J0/H2J1': s0, 'Fix_H2J0/H2J1': s1}
	global fix_H2J0_col_den
	fix_H2J0_col_den = hzdict5[label5]

	if (fix_H2J0_col_den==1):
		print ('H2J0/H2J1 col_den fixed')
	elif (fix_H2J0_col_den==0):
		print ('H2J0/H2J1 col_den free')

	
radio5.on_clicked(hzfunc5)


rax6 = plt.axes([0.905, 0.42, 0.07, 0.095])
radio6 = RadioButtons(rax6, ('Vary_H2J2+', 'Fix_H2J2+'))

def hzfunc6(label6):
	s0 = 0
	s1 = 1
	hzdict6 = {'Vary_H2J2+': s0, 'Fix_H2J2+': s1}
	global fix_H2J2_col_den
	fix_H2J2_col_den = hzdict6[label6]

	if (fix_H2J2_col_den==1):
		print ('H2J2+ col_den fixed')
	elif (fix_H2J2_col_den==0):
		print ('H2J2+ col_den free')

	
radio6.on_clicked(hzfunc6)


rax7 = plt.axes([0.905, 0.32, 0.07, 0.095])
radio7 = RadioButtons(rax7, ('Vary_HD', 'Fix_HD'))

def hzfunc7(label7):
	s0 = 0
	s1 = 1
	hzdict7 = {'Vary_HD': s0, 'Fix_HD': s1}
	global fix_HD_col_den
	fix_HD_col_den = hzdict7[label7]

	if (fix_HD_col_den==1):
		print ('HD col_den fixed')
	elif (fix_HD_col_den==0):
		print ('HD col_den free')

	
radio7.on_clicked(hzfunc7)








#DEFINITION FOR THE 4 POINTS ON PLT.AXES 
#([1 - SIGNIFIES WHERE THE BAR WILL START ON THE LEFT (0.0 - EXTREME LEFT, 1.0 - EXTREME RIGHT(OUT OF SIGHT))])
#([2 - SIGNIFIES WHERE THE BAR WILL START FROM THE BOTTOM (0.0 - EXTREME BOTTOM, 1.0 - EXTREME TOP(OUT OF SIGHT))])
#([3 - LENGTH OF THE BAR (0.0 - NO LENGTH, 1.0 - FULL MATPLOTLIB WINDOW])
#([4 - WIDTH OF THE BAR (0.0 - ERROR(NO-WIDTH), 0.05 - FAT BAR, 0.01/0.02 - IDEAL VALUES])

#Sliders defined for varying the column densities of HI, H2 and HD.  

axcolor = 'lightgoldenrodyellow'
axlogN_HI = plt.axes([0.025, 0.010, 0.30, 0.01])
slogN_HI = Slider(axlogN_HI, 'HI', logN_HI-15, logN_HI+15, valinit=logN_HI)

axlogN_H2J0 = plt.axes([0.025, 0.025, 0.30, 0.01])
slogN_H2J0 = Slider(axlogN_H2J0, 'H2J0', logN_H2J0-15, logN_H2J0+15, valinit=logN_H2J0)

axlogN_H2J1 = plt.axes([0.025, 0.040, 0.30, 0.01])
slogN_H2J1 = Slider(axlogN_H2J1, 'H2J1', logN_H2J1-15, logN_H2J1+15, valinit=logN_H2J1)

axlogN_H2J2 = plt.axes([0.020, 0.055, 0.13, 0.01])
slogN_H2J2 = Slider(axlogN_H2J2, 'H2J2', logN_H2J2-15, logN_H2J2+15, valinit=logN_H2J2)

axlogN_H2J3 = plt.axes([0.195, 0.055, 0.13, 0.01])
slogN_H2J3 = Slider(axlogN_H2J3, 'H2J3', logN_H2J3-15, logN_H2J3+15, valinit=logN_H2J3)

axlogN_H2J4 = plt.axes([0.375, 0.010, 0.13, 0.01])
slogN_H2J4 = Slider(axlogN_H2J4, 'H2J4', logN_H2J4-15, logN_H2J4+15, valinit=logN_H2J4)

axlogN_H2J5 = plt.axes([0.55, 0.010, 0.13, 0.01])
slogN_H2J5 = Slider(axlogN_H2J5, 'H2J5', logN_H2J5-15, logN_H2J5+15, valinit=logN_H2J5)

axlogN_H2J6 = plt.axes([0.375, 0.025, 0.13, 0.01])
slogN_H2J6 = Slider(axlogN_H2J6, 'H2J6', logN_H2J6-15, logN_H2J6+15, valinit=logN_H2J6)

axlogN_H2J7 = plt.axes([0.55, 0.025, 0.13, 0.01])
slogN_H2J7 = Slider(axlogN_H2J7, 'H2J7', logN_H2J7-15, logN_H2J7+15, valinit=logN_H2J7)

axlogN_HDJ0 = plt.axes([0.375, 0.040, 0.13, 0.01])
slogN_HDJ0 = Slider(axlogN_HDJ0, 'HDJ0', logN_HDJ0-15, logN_HDJ0+15, valinit=logN_HDJ0)

axlogN_HDJ1 = plt.axes([0.55, 0.040, 0.13, 0.01])
slogN_HDJ1 = Slider(axlogN_HDJ1, 'HDJ1', logN_HDJ1-15, logN_HDJ1+15, valinit=logN_HDJ1)

axlogN_HDJ2 = plt.axes([0.375, 0.055, 0.13, 0.01])
slogN_HDJ2 = Slider(axlogN_HDJ2, 'HDJ2', logN_HDJ2-15, logN_HDJ2+15, valinit=logN_HDJ2)




#Sliders defined for varying the Redshift of the QSO, DLA-HI and DLA-H2. HD follows H2 redshfit  

axqso = plt.axes([0.065, 0.97, 0.3, 0.02])
sqso = Slider(axqso, 'QSO-Redshift', (z_qso-1), (z_qso+1), valinit=z_qso)

axHIdla = plt.axes([0.065, 0.94, 0.3, 0.02])
sHIdla = Slider(axHIdla, 'DLA-HI-Redshift', (z_HI-0.1), (z_HI+0.1), valinit=z_HI)

axH2dla = plt.axes([0.065, 0.91, 0.3, 0.02])
sH2dla = Slider(axH2dla, 'DLA-H2-Redshift', (z_H2-0.001), (z_H2+0.001), valinit=z_H2)



#Sliders defined for varying the flux scale of the QSO template, instrument resolution

axflux_red = plt.axes([0.60, 0.97, 0.3, 0.02])
sflux_red = Slider(axflux_red, 'Flux_reduction', (flux_red-0.9), (flux_red+25), valinit=flux_red)

ax_resol = plt.axes([0.60, 0.94, 0.3, 0.02])
sresol = Slider(ax_resol, 'Resolution', (100), (resolution+20000), valinit=resolution)




# Smoothing parameter defined for dataset. Please note that the fitting functions will not take the smoothing into account. This is not graphical purpose only 
global data_smooth
data_smooth = 1.0
axdata_smooth = plt.axes([0.87, 0.055, 0.1, 0.02])
sdata_smooth = Slider(axdata_smooth, 'Smoothing', (data_smooth+0), (data_smooth+101), valinit=data_smooth)




# Update values

def update(val):
	logN_HI = slogN_HI.val
	logN_H2J0 = slogN_H2J0.val
    	logN_H2J1 = slogN_H2J1.val
    	logN_H2J2 = slogN_H2J2.val
    	logN_H2J3 = slogN_H2J3.val
    	logN_H2J4 = slogN_H2J4.val
    	logN_H2J5 = slogN_H2J5.val
    	logN_H2J6 = slogN_H2J6.val
    	logN_H2J7 = slogN_H2J7.val
    	logN_HDJ0 = slogN_HDJ0.val
    	logN_HDJ1 = slogN_HDJ1.val
    	logN_HDJ2 = slogN_HDJ2.val
	flux_red = sflux_red.val
	resolution = sresol.val
	box_pts = int(sdata_smooth.val)
	z_HI = sHIdla.val
	z_H2 = sH2dla.val
   	z_qso = sqso.val    
    	global profile_updated
    	name, wave, osc, tau, logN, z_ind, b_ind = parameter_declarations2(atom_name, z_HI, z_H2, logN_HI, logN_H2J0, logN_H2J1, logN_H2J2, logN_H2J3, logN_H2J4, logN_H2J5, logN_H2J6, logN_H2J7, logN_HDJ0, logN_HDJ1, logN_HDJ2)
	
	if 'cont' in dir():
		profile_updated = group_voigt2_H2(wave, osc, tau, xdata, logN, b_ind, z_ind, resolution)*cont
	else:
		profile_updated = group_voigt2_H2(wave, osc, tau, xdata, logN, b_ind, z_ind, resolution)*contdata
    	
    	line1.set_ydata(smooth(ydata,box_pts))
	line3.set_ydata(smooth(profile_updated, box_pts)) 
    	line7.set_xdata(qso_template_wave_rest*(1+z_qso))
    	line7.set_ydata(qso_template_flux/flux_red)


fig.canvas.draw_idle()
slogN_HI.on_changed(update)
slogN_H2J0.on_changed(update)
slogN_H2J1.on_changed(update)
slogN_H2J2.on_changed(update)
slogN_H2J3.on_changed(update)
slogN_H2J4.on_changed(update)
slogN_H2J5.on_changed(update)
slogN_H2J6.on_changed(update)
slogN_H2J7.on_changed(update)
slogN_HDJ0.on_changed(update)
slogN_HDJ1.on_changed(update)
slogN_HDJ2.on_changed(update)
sqso.on_changed(update)
sHIdla.on_changed(update)
sH2dla.on_changed(update)
sflux_red.on_changed(update)
sdata_smooth.on_changed(update)
sresol.on_changed(update)






#Function to select the points which can be used for fitting

text = ax.text(0, 0, "")
xnew = xdata
ynew = np.full([len(ydata)], np.nan)
ynew_err = np.full([len(ydata_err)], np.nan)
contnew = np.full([len(contdata)], np.nan)
str1 = ''


def line_select_callback(eclick, erelease):
	global str1
	global xnew, ynew, ynew_err
	x1, y1 = eclick.xdata, eclick.ydata
	x2, y2 = erelease.xdata, erelease.ydata
	idx5 = np.searchsorted(xdata, min(x1, x2))
	idx6 = np.searchsorted(xdata, max(x1, x2))
	min_y = min(y1, y2)
	max_y = max(y1, y2)

	if (str1 == 'add'):
		for i in range(idx5, idx6):
        		if (min_y < ydata[i] < max_y):
                		xnew[i] = xdata[i]
                		ynew[i] = ydata[i]
                		ynew_err[i] = ydata_err[i]
                		contnew[i] = contdata[i]

	elif (str1 == 'rem'):
        	for i in range(idx5, idx6):
        		for j in range(len(xnew)):
                		if (min_y < ydata[i] < max_y) and (xdata[i] == xnew[j]) and (ydata[i] == ynew[j]):
                    			xnew[j] = np.NAN
                    			ynew[j] = np.NAN
                    			ynew_err[j] = np.NAN
                    			contnew[j] = np.NAN

	else:
        	print ("Function Inactive.....")

	line2.set_data(xnew, ynew)
	fig.canvas.draw()


def toggle_selector(event):
	global str1
	print(' Key pressed.')
	if event.key in ['T', 't'] and toggle_selector.RS.active:
        	print(' RectangleSelector deactivated.')
        	toggle_selector.RS.set_active(False)
	if event.key in ['Y', 'y'] and not toggle_selector.RS.active:
        	print(' RectangleSelector activated.')
        	toggle_selector.RS.set_active(True)
	if event.key in ['H', 'h'] and toggle_selector.RS.active:
        	print('Add function activated')
        	str1 = 'add'
        	toggle_selector.RS.set_active(True)
	if event.key in ['J', 'j'] and toggle_selector.RS.active:
        	print('Remove function activated')
        	str1 = 'rem'
        	toggle_selector.RS.set_active(True)


toggle_selector.RS = RectangleSelector(ax, line_select_callback, drawtype='box', useblit=False, button=[1], minspanx=5, minspany=5, spancoords='pixels', interactive=True)

plt.connect('key_press_event', toggle_selector)










################LOAD_AND_SAVE_SELECTION_POINTS####################

save_selec_ax = plt.axes([0.710, 0.025, 0.1, 0.02])
button6 = Button(save_selec_ax, 'save_data')
def save_selec_ax(event):
	x2 = np.nan_to_num(xnew)
	y2 = np.nan_to_num(ynew)
	err2 = np.nan_to_num(ynew_err)
	cont2 = np.nan_to_num(contnew)
	data2 = np.transpose(np.array([x2, y2, err2, cont2]))
	if (os.path.isfile(file2)):
		save_prompt2=raw_input("File already exists. Overwrite?(y/n) : ")
        	if (save_prompt2=='y'):
        		np.savetxt(file2, data2)
        		print ('saved')
        	else:
        		print ('data file not saved')
	else:
       		np.savetxt(file2, data2)
       		print ('saved')
button6.on_clicked(save_selec_ax)


get_selec_ax = plt.axes([0.710, 0.0035, 0.1, 0.02])
button4 = Button(get_selec_ax, 'get_data')
def get_selec_ax(event):
	x, y, err2, cont2 = np.loadtxt(file2, unpack=True)
	x[x==0] = np.nan
	y[y==0] = np.nan
	err2[err2==0] = np.nan
	cont2[cont2==0] = np.nan
	global xnew, ynew, err, contnew
	xnew = x
	ynew = y
	err = err2
	contnew = cont2
	line2.set_data(x, y)
	fig.canvas.draw()   
button4.on_clicked(get_selec_ax)

################LOAD_AND_SAVE_SELECTION_POINTS####################









################SAVE_AND_LOAD_CONTINUUM_POINTS####################

save_cont_ax = plt.axes([0.81, 0.025, 0.1, 0.02])
button2 = Button(save_cont_ax, 'save_cont')
def save_cont_ax(event):
	cont_points = np.transpose(np.array([coordsx, coordsy]))  
	np.savetxt(file3, cont_points)
	print ("continuum points saved")
button2.on_clicked(save_cont_ax)


get_cont_ax = plt.axes([0.81, 0.0035, 0.1, 0.02])
button3 = Button(get_cont_ax, 'get_cont')
def get_cont_ax(event):
	cont_points = np.loadtxt(file3) 
	global coordsx, coordsy
	coordsx = np.append(coordsx, cont_points[:,0])
	coordsy = np.append(coordsy, cont_points[:,1])
	print ("continuum points loaded")
	points = zip(coordsx, coordsy)
	points = sorted(points, key=lambda point: point[0])
	x1, y1 = zip(*points)
	new_length = len(wave)
	l1 = np.searchsorted(wave, min(x1))
	l2 = np.searchsorted(wave, max(x1))
	new_x = np.linspace(min(x1), max(x1), (l2-l1))
	new_y = sp.interpolate.interp1d(x1, y1, kind='cubic')(new_x)
	global cont
	cont = np.zeros([(len(wave))])
	cont[l1:l2] = new_y
	line5.set_data(x1, y1)
	line6.set_data(wave, cont)
	fig.canvas.draw()   
button3.on_clicked(get_cont_ax)

################SAVE_AND_LOAD_CONTINUUM_POINTS####################





################SAVE_AND_LOAD_PARAMETERS####################

save_params_ax = plt.axes([0.905, 0.125, 0.07, 0.02])
button9 = Button(save_params_ax, 'save_params')
def save_params_ax(event):
	params_fitted = []

	if 'logN_HI_fit' in globals():

		print ('Found fitted parameters....')
		print ('Saving Paramters....')

		z_qso = sqso.val   
		resolution = sresol.val
		flux_red = sflux_red.val

		global logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err

		logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err = logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err

		params_fitted = np.append(params_fitted, [logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err, z_qso, resolution, flux_red]) 

	else:

		print ('Could not find fitted parameters....')
		print ('So saving paramters from GUI without error details....')

		logN_HI_fit = slogN_HI.val
		logN_HI_fit_err = 0.0
		logN_H2J0_fit = slogN_H2J0.val
		logN_H2J0_fit_err = 0.0
		logN_H2J1_fit = slogN_H2J1.val
		logN_H2J1_fit_err = 0.0
		logN_H2J2_fit = slogN_H2J2.val
		logN_H2J2_fit_err = 0.0
		logN_H2J3_fit = slogN_H2J3.val
		logN_H2J3_fit_err = 0.0
		logN_H2J4_fit = slogN_H2J4.val
		logN_H2J4_fit_err = 0.0
		logN_H2J5_fit = slogN_H2J5.val
		logN_H2J5_fit_err = 0.0
		logN_H2J6_fit = slogN_H2J6.val
		logN_H2J6_fit_err = 0.0
		logN_H2J7_fit = slogN_H2J7.val
		logN_H2J7_fit_err = 0.0
		logN_HDJ0_fit = slogN_HDJ0.val
		logN_HDJ0_fit_err = 0.0
		logN_HDJ1_fit = slogN_HDJ1.val
		logN_HDJ1_fit_err = 0.0
		logN_HDJ2_fit = slogN_HDJ2.val
		logN_HDJ2_fit_err = 0.0
		flux_red = sflux_red.val
		resolution = sresol.val
		z_HI_fit = sHIdla.val
		z_HI_fit_err = 0.0
		z_H2_fit = sH2dla.val
		z_H2_fit_err = 0.0
		z_qso = sqso.val  

		params_fitted = np.append(params_fitted, [logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err, z_qso, resolution, flux_red])
	

	
	np.savetxt(file9, np.transpose(params_fitted))
	print ("Fitted Parameters Saved")
button9.on_clicked(save_params_ax)




get_params_ax = plt.axes([0.905, 0.105, 0.07, 0.02])
button10 = Button(get_params_ax, 'get_params')

def get_params_ax(event):

	params_fitted_new = np.loadtxt(file9)

	logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err, z_qso, resolution, flux_red = params_fitted_new

	

	slogN_HI.set_val(logN_HI_fit)
	slogN_H2J0.set_val(logN_H2J0_fit)
	slogN_H2J1.set_val(logN_H2J1_fit)
	slogN_H2J2.set_val(logN_H2J2_fit)
	slogN_H2J3.set_val(logN_H2J3_fit)
	slogN_H2J4.set_val(logN_H2J4_fit)
	slogN_H2J5.set_val(logN_H2J5_fit)
	slogN_H2J6.set_val(logN_H2J6_fit)
	slogN_H2J7.set_val(logN_H2J7_fit)
	slogN_HDJ0.set_val(logN_HDJ0_fit)
	slogN_HDJ1.set_val(logN_HDJ1_fit)
	slogN_HDJ2.set_val(logN_HDJ2_fit)
	sHIdla.set_val(z_HI_fit)
	sH2dla.set_val(z_H2_fit)
	sqso.set_val(z_qso)
	sresol.set_val(resolution)
	sflux_red.set_val(flux_red)

	print ("Fitted Parameters Loaded")
button10.on_clicked(get_params_ax)

################SAVE_AND_LOAD_PARAMETERS####################






################EXECUTE_THE_FIT####################

get_fitter_ax = plt.axes([0.905, 0.200, 0.07, 0.02])
button5 = Button(get_fitter_ax, 'Fit')
def get_fitter_ax(event):

	global idx1, idx2, xdata2, ydata2, yerr2, ynew2, contnew2, x, y, err, cont, order_test, coeff
	idx1 = np.searchsorted(xdata, region1_start)
	idx2 = np.searchsorted(xdata, region1_end)
	xdata2 = xdata[idx1:idx2]
	ydata2 = ydata[idx1:idx2]
	yerr2 = flux_err[idx1:idx2]
	ynew2 = ynew[idx1:idx2]
	contnew2 = contdata[idx1:idx2]
	x = xdata2[numpy.logical_not(numpy.isnan(ynew2))]
	y = ydata2[numpy.logical_not(numpy.isnan(ynew2))]
	err = yerr2[numpy.logical_not(numpy.isnan(ynew2))]
	cont = contnew2[numpy.logical_not(numpy.isnan(ynew2))]
	order_test, coeff = chebyshev_order(xdata2, contnew2, stopping_number=1.)
	global cheb_order_HI
	cheb_order_HI = order_test


	global idx3, idx4, xdata3, ydata3, yerr3, ynew3, contnew3, x_H2, y_H2, err_H2, cont_H2, order_test_H2, coeff_H2	
	idx3 = np.searchsorted(xdata, region2_start)
	idx4 = np.searchsorted(xdata, region2_end)
	xdata3 = xdata[idx3:idx4]
	ydata3 = ydata[idx3:idx4]
	yerr3 = flux_err[idx3:idx4]
	ynew3 = ynew[idx3:idx4]
	contnew3 = contdata[idx3:idx4]
	x_H2 = xdata3[numpy.logical_not(numpy.isnan(ynew3))]
	y_H2 = ydata3[numpy.logical_not(numpy.isnan(ynew3))]
	err_H2 = yerr3[numpy.logical_not(numpy.isnan(ynew3))]	
	cont_H2 = contnew3[numpy.logical_not(numpy.isnan(ynew3))]
	order_test_H2, coeff_H2 = chebyshev_order(xdata3, contnew3, stopping_number=1.)
	global cheb_order_H2
	cheb_order_H2 = order_test_H2


	global idx5, idx6, xdata4, ydata4, yerr4, ynew4, contnew4, x_H2_2, y_H2_2, err_H2_2, cont_H2_2, order_test_H2_2, coeff_H2_2
	idx5 = np.searchsorted(xdata, region3_start)
	idx6 = np.searchsorted(xdata, region3_end)	
	xdata4 = xdata[idx5:idx6]
	ydata4 = ydata[idx5:idx6]
	yerr4 = flux_err[idx5:idx6]
	ynew4 = ynew[idx5:idx6]
	contnew4 = contdata[idx5:idx6]
	x_H2_2 = xdata4[numpy.logical_not(numpy.isnan(ynew4))]	
	y_H2_2 = ydata4[numpy.logical_not(numpy.isnan(ynew4))]
	err_H2_2 = yerr4[numpy.logical_not(numpy.isnan(ynew4))]
	cont_H2_2 = contnew4[numpy.logical_not(numpy.isnan(ynew4))]
	order_test_H2_2, coeff_H2_2 = chebyshev_order(xdata4, contnew4, stopping_number=1.)
	global cheb_order_H2_2
	cheb_order_H2_2 = order_test_H2_2


	pars_complete = []
	logN_HI = slogN_HI.val
	logN_H2J0 = slogN_H2J0.val
	logN_H2J1 = slogN_H2J1.val
	logN_H2J2 = slogN_H2J2.val
	logN_H2J3 = slogN_H2J3.val
	logN_H2J4 = slogN_H2J4.val
	logN_H2J5 = slogN_H2J5.val
	logN_H2J6 = slogN_H2J6.val
	logN_H2J7 = slogN_H2J7.val
	logN_HDJ0 = slogN_HDJ0.val
	logN_HDJ1 = slogN_HDJ1.val
	logN_HDJ2 = slogN_HDJ2.val
	flux_red = sflux_red.val
	resolution = sresol.val
	z_HI = sHIdla.val
	z_H2 = sH2dla.val
	z_qso = sqso.val    
	pars_complete = np.append(pars_complete, [logN_HI, logN_H2J0, logN_H2J1, logN_H2J2, logN_H2J3, logN_H2J4, logN_H2J5, logN_H2J6, logN_H2J7, logN_HDJ0, logN_HDJ1, logN_HDJ2, z_HI*1e5, z_H2*1e5])
	pars_complete = np.append(pars_complete, coeff*1e5)
	pars_complete = np.append(pars_complete, coeff_H2*1e5)
	pars_complete = np.append(pars_complete, coeff_H2_2*1e5)
	global res
	res = resolution

	xdata_total = []
	xdata_total = np.append(xdata_total, xdata2)
	xdata_total = np.append(xdata_total, xdata3)
	xdata_total = np.append(xdata_total, xdata4)

	ydata_total = []
	ydata_total = np.append(ydata_total, y)
	ydata_total = np.append(ydata_total, y_H2)
	ydata_total = np.append(ydata_total, y_H2_2)

	errdata_total = []
	errdata_total = np.append(errdata_total, err)
	errdata_total = np.append(errdata_total, err_H2)
	errdata_total = np.append(errdata_total, err_H2_2)


	print ("Fitting.............")
	print ("Press Ctrl+x to stop fitting!")
	method_str = 'trf'
	#method_str = 'dogbox'

	global logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err

	popt, perr = fit_curvefit_complete(pars_complete, xdata_total, ydata_total, errdata_total, fitting_func_complete, method_str)
	cheb_order = cheb_order_HI+cheb_order_H2+cheb_order_H2_2
	cheb_coeff_fit = popt[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):]
	cheb_coeff_HI_fit = cheb_coeff_fit[0:cheb_order_HI]/1e5
	cheb_coeff_H2_fit = cheb_coeff_fit[cheb_order_HI:cheb_order_HI+cheb_order_H2]/1e5
	cheb_coeff_H2_2_fit = cheb_coeff_fit[cheb_order_HI+cheb_order_H2:cheb_order_HI+cheb_order_H2+cheb_order_H2_2]/1e5
	logN_HI_fit, logN_H2J0_fit, logN_H2J1_fit, logN_H2J2_fit, logN_H2J3_fit, logN_H2J4_fit, logN_H2J5_fit, logN_H2J6_fit, logN_H2J7_fit, logN_HDJ0_fit, logN_HDJ1_fit, logN_HDJ2_fit, z_HI_fit, z_H2_fit = (popt[:(len(popt)-cheb_order)])
	z_HI_fit = z_HI_fit/1e5
	z_H2_fit = z_H2_fit/1e5

	cheb_coeff_fit_err = perr[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):]
	cheb_coeff_HI_fit_err = cheb_coeff_fit_err[0:cheb_order_HI]/1e5
	cheb_coeff_H2_fit_err = cheb_coeff_fit_err[cheb_order_HI:cheb_order_HI+cheb_order_H2]/1e5
	cheb_coeff_H2_2_fit_err = cheb_coeff_fit_err[cheb_order_HI+cheb_order_H2:cheb_order_HI+cheb_order_H2+cheb_order_H2_2]/1e5
	logN_HI_fit_err, logN_H2J0_fit_err, logN_H2J1_fit_err, logN_H2J2_fit_err, logN_H2J3_fit_err, logN_H2J4_fit_err, logN_H2J5_fit_err, logN_H2J6_fit_err, logN_H2J7_fit_err, logN_HDJ0_fit_err, logN_HDJ1_fit_err, logN_HDJ2_fit_err, z_HI_fit_err, z_H2_fit_err = (perr[:(len(perr)-cheb_order)])
	z_HI_fit_err = z_HI_fit_err/1e5
	z_H2_fit_err = z_H2_fit_err/1e5

	print ('Coefficients for HI', cheb_coeff_fit)
	print ('Error on coefficients for HI', cheb_coeff_fit_err)
	print ('Coefficients for H2 main region', cheb_coeff_H2_fit)
	print ('Error on coefficients for H2 main region', cheb_coeff_H2_fit_err)
	print ('Coefficients for H2 secondary region', cheb_coeff_H2_2_fit)
	print ('Error on coefficients for H2 secondary region', cheb_coeff_H2_2_fit_err)
	print ('log N(HI) - ', logN_HI_fit, '+/-', logN_HI_fit_err)
	print ('log N(H2-J0) - ', logN_H2J0_fit, '+/-', logN_H2J0_fit_err)
	print ('log N(H2-J1) - ', logN_H2J1_fit, '+/-', logN_H2J1_fit_err)
	print ('log N(H2-J2) - ', logN_H2J2_fit, '+/-', logN_H2J2_fit_err)
	print ('log N(H2-J3) - ', logN_H2J3_fit, '+/-', logN_H2J3_fit_err)
	print ('log N(H2-J4) - ', logN_H2J4_fit, '+/-', logN_H2J4_fit_err)
	print ('log N(H2-J5) - ', logN_H2J5_fit, '+/-', logN_H2J5_fit_err)
	print ('log N(H2-J6) - ', logN_H2J6_fit, '+/-', logN_H2J6_fit_err)
	print ('log N(HD-J7) - ', logN_H2J7_fit, '+/-', logN_H2J7_fit_err)
	print ('log N(HD-J0) - ', logN_HDJ0_fit, '+/-', logN_HDJ0_fit_err)
	print ('log N(HD-J1) - ', logN_HDJ1_fit, '+/-', logN_HDJ1_fit_err)
	print ('log N(HD-J2) - ', logN_HDJ2_fit, '+/-', logN_HDJ2_fit_err)
	print ('N(HI) Redshift - ', z_HI_fit, '+/-', z_HI_fit_err)
	print ('N(H2) Redshift - ', z_H2_fit, '+/-', z_H2_fit_err)


	slogN_HI.set_val(logN_HI_fit)
	slogN_H2J0.set_val(logN_H2J0_fit)
	slogN_H2J1.set_val(logN_H2J1_fit)
	slogN_H2J2.set_val(logN_H2J2_fit)
	slogN_H2J3.set_val(logN_H2J3_fit)
	slogN_H2J4.set_val(logN_H2J4_fit)
	slogN_H2J5.set_val(logN_H2J5_fit)
	slogN_H2J6.set_val(logN_H2J6_fit)
	slogN_H2J7.set_val(logN_H2J7_fit)
	slogN_HDJ0.set_val(logN_HDJ0_fit)
	slogN_HDJ1.set_val(logN_HDJ1_fit)
	slogN_HDJ2.set_val(logN_HDJ2_fit)
	sHIdla.set_val(z_HI_fit)
	sH2dla.set_val(z_H2_fit)

	fitted_profile, fitted_continuum = viewing_func_complete3(wave, *popt)

	line4.set_data(xdata_total, fitted_continuum)

button5.on_clicked(get_fitter_ax)

################EXECUTE_THE_FIT####################





#Button for saving the fit

save_fitted_data_ax = plt.axes([0.905, 0.175, 0.07, 0.02])
button7 = Button(save_fitted_data_ax, 'save_fit')
def save_fitted_data_ax(event):
	x2 = xdata
	y2 = ydata
	err2 = ydata_err
	cont2 = contdata
	fit2 = profile_updated
	data2 = np.transpose(np.array([x2, y2, err2, cont2, fit2]))


	x3_selected = xdata[numpy.logical_not(numpy.isnan(ynew))]
	y3_selected = ydata[numpy.logical_not(numpy.isnan(ynew))]
	err3_selected = ydata_err[numpy.logical_not(numpy.isnan(ynew))]
	cont3_selected = contdata[numpy.logical_not(numpy.isnan(ynew))]
	data3_selected = np.transpose(np.array([x3_selected, y3_selected, err3_selected, cont3_selected]))
	

	if (os.path.isfile(file7)):
		save_prompt2=raw_input("File already exists. Overwrite?(y/n) : ")
        	if (save_prompt2=='y'):
        		np.savetxt(file7, data2)
			print ('saved Main File')
        		np.savetxt(file8, data3_selected)
			print ('saved Selected Datapoint File')
			
        		print ('saved')
        	else:
        		print ('data file not saved')
	else:
       		np.savetxt(file7, data2)
       		np.savetxt(file8, data2)
       		print ('saved')
button7.on_clicked(save_fitted_data_ax)

















###########################################PLOTTING_AND_ANALYSING_RESULTS###########################################


#Button for plotting the fit of HI and H2. The function plots the fitted data and saves it in a file. Please check out the plotter file 'plotting_hydrogen2.py' for details

plot_fitted_data_ax = plt.axes([0.005, 0.650, 0.07, 0.02])
button8 = Button(plot_fitted_data_ax, 'plot_total_fit')

def plot_fitted_data_ax(event):
	print ('Plotting....')
	z_H2 = sH2dla.val
	str_prog_code = 'python ' + 'plotting_files/plotting_hydrogen2.py ' + str(file7) + ' ' + str(file8) + ' ' + str(float(z_H2))
	print (str_prog_code)
	os.system(str_prog_code)
	print ('Plotted')


button8.on_clicked(plot_fitted_data_ax)








#Button for plotting the fit of H2 in velocity space. The function plots the fitted data and saves it in a file. Please check out the plotter file 'plotting_H2J0_H2J1_like_metals.py' for details

plot_fitted_H2J0_data_ax = plt.axes([0.005, 0.625, 0.07, 0.02])
button11 = Button(plot_fitted_H2J0_data_ax, 'plot_H2J0/H2J1')

def plot_fitted_H2J0_data_ax(event):
	print ('Plotting....')
	z_HI = sHIdla.val
	z_H2 = sH2dla.val
	str_prog_code = 'python ' + 'plotting_files/plotting_H2J0_H2J1_like_metals.py ' + str(file7) + ' ' + str(file8) + ' ' + str(float(z_HI)) + ' ' + str(float(z_H2)) + ' ' + str('subplot_grid_custom')
	print (str_prog_code)
	os.system(str_prog_code)
	print ('Plotted')


button11.on_clicked(plot_fitted_H2J0_data_ax)






#Button for plotting the fit of HI in velocity space. The function plots the fitted data and saves it in a file. Please check out the plotter file 'plotting_HI.py' for details

plot_fitted_HI_data_ax = plt.axes([0.005, 0.600, 0.07, 0.02])
button12 = Button(plot_fitted_HI_data_ax, 'plot_HI')

def plot_fitted_HI_data_ax(event):
	print ('Plotting....')
	z_HI = sHIdla.val
	z_H2 = sH2dla.val
	str_prog_code = 'python ' + 'plotting_files/plotting_HI.py ' + str(file7) + ' ' + str(file8) + ' ' + str(float(z_HI)) + ' ' + str(float(z_H2)) + ' ' + str('subplot_grid_custom')
	print (str_prog_code)
	os.system(str_prog_code)
	print ('Plotted')


button12.on_clicked(plot_fitted_HI_data_ax)










#Button for making an excel table with the fitted information.

save_binary_table_ax = plt.axes([0.005, 0.525, 0.07, 0.02])
button15 = Button(save_binary_table_ax, 'Save Binary Table')

def save_binary_table_ax(event):

	print ("Saving Binary Table in Excel")
	print ("##################################")
	print ("\n" "\n")


	params_fitted = []

	if 'logN_HI_fit' in globals():

		print ('Found fitted parameters....')
		print ('Printing Paramters in latex table form....')
		print ("\n" "\n")

		z_qso = sqso.val   
		resolution = sresol.val
		flux_red = sflux_red.val

		global logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err

		logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err = logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err


		lin_N_H2J0_fit = 10**logN_H2J0_fit
		lin_N_H2J1_fit = 10**logN_H2J1_fit
		lin_N_H2J2_fit = 10**logN_H2J2_fit
		lin_N_H2J3_fit = 10**logN_H2J3_fit
		lin_N_H2J4_fit = 10**logN_H2J4_fit
		lin_N_H2J5_fit = 10**logN_H2J5_fit
		lin_N_H2J6_fit = 10**logN_H2J6_fit
		lin_N_H2J7_fit = 10**logN_H2J7_fit
		
		linN_H2J0_fit_err = log_to_real_err(logN_H2J0_fit, logN_H2J0_fit_err)
		linN_H2J1_fit_err = log_to_real_err(logN_H2J1_fit, logN_H2J1_fit_err)
		linN_H2J2_fit_err = log_to_real_err(logN_H2J2_fit, logN_H2J2_fit_err)
		linN_H2J3_fit_err = log_to_real_err(logN_H2J3_fit, logN_H2J3_fit_err)
		linN_H2J4_fit_err = log_to_real_err(logN_H2J4_fit, logN_H2J4_fit_err)
		linN_H2J5_fit_err = log_to_real_err(logN_H2J5_fit, logN_H2J5_fit_err)
		linN_H2J6_fit_err = log_to_real_err(logN_H2J6_fit, logN_H2J6_fit_err)
		linN_H2J7_fit_err = log_to_real_err(logN_H2J7_fit, logN_H2J7_fit_err)

		total_lin_H2_col_den = lin_N_H2J0_fit + lin_N_H2J1_fit + lin_N_H2J2_fit +lin_N_H2J3_fit + lin_N_H2J4_fit + lin_N_H2J5_fit + lin_N_H2J6_fit + lin_N_H2J7_fit
		total_lin_H2_col_den_err = np.sqrt(linN_H2J0_fit_err**2 + linN_H2J1_fit_err**2 + linN_H2J2_fit_err**2 +linN_H2J3_fit_err**2 + linN_H2J4_fit_err**2 + linN_H2J5_fit_err**2 + linN_H2J6_fit_err**2 + linN_H2J7_fit_err**2)
		total_log_H2_col_den = np.log10(total_lin_H2_col_den)
		total_log_H2_col_den_err = real_to_log_err(total_lin_H2_col_den, total_lin_H2_col_den_err)


		params_fitted = np.append(params_fitted, [logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, total_log_H2_col_den, total_log_H2_col_den_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err, z_qso, resolution, flux_red]) 



		



	else:

		print ('Could not find fitted parameters....')
		print ('So Printing Paramters in latex table form from GUI without error details....')
		print ("\n" "\n")
	

		logN_HI_fit = slogN_HI.val
		logN_HI_fit_err = 0.0
		logN_H2J0_fit = slogN_H2J0.val
		logN_H2J0_fit_err = 0.0
		logN_H2J1_fit = slogN_H2J1.val
		logN_H2J1_fit_err = 0.0
		logN_H2J2_fit = slogN_H2J2.val
		logN_H2J2_fit_err = 0.0
		logN_H2J3_fit = slogN_H2J3.val
		logN_H2J3_fit_err = 0.0
		logN_H2J4_fit = slogN_H2J4.val
		logN_H2J4_fit_err = 0.0
		logN_H2J5_fit = slogN_H2J5.val
		logN_H2J5_fit_err = 0.0
		logN_H2J6_fit = slogN_H2J6.val
		logN_H2J6_fit_err = 0.0
		logN_H2J7_fit = slogN_H2J7.val
		logN_H2J7_fit_err = 0.0
		logN_HDJ0_fit = slogN_HDJ0.val
		logN_HDJ0_fit_err = 0.0
		logN_HDJ1_fit = slogN_HDJ1.val
		logN_HDJ1_fit_err = 0.0
		logN_HDJ2_fit = slogN_HDJ2.val
		logN_HDJ2_fit_err = 0.0
		flux_red = sflux_red.val
		resolution = sresol.val
		z_HI_fit = sHIdla.val
		z_HI_fit_err = 0.0
		z_H2_fit = sH2dla.val
		z_H2_fit_err = 0.0
		z_qso = sqso.val  



		lin_N_H2J0_fit = 10**logN_H2J0_fit
		lin_N_H2J1_fit = 10**logN_H2J1_fit
		lin_N_H2J2_fit = 10**logN_H2J2_fit
		lin_N_H2J3_fit = 10**logN_H2J3_fit
		lin_N_H2J4_fit = 10**logN_H2J4_fit
		lin_N_H2J5_fit = 10**logN_H2J5_fit
		lin_N_H2J6_fit = 10**logN_H2J6_fit
		lin_N_H2J7_fit = 10**logN_H2J7_fit
		
		linN_H2J0_fit_err = log_to_real_err(logN_H2J0_fit, logN_H2J0_fit_err)
		linN_H2J1_fit_err = log_to_real_err(logN_H2J1_fit, logN_H2J1_fit_err)
		linN_H2J2_fit_err = log_to_real_err(logN_H2J2_fit, logN_H2J2_fit_err)
		linN_H2J3_fit_err = log_to_real_err(logN_H2J3_fit, logN_H2J3_fit_err)
		linN_H2J4_fit_err = log_to_real_err(logN_H2J4_fit, logN_H2J4_fit_err)
		linN_H2J5_fit_err = log_to_real_err(logN_H2J5_fit, logN_H2J5_fit_err)
		linN_H2J6_fit_err = log_to_real_err(logN_H2J6_fit, logN_H2J6_fit_err)
		linN_H2J7_fit_err = log_to_real_err(logN_H2J7_fit, logN_H2J7_fit_err)

		total_lin_H2_col_den = lin_N_H2J0_fit + lin_N_H2J1_fit + lin_N_H2J2_fit +lin_N_H2J3_fit + lin_N_H2J4_fit + lin_N_H2J5_fit + lin_N_H2J6_fit + lin_N_H2J7_fit
		total_lin_H2_col_den_err = np.sqrt(linN_H2J0_fit_err**2 + linN_H2J1_fit_err**2 + linN_H2J2_fit_err**2 +linN_H2J3_fit_err**2 + linN_H2J4_fit_err**2 + linN_H2J5_fit_err**2 + linN_H2J6_fit_err**2 + linN_H2J7_fit_err**2)
		total_log_H2_col_den = np.log10(total_lin_H2_col_den)
		total_log_H2_col_den_err = real_to_log_err(total_lin_H2_col_den, total_lin_H2_col_den_err)



		params_fitted = np.append(params_fitted, [logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, total_log_H2_col_den, total_log_H2_col_den_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err, z_qso, resolution, flux_red])


	params_fitted_name = np.chararray([], itemsize=35)
	params_fitted_name = np.delete(params_fitted_name, 0)
	params_fitted_name = np.append(params_fitted_name, ['log(N(HI))', 'log(N(HI)) err', 'log(N(H2J0))', 'log(N(H2J0)) err', 'log(N(H2J1))', 'log(N(H2J1)) err', 'log(N(H2J2))', 'log(N(H2J2)) err', 'log(N(H2J3))', 'log(N(H2J3)) err', 'log(N(H2J4))', 'log(N(H2J4)) err', 'log(N(H2J5))', 'log(N(H2J5)) err', 'log(N(H2J6))', 'log(N(H2J6)) err', 'log(N(H2J7))', 'log(N(H2J7)) err', 'Total(H2)', 'Total(H2 err)', 'log(N(HDJ0))', 'log(N(HDJ0)) err', 'log(N(HDJ1))', 'log(N(HDJ1)) err', 'log(N(HDJ2))', 'log(N(HDJ2)) err', 'HI Redshift', 'HI Redshift err', 'H2 Redshift', 'H2 Redshift err', 'QSO Redshift', 'Resolution', 'Flux Reduction index'])
		



	table_array = np.chararray([len(params_fitted)], itemsize=35)
	table_array = params_fitted
	np.savetxt(str(file10), table_array, fmt='%s')
	
	print ("\n" "\n")
	print ("##################################")
	print ("Parameter Table for Latex Created")


button15.on_clicked(save_binary_table_ax)








#Button for making a latex table with the fitted information. The function prints the script for latex that can be directly used.

make_latex_table_ax = plt.axes([0.005, 0.575, 0.07, 0.02])
button13 = Button(make_latex_table_ax, 'latex_table')

def make_latex_table_ax(event):

	print ("Creating Parameter Table for Latex")
	print ("##################################")
	print ("\n" "\n")


	params_fitted = []

	if 'logN_HI_fit' in globals():

		print ('Found fitted parameters....')
		print ('Printing Paramters in latex table form....')
		print ("\n" "\n")

		z_qso = sqso.val   
		resolution = sresol.val
		flux_red = sflux_red.val

		global logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err

		logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err = logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err


		lin_N_H2J0_fit = 10**logN_H2J0_fit
		lin_N_H2J1_fit = 10**logN_H2J1_fit
		lin_N_H2J2_fit = 10**logN_H2J2_fit
		lin_N_H2J3_fit = 10**logN_H2J3_fit
		lin_N_H2J4_fit = 10**logN_H2J4_fit
		lin_N_H2J5_fit = 10**logN_H2J5_fit
		lin_N_H2J6_fit = 10**logN_H2J6_fit
		lin_N_H2J7_fit = 10**logN_H2J7_fit
		
		linN_H2J0_fit_err = log_to_real_err(logN_H2J0_fit, logN_H2J0_fit_err)
		linN_H2J1_fit_err = log_to_real_err(logN_H2J1_fit, logN_H2J1_fit_err)
		linN_H2J2_fit_err = log_to_real_err(logN_H2J2_fit, logN_H2J2_fit_err)
		linN_H2J3_fit_err = log_to_real_err(logN_H2J3_fit, logN_H2J3_fit_err)
		linN_H2J4_fit_err = log_to_real_err(logN_H2J4_fit, logN_H2J4_fit_err)
		linN_H2J5_fit_err = log_to_real_err(logN_H2J5_fit, logN_H2J5_fit_err)
		linN_H2J6_fit_err = log_to_real_err(logN_H2J6_fit, logN_H2J6_fit_err)
		linN_H2J7_fit_err = log_to_real_err(logN_H2J7_fit, logN_H2J7_fit_err)

		total_lin_H2_col_den = lin_N_H2J0_fit + lin_N_H2J1_fit + lin_N_H2J2_fit +lin_N_H2J3_fit + lin_N_H2J4_fit + lin_N_H2J5_fit + lin_N_H2J6_fit + lin_N_H2J7_fit
		total_lin_H2_col_den_err = np.sqrt(linN_H2J0_fit_err**2 + linN_H2J1_fit_err**2 + linN_H2J2_fit_err**2 +linN_H2J3_fit_err**2 + linN_H2J4_fit_err**2 + linN_H2J5_fit_err**2 + linN_H2J6_fit_err**2 + linN_H2J7_fit_err**2)
		total_log_H2_col_den = np.log10(total_lin_H2_col_den)
		total_log_H2_col_den_err = real_to_log_err(total_lin_H2_col_den, total_lin_H2_col_den_err)


		params_fitted = np.append(params_fitted, [logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, total_log_H2_col_den, total_log_H2_col_den_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err, z_qso, resolution, flux_red]) 



		



	else:

		print ('Could not find fitted parameters....')
		print ('So Printing Paramters in latex table form from GUI without error details....')
		print ("\n" "\n")
	

		logN_HI_fit = slogN_HI.val
		logN_HI_fit_err = 0.0
		logN_H2J0_fit = slogN_H2J0.val
		logN_H2J0_fit_err = 0.0
		logN_H2J1_fit = slogN_H2J1.val
		logN_H2J1_fit_err = 0.0
		logN_H2J2_fit = slogN_H2J2.val
		logN_H2J2_fit_err = 0.0
		logN_H2J3_fit = slogN_H2J3.val
		logN_H2J3_fit_err = 0.0
		logN_H2J4_fit = slogN_H2J4.val
		logN_H2J4_fit_err = 0.0
		logN_H2J5_fit = slogN_H2J5.val
		logN_H2J5_fit_err = 0.0
		logN_H2J6_fit = slogN_H2J6.val
		logN_H2J6_fit_err = 0.0
		logN_H2J7_fit = slogN_H2J7.val
		logN_H2J7_fit_err = 0.0
		logN_HDJ0_fit = slogN_HDJ0.val
		logN_HDJ0_fit_err = 0.0
		logN_HDJ1_fit = slogN_HDJ1.val
		logN_HDJ1_fit_err = 0.0
		logN_HDJ2_fit = slogN_HDJ2.val
		logN_HDJ2_fit_err = 0.0
		flux_red = sflux_red.val
		resolution = sresol.val
		z_HI_fit = sHIdla.val
		z_HI_fit_err = 0.0
		z_H2_fit = sH2dla.val
		z_H2_fit_err = 0.0
		z_qso = sqso.val  



		lin_N_H2J0_fit = 10**logN_H2J0_fit
		lin_N_H2J1_fit = 10**logN_H2J1_fit
		lin_N_H2J2_fit = 10**logN_H2J2_fit
		lin_N_H2J3_fit = 10**logN_H2J3_fit
		lin_N_H2J4_fit = 10**logN_H2J4_fit
		lin_N_H2J5_fit = 10**logN_H2J5_fit
		lin_N_H2J6_fit = 10**logN_H2J6_fit
		lin_N_H2J7_fit = 10**logN_H2J7_fit
		
		linN_H2J0_fit_err = log_to_real_err(logN_H2J0_fit, logN_H2J0_fit_err)
		linN_H2J1_fit_err = log_to_real_err(logN_H2J1_fit, logN_H2J1_fit_err)
		linN_H2J2_fit_err = log_to_real_err(logN_H2J2_fit, logN_H2J2_fit_err)
		linN_H2J3_fit_err = log_to_real_err(logN_H2J3_fit, logN_H2J3_fit_err)
		linN_H2J4_fit_err = log_to_real_err(logN_H2J4_fit, logN_H2J4_fit_err)
		linN_H2J5_fit_err = log_to_real_err(logN_H2J5_fit, logN_H2J5_fit_err)
		linN_H2J6_fit_err = log_to_real_err(logN_H2J6_fit, logN_H2J6_fit_err)
		linN_H2J7_fit_err = log_to_real_err(logN_H2J7_fit, logN_H2J7_fit_err)

		total_lin_H2_col_den = lin_N_H2J0_fit + lin_N_H2J1_fit + lin_N_H2J2_fit +lin_N_H2J3_fit + lin_N_H2J4_fit + lin_N_H2J5_fit + lin_N_H2J6_fit + lin_N_H2J7_fit
		total_lin_H2_col_den_err = np.sqrt(linN_H2J0_fit_err**2 + linN_H2J1_fit_err**2 + linN_H2J2_fit_err**2 +linN_H2J3_fit_err**2 + linN_H2J4_fit_err**2 + linN_H2J5_fit_err**2 + linN_H2J6_fit_err**2 + linN_H2J7_fit_err**2)
		total_log_H2_col_den = np.log10(total_lin_H2_col_den)
		total_log_H2_col_den_err = real_to_log_err(total_lin_H2_col_den, total_lin_H2_col_den_err)



		params_fitted = np.append(params_fitted, [logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, total_log_H2_col_den, total_log_H2_col_den_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err, z_qso, resolution, flux_red])


	params_fitted_name = np.chararray([], itemsize=35)
	params_fitted_name = np.delete(params_fitted_name, 0)
	params_fitted_name = np.append(params_fitted_name, ['log(N(HI))', 'log(N(HI)) err', 'log(N(H2J0))', 'log(N(H2J0)) err', 'log(N(H2J1))', 'log(N(H2J1)) err', 'log(N(H2J2))', 'log(N(H2J2)) err', 'log(N(H2J3))', 'log(N(H2J3)) err', 'log(N(H2J4))', 'log(N(H2J4)) err', 'log(N(H2J5))', 'log(N(H2J5)) err', 'log(N(H2J6))', 'log(N(H2J6)) err', 'log(N(H2J7))', 'log(N(H2J7)) err', 'Total(H2)', 'Total(H2 err)', 'log(N(HDJ0))', 'log(N(HDJ0)) err', 'log(N(HDJ1))', 'log(N(HDJ1)) err', 'log(N(HDJ2))', 'log(N(HDJ2)) err', 'HI Redshift', 'HI Redshift err', 'H2 Redshift', 'H2 Redshift err', 'QSO Redshift', 'Resolution', 'Flux Reduction index'])
		


	table_array = np.chararray([len(params_fitted), 2], itemsize=35)
	table_array[:,0] = params_fitted_name
	table_array[:,1] = params_fitted

	print(tabulate(table_array, tablefmt="latex", floatfmt="2.6f"))

	
	print ("\n" "\n")
	print ("##################################")
	print ("Parameter Table for Latex Created")


button13.on_clicked(make_latex_table_ax)





#Button for making an H2 excitation diagram comparing with a DLA, a star and a MW ISM cloud. It will also fit the lower rotational levels, J=0,1 and 2 to give a gas temperature. The resultant figure will be saved in a file. Please check out the code - excitation_diagram2.py for further details. 

plot_excitation_diagram_ax = plt.axes([0.005, 0.550, 0.07, 0.02])
button14 = Button(plot_excitation_diagram_ax, 'excitation')

def plot_excitation_diagram_ax(event):

	print ('Plotting excitation diagram....')

	params_fitted = []

	if 'logN_HI_fit' in globals():

		print ('Found fitted parameters....')

		z_qso = sqso.val   
		resolution = sresol.val
		flux_red = sflux_red.val

		global logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err

		logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err = logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err

		params_fitted = np.append(params_fitted, [logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err, z_qso, resolution, flux_red]) 

	else:

		print ('Could not find fitted parameters....')
		print ('So Printing Paramters in latex table form from GUI without error details....')
		print ("\n" "\n")

		logN_HI_fit = slogN_HI.val
		logN_HI_fit_err = 0.0
		logN_H2J0_fit = slogN_H2J0.val
		logN_H2J0_fit_err = 0.0
		logN_H2J1_fit = slogN_H2J1.val
		logN_H2J1_fit_err = 0.0
		logN_H2J2_fit = slogN_H2J2.val
		logN_H2J2_fit_err = 0.0
		logN_H2J3_fit = slogN_H2J3.val
		logN_H2J3_fit_err = 0.0
		logN_H2J4_fit = slogN_H2J4.val
		logN_H2J4_fit_err = 0.0
		logN_H2J5_fit = slogN_H2J5.val
		logN_H2J5_fit_err = 0.0
		logN_H2J6_fit = slogN_H2J6.val
		logN_H2J6_fit_err = 0.0
		logN_H2J7_fit = slogN_H2J7.val
		logN_H2J7_fit_err = 0.0
		logN_HDJ0_fit = slogN_HDJ0.val
		logN_HDJ0_fit_err = 0.0
		logN_HDJ1_fit = slogN_HDJ1.val
		logN_HDJ1_fit_err = 0.0
		logN_HDJ2_fit = slogN_HDJ2.val
		logN_HDJ2_fit_err = 0.0
		flux_red = sflux_red.val
		resolution = sresol.val
		z_HI_fit = sHIdla.val
		z_HI_fit_err = 0.0
		z_H2_fit = sH2dla.val
		z_H2_fit_err = 0.0
		z_qso = sqso.val  

		params_fitted = np.append(params_fitted, [logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err, z_qso, resolution, flux_red])

	

	file_name_excitation = file_name1[:-9] + '_excitation_diagram' + '.pdf'

	str_prog_code = 'python ' + 'h2_excitation_diagram_files/excitation_diagram2.py' + ' ' + str(logN_H2J0_fit) + ' ' + str(logN_H2J0_fit_err) + ' ' + str(logN_H2J1_fit) + ' ' + str(logN_H2J1_fit_err) + ' ' + str(logN_H2J2_fit) + ' ' + str(logN_H2J2_fit_err) + ' ' + str(logN_H2J3_fit) + ' ' + str(logN_H2J3_fit_err) + ' ' + str(logN_H2J4_fit) + ' ' + str(logN_H2J4_fit_err) + ' ' + str(logN_H2J5_fit) + ' ' + str(logN_H2J5_fit_err) + ' ' + str(logN_H2J6_fit) + ' ' + str(logN_H2J6_fit_err) + ' ' + str(logN_H2J7_fit) + ' ' + str(logN_H2J7_fit_err) + ' ' + str(file_name_excitation)
	print (str_prog_code)
	os.system(str_prog_code)


	print ('Plotted excitation diagram....')


button14.on_clicked(plot_excitation_diagram_ax)


###########################################PLOTTING_AND_ANALYSING_RESULTS###########################################

























#This function can be used to find the b-value of H2 in cases where the resolution does not allow identifying the proper b-value. The function fits for column density while varying the Doppler parameter from 1km/s to 10km/s. However this range can easily be changed from inside this function.  

fit_b_val_ax = plt.axes([0.005, 0.425, 0.07, 0.02])
button16 = Button(fit_b_val_ax, 'fit_b_val')

def fit_b_val_ax(event):
	print ('Test1....')

	
	global idx1, idx2, xdata2, ydata2, yerr2, ynew2, contnew2, x, y, err, cont, order_test, coeff
	idx1 = np.searchsorted(xdata, region1_start)
	idx2 = np.searchsorted(xdata, region1_end)
	xdata2 = xdata[idx1:idx2]
	ydata2 = ydata[idx1:idx2]
	yerr2 = flux_err[idx1:idx2]
	ynew2 = ynew[idx1:idx2]
	contnew2 = contdata[idx1:idx2]
	x = xdata2[numpy.logical_not(numpy.isnan(ynew2))]
	y = ydata2[numpy.logical_not(numpy.isnan(ynew2))]
	err = yerr2[numpy.logical_not(numpy.isnan(ynew2))]
	cont = contnew2[numpy.logical_not(numpy.isnan(ynew2))]
	order_test, coeff = chebyshev_order(xdata2, contnew2, stopping_number=1.)
	global cheb_order_HI
	cheb_order_HI = order_test


	global idx3, idx4, xdata3, ydata3, yerr3, ynew3, contnew3, x_H2, y_H2, err_H2, cont_H2, order_test_H2, coeff_H2	
	idx3 = np.searchsorted(xdata, region2_start)
	idx4 = np.searchsorted(xdata, region2_end)
	xdata3 = xdata[idx3:idx4]
	ydata3 = ydata[idx3:idx4]
	yerr3 = flux_err[idx3:idx4]
	ynew3 = ynew[idx3:idx4]
	contnew3 = contdata[idx3:idx4]
	x_H2 = xdata3[numpy.logical_not(numpy.isnan(ynew3))]
	y_H2 = ydata3[numpy.logical_not(numpy.isnan(ynew3))]
	err_H2 = yerr3[numpy.logical_not(numpy.isnan(ynew3))]	
	cont_H2 = contnew3[numpy.logical_not(numpy.isnan(ynew3))]
	order_test_H2, coeff_H2 = chebyshev_order(xdata3, contnew3, stopping_number=1.)
	global cheb_order_H2
	cheb_order_H2 = order_test_H2


	global idx5, idx6, xdata4, ydata4, yerr4, ynew4, contnew4, x_H2_2, y_H2_2, err_H2_2, cont_H2_2, order_test_H2_2, coeff_H2_2
	idx5 = np.searchsorted(xdata, region3_start)
	idx6 = np.searchsorted(xdata, region3_end)	
	xdata4 = xdata[idx5:idx6]
	ydata4 = ydata[idx5:idx6]
	yerr4 = flux_err[idx5:idx6]
	ynew4 = ynew[idx5:idx6]
	contnew4 = contdata[idx5:idx6]
	x_H2_2 = xdata4[numpy.logical_not(numpy.isnan(ynew4))]	
	y_H2_2 = ydata4[numpy.logical_not(numpy.isnan(ynew4))]
	err_H2_2 = yerr4[numpy.logical_not(numpy.isnan(ynew4))]
	cont_H2_2 = contnew4[numpy.logical_not(numpy.isnan(ynew4))]
	order_test_H2_2, coeff_H2_2 = chebyshev_order(xdata4, contnew4, stopping_number=1.)
	global cheb_order_H2_2
	cheb_order_H2_2 = order_test_H2_2


	pars_complete = []
	logN_HI = slogN_HI.val
	logN_H2J0 = slogN_H2J0.val
	logN_H2J1 = slogN_H2J1.val
	logN_H2J2 = slogN_H2J2.val
	logN_H2J3 = slogN_H2J3.val
	logN_H2J4 = slogN_H2J4.val
	logN_H2J5 = slogN_H2J5.val
	logN_H2J6 = slogN_H2J6.val
	logN_H2J7 = slogN_H2J7.val
	logN_HDJ0 = slogN_HDJ0.val
	logN_HDJ1 = slogN_HDJ1.val
	logN_HDJ2 = slogN_HDJ2.val
	flux_red = sflux_red.val
	resolution = sresol.val
	z_HI = sHIdla.val
	z_H2 = sH2dla.val
	z_qso = sqso.val    
	pars_complete = np.append(pars_complete, [logN_HI, logN_H2J0, logN_H2J1, logN_H2J2, logN_H2J3, logN_H2J4, logN_H2J5, logN_H2J6, logN_H2J7, logN_HDJ0, logN_HDJ1, logN_HDJ2, z_HI*1e5, z_H2*1e5])


	b_val_new = np.array([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.])

	global res
	res = resolution


	global logN_HI_fit, logN_HI_fit_err, logN_H2J0_fit, logN_H2J0_fit_err, logN_H2J1_fit, logN_H2J1_fit_err, logN_H2J2_fit, logN_H2J2_fit_err, logN_H2J3_fit, logN_H2J3_fit_err, logN_H2J4_fit, logN_H2J4_fit_err, logN_H2J5_fit, logN_H2J5_fit_err, logN_H2J6_fit, logN_H2J6_fit_err, logN_H2J7_fit, logN_H2J7_fit_err, logN_HDJ0_fit, logN_HDJ0_fit_err, logN_HDJ1_fit, logN_HDJ1_fit_err, logN_HDJ2_fit, logN_HDJ2_fit_err, z_HI_fit, z_HI_fit_err, z_H2_fit, z_H2_fit_err


	logN_HI_fit_array = []
	logN_HI_fit_err_array = [] 
	logN_H2J0_fit_array = []
	logN_H2J0_fit_err_array = [] 
	logN_H2J1_fit_array = []
	logN_H2J1_fit_err_array = []
	logN_H2J2_fit_array = []
	logN_H2J2_fit_err_array = []
	logN_H2J3_fit_array = [] 
	logN_H2J3_fit_err_array = []
	logN_H2J4_fit_array = []
	logN_H2J4_fit_err_array = []
	logN_H2J5_fit_array = [] 
	logN_H2J5_fit_err_array = [] 
	logN_H2J6_fit_array = [] 
	logN_H2J6_fit_err_array = [] 
	logN_H2J7_fit_array = [] 
	logN_H2J7_fit_err_array = [] 
	logN_HDJ0_fit_array = [] 
	logN_HDJ0_fit_err_array = [] 
	logN_HDJ1_fit_array = [] 
	logN_HDJ1_fit_err_array = [] 
	logN_HDJ2_fit_array = [] 
	logN_HDJ2_fit_err_array = [] 
	b_val_array = []
	b_val_array_err = []
	

	for i in range(len(b_val_new)):
		pars_complete_new = pars_complete
		pars_complete_new = np.append(pars_complete_new, b_val_new[i])
		pars_complete_new = np.append(pars_complete_new, coeff*1e5)
		pars_complete_new = np.append(pars_complete_new, coeff_H2*1e5)
		pars_complete_new = np.append(pars_complete_new, coeff_H2_2*1e5)
		xdata_total = []
		xdata_total = np.append(xdata_total, xdata2)
		xdata_total = np.append(xdata_total, xdata3)
		xdata_total = np.append(xdata_total, xdata4)
		ydata_total = []
		ydata_total = np.append(ydata_total, y)
		ydata_total = np.append(ydata_total, y_H2)
		ydata_total = np.append(ydata_total, y_H2_2)
		errdata_total = []
		errdata_total = np.append(errdata_total, err)
		errdata_total = np.append(errdata_total, err_H2)
		errdata_total = np.append(errdata_total, err_H2_2)
		print ("Fitting.............")
		print ("Press Ctrl+x to stop fitting!")
		method_str = 'trf'
		#method_str = 'dogbox'
		popt, perr = fit_curvefit_complete_new(pars_complete_new, xdata_total, ydata_total, errdata_total, fitting_func_complete_new, method_str)
		cheb_order = cheb_order_HI+cheb_order_H2+cheb_order_H2_2
		cheb_coeff_fit = popt[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):]
		cheb_coeff_HI_fit = cheb_coeff_fit[0:cheb_order_HI]/1e5
		cheb_coeff_H2_fit = cheb_coeff_fit[cheb_order_HI:cheb_order_HI+cheb_order_H2]/1e5
		cheb_coeff_H2_2_fit = cheb_coeff_fit[cheb_order_HI+cheb_order_H2:cheb_order_HI+cheb_order_H2+cheb_order_H2_2]/1e5
		logN_HI_fit, logN_H2J0_fit, logN_H2J1_fit, logN_H2J2_fit, logN_H2J3_fit, logN_H2J4_fit, logN_H2J5_fit, logN_H2J6_fit, logN_H2J7_fit, logN_HDJ0_fit, logN_HDJ1_fit, logN_HDJ2_fit, z_HI_fit, z_H2_fit, b_val_fitted_new = (popt[:(len(popt)-cheb_order)])
		z_HI_fit = z_HI_fit/1e5
		z_H2_fit = z_H2_fit/1e5
		cheb_coeff_fit_err = perr[-(cheb_order_HI+cheb_order_H2+cheb_order_H2_2):]
		cheb_coeff_HI_fit_err = cheb_coeff_fit_err[0:cheb_order_HI]/1e5
		cheb_coeff_H2_fit_err = cheb_coeff_fit_err[cheb_order_HI:cheb_order_HI+cheb_order_H2]/1e5
		cheb_coeff_H2_2_fit_err = cheb_coeff_fit_err[cheb_order_HI+cheb_order_H2:cheb_order_HI+cheb_order_H2+cheb_order_H2_2]/1e5
		logN_HI_fit_err, logN_H2J0_fit_err, logN_H2J1_fit_err, logN_H2J2_fit_err, logN_H2J3_fit_err, logN_H2J4_fit_err, logN_H2J5_fit_err, logN_H2J6_fit_err, logN_H2J7_fit_err, logN_HDJ0_fit_err, logN_HDJ1_fit_err, logN_HDJ2_fit_err, z_HI_fit_err, z_H2_fit_err, b_val_fitted_new_err = (perr[:(len(perr)-cheb_order)])
		z_HI_fit_err = z_HI_fit_err/1e5
		z_H2_fit_err = z_H2_fit_err/1e5
		print (b_val_new[i])

		#print (logN_HI_fit)
		#print (logN_HI_fit_err)

		logN_HI_fit_array = np.append(logN_HI_fit_array, 10**logN_HI_fit)
		logN_HI_fit_err_array = np.append(logN_HI_fit_err_array, log_to_real_err(logN_HI_fit, logN_HI_fit_err)) 
		logN_H2J0_fit_array = np.append(logN_H2J0_fit_array, 10**logN_H2J0_fit)
		logN_H2J0_fit_err_array = np.append(logN_H2J0_fit_err_array, log_to_real_err(logN_H2J0_fit, logN_H2J0_fit_err)) 
		logN_H2J1_fit_array = np.append(logN_H2J1_fit_array, 10**logN_H2J1_fit)
		logN_H2J1_fit_err_array = np.append(logN_H2J1_fit_err_array, log_to_real_err(logN_H2J1_fit, logN_H2J1_fit_err))
		logN_H2J2_fit_array = np.append(logN_H2J2_fit_array, 10**logN_H2J2_fit)
		logN_H2J2_fit_err_array = np.append(logN_H2J2_fit_err_array, log_to_real_err(logN_H2J2_fit, logN_H2J2_fit_err))
		logN_H2J3_fit_array = np.append(logN_H2J3_fit_array, 10**logN_H2J3_fit) 
		logN_H2J3_fit_err_array = np.append(logN_H2J3_fit_err_array, log_to_real_err(logN_H2J3_fit, logN_H2J3_fit_err))
		logN_H2J4_fit_array = np.append(logN_H2J4_fit_array, 10**logN_H2J4_fit)
		logN_H2J4_fit_err_array = np.append(logN_H2J4_fit_err_array, log_to_real_err(logN_H2J4_fit, logN_H2J4_fit_err))
		logN_H2J5_fit_array = np.append(logN_H2J5_fit_array, 10**logN_H2J5_fit) 
		logN_H2J5_fit_err_array = np.append(logN_H2J5_fit_err_array, log_to_real_err(logN_H2J5_fit, logN_H2J5_fit_err)) 
		logN_H2J6_fit_array = np.append(logN_H2J6_fit_array, 10**logN_H2J6_fit) 
		logN_H2J6_fit_err_array = np.append(logN_H2J6_fit_err_array, log_to_real_err(logN_H2J6_fit, logN_H2J6_fit_err)) 
		logN_H2J7_fit_array = np.append(logN_H2J7_fit_array, 10**logN_H2J7_fit) 
		logN_H2J7_fit_err_array = np.append(logN_H2J7_fit_err_array, log_to_real_err(logN_H2J7_fit, logN_H2J7_fit_err)) 
		logN_HDJ0_fit_array = np.append(logN_HDJ0_fit_array, 10**logN_HDJ0_fit) 
		logN_HDJ0_fit_err_array = np.append(logN_HDJ0_fit_err_array, log_to_real_err(logN_HDJ0_fit, logN_HDJ0_fit_err)) 
		logN_HDJ1_fit_array = np.append(logN_HDJ1_fit_array, 10**logN_HDJ1_fit) 
		logN_HDJ1_fit_err_array = np.append(logN_HDJ1_fit_err_array, log_to_real_err(logN_HDJ1_fit, logN_HDJ1_fit_err))
		logN_HDJ2_fit_array = np.append(logN_HDJ2_fit_array, 10**logN_HDJ2_fit) 
		logN_HDJ2_fit_err_array = np.append(logN_HDJ2_fit_err_array, log_to_real_err(logN_HDJ2_fit, logN_HDJ2_fit_err)) 
		b_val_array = np.append(b_val_array, b_val_fitted_new)
		b_val_array_err = np.append(b_val_array_err, b_val_fitted_new_err)


	
	global b_val_new_array, b_val_new_array_err
	b_val_new_array = np.zeros([11, 13])
	b_val_new_array_err = np.zeros([11, 13])


	for i in range(len(b_val_new)):
		b_val_new_array[i,:] = np.array([np.log10(logN_HI_fit_array[i]), np.log10(logN_H2J0_fit_array[i]), np.log10(logN_H2J1_fit_array[i]), np.log10(logN_H2J2_fit_array[i]), np.log10(logN_H2J3_fit_array[i]), np.log10(logN_H2J4_fit_array[i]), np.log10(logN_H2J5_fit_array[i]), np.log10(logN_H2J6_fit_array[i]), np.log10(logN_H2J7_fit_array[i]), np.log10(logN_HDJ0_fit_array[i]), np.log10(logN_HDJ1_fit_array[i]), np.log10(logN_HDJ2_fit_array[i]), b_val_array[i]])

		b_val_new_array_err[i,:] = np.array([real_to_log_err(logN_HI_fit_array[i], logN_HI_fit_err_array[i]), real_to_log_err(logN_H2J0_fit_array[i], logN_H2J0_fit_err_array[i]), real_to_log_err(logN_H2J1_fit_array[i], logN_H2J1_fit_err_array[i]), real_to_log_err(logN_H2J2_fit_array[i], logN_H2J2_fit_err_array[i]), real_to_log_err(logN_H2J3_fit_array[i], logN_H2J3_fit_err_array[i]), real_to_log_err(logN_H2J4_fit_array[i], logN_H2J4_fit_err_array[i]), real_to_log_err(logN_H2J5_fit_array[i], logN_H2J5_fit_err_array[i]), real_to_log_err(logN_H2J6_fit_array[i], logN_H2J6_fit_err_array[i]), real_to_log_err(logN_H2J7_fit_array[i], logN_H2J7_fit_err_array[i]), real_to_log_err(logN_HDJ0_fit_array[i], logN_HDJ0_fit_err_array[i]), real_to_log_err(logN_HDJ1_fit_array[i], logN_HDJ1_fit_err_array[i]), real_to_log_err(logN_HDJ2_fit_array[i], logN_HDJ2_fit_err_array[i]), b_val_array_err[i]])


	

	
	logN_HI_fit_min_idx = np.where(logN_HI_fit_array == np.min(logN_HI_fit_array))[0][0]
	
	logN_HI_fit_min = logN_HI_fit_array[logN_HI_fit_min_idx] - logN_HI_fit_err_array[logN_HI_fit_min_idx]

	logN_H2J0_fit_array_idx = np.where(logN_H2J0_fit_array == np.min(logN_H2J0_fit_array))[0][0]
	logN_H2J0_fit_array_min = logN_H2J0_fit_array[logN_H2J0_fit_array_idx] - logN_H2J0_fit_err_array[logN_H2J0_fit_array_idx]

	logN_H2J1_fit_array_idx = np.where(logN_H2J1_fit_array == np.min(logN_H2J1_fit_array))[0][0]
	logN_H2J1_fit_array_min = logN_H2J1_fit_array[logN_H2J1_fit_array_idx] - logN_H2J1_fit_err_array[logN_H2J1_fit_array_idx]
	
	logN_H2J2_fit_array_idx = np.where(logN_H2J2_fit_array == np.min(logN_H2J2_fit_array))[0][0]
	logN_H2J2_fit_array_min = logN_H2J2_fit_array[logN_H2J2_fit_array_idx] - logN_H2J2_fit_err_array[logN_H2J2_fit_array_idx]

	logN_H2J3_fit_array_idx = np.where(logN_H2J3_fit_array == np.min(logN_H2J3_fit_array))[0][0]
	logN_H2J3_fit_array_min = logN_H2J3_fit_array[logN_H2J3_fit_array_idx] - logN_H2J3_fit_err_array[logN_H2J3_fit_array_idx]

	logN_H2J4_fit_array_idx = np.where(logN_H2J4_fit_array == np.min(logN_H2J4_fit_array))[0][0]
	logN_H2J4_fit_array_min = logN_H2J4_fit_array[logN_H2J4_fit_array_idx] - logN_H2J4_fit_err_array[logN_H2J4_fit_array_idx]

	logN_H2J5_fit_array_idx = np.where(logN_H2J5_fit_array == np.min(logN_H2J5_fit_array))[0][0]
	logN_H2J5_fit_array_min = logN_H2J5_fit_array[logN_H2J5_fit_array_idx] - logN_H2J5_fit_err_array[logN_H2J5_fit_array_idx]

	logN_H2J6_fit_array_idx = np.where(logN_H2J6_fit_array == np.min(logN_H2J6_fit_array))[0][0]
	logN_H2J6_fit_array_min = logN_H2J6_fit_array[logN_H2J6_fit_array_idx] - logN_H2J6_fit_err_array[logN_H2J6_fit_array_idx]

	logN_H2J7_fit_array_idx = np.where(logN_H2J7_fit_array == np.min(logN_H2J7_fit_array))[0][0]
	logN_H2J7_fit_array_min = logN_H2J7_fit_array[logN_H2J7_fit_array_idx] - logN_H2J7_fit_err_array[logN_H2J7_fit_array_idx]

	logN_HDJ0_fit_array_idx = np.where(logN_HDJ0_fit_array == np.min(logN_HDJ0_fit_array))[0][0]
	logN_HDJ0_fit_array_min = logN_HDJ0_fit_array[logN_HDJ0_fit_array_idx] - logN_HDJ0_fit_err_array[logN_HDJ0_fit_array_idx]

	logN_HDJ1_fit_array_idx = np.where(logN_HDJ1_fit_array == np.min(logN_HDJ1_fit_array))[0][0]
	logN_HDJ1_fit_array_min = logN_HDJ1_fit_array[logN_HDJ1_fit_array_idx] - logN_HDJ1_fit_err_array[logN_HDJ1_fit_array_idx]


	logN_HDJ2_fit_array_idx = np.where(logN_HDJ2_fit_array == np.min(logN_HDJ2_fit_array))[0][0]

	if ((np.abs(logN_HDJ2_fit_err_array[logN_HDJ2_fit_array_idx] - logN_HDJ2_fit_array[logN_HDJ2_fit_array_idx]))<22):
		logN_HDJ2_fit_array_min = logN_HDJ2_fit_array[logN_HDJ2_fit_array_idx] - logN_HDJ2_fit_err_array[logN_HDJ2_fit_array_idx]
	else:
		logN_HDJ2_fit_array_min = logN_HDJ2_fit_array[logN_HDJ2_fit_array_idx] - 1e-5


	
	logN_HI_fit_max_idx2 = np.where(logN_HI_fit_array == np.max(logN_HI_fit_array))[0][0]
	logN_HI_fit_max = logN_HI_fit_array[logN_HI_fit_max_idx2] + logN_HI_fit_err_array[logN_HI_fit_max_idx2]
	
	logN_H2J0_fit_array_idx2 = np.where(logN_H2J0_fit_array == np.max(logN_H2J0_fit_array))[0][0]
	logN_H2J0_fit_array_max = logN_H2J0_fit_array[logN_H2J0_fit_array_idx2] + logN_H2J0_fit_err_array[logN_H2J0_fit_array_idx2]

	logN_H2J1_fit_array_idx2 = np.where(logN_H2J1_fit_array == np.max(logN_H2J1_fit_array))[0][0]
	logN_H2J1_fit_array_max = logN_H2J1_fit_array[logN_H2J1_fit_array_idx2] + logN_H2J1_fit_err_array[logN_H2J1_fit_array_idx2]
	
	logN_H2J2_fit_array_idx2 = np.where(logN_H2J2_fit_array == np.max(logN_H2J2_fit_array))[0][0]
	logN_H2J2_fit_array_max = logN_H2J2_fit_array[logN_H2J2_fit_array_idx2] + logN_H2J2_fit_err_array[logN_H2J2_fit_array_idx2]

	logN_H2J3_fit_array_idx2 = np.where(logN_H2J3_fit_array == np.max(logN_H2J3_fit_array))[0][0]
	logN_H2J3_fit_array_max = logN_H2J3_fit_array[logN_H2J3_fit_array_idx2] + logN_H2J3_fit_err_array[logN_H2J3_fit_array_idx2]

	logN_H2J4_fit_array_idx2 = np.where(logN_H2J4_fit_array == np.max(logN_H2J4_fit_array))[0][0]
	logN_H2J4_fit_array_max = logN_H2J4_fit_array[logN_H2J4_fit_array_idx2] + logN_H2J4_fit_err_array[logN_H2J4_fit_array_idx2]

	logN_H2J5_fit_array_idx2 = np.where(logN_H2J5_fit_array == np.max(logN_H2J5_fit_array))[0][0]
	logN_H2J5_fit_array_max = logN_H2J5_fit_array[logN_H2J5_fit_array_idx2] + logN_H2J5_fit_err_array[logN_H2J5_fit_array_idx2]

	logN_H2J6_fit_array_idx2 = np.where(logN_H2J6_fit_array == np.max(logN_H2J6_fit_array))[0][0]
	logN_H2J6_fit_array_max = logN_H2J6_fit_array[logN_H2J6_fit_array_idx2] + logN_H2J6_fit_err_array[logN_H2J6_fit_array_idx2]

	logN_H2J7_fit_array_idx2 = np.where(logN_H2J7_fit_array == np.max(logN_H2J7_fit_array))[0][0]
	logN_H2J7_fit_array_max = logN_H2J7_fit_array[logN_H2J7_fit_array_idx2] + logN_H2J7_fit_err_array[logN_H2J7_fit_array_idx2]

	logN_HDJ0_fit_array_idx2 = np.where(logN_HDJ0_fit_array == np.max(logN_HDJ0_fit_array))[0][0]
	logN_HDJ0_fit_array_max = logN_HDJ0_fit_array[logN_HDJ0_fit_array_idx2] + logN_HDJ0_fit_err_array[logN_HDJ0_fit_array_idx2]

	logN_HDJ1_fit_array_idx2 = np.where(logN_HDJ1_fit_array == np.max(logN_HDJ1_fit_array))[0][0]
	logN_HDJ1_fit_array_max = logN_HDJ1_fit_array[logN_HDJ1_fit_array_idx2] + logN_HDJ1_fit_err_array[logN_HDJ1_fit_array_idx2]

	logN_HDJ2_fit_array_idx2 = np.where(logN_HDJ2_fit_array == np.max(logN_HDJ2_fit_array))[0][0]

	if ((np.abs(logN_HDJ2_fit_err_array[logN_HDJ2_fit_array_idx2] - logN_HDJ2_fit_array[logN_HDJ2_fit_array_idx2]))<22):
		logN_HDJ2_fit_array_max = logN_HDJ2_fit_array[logN_HDJ2_fit_array_idx2] + logN_HDJ2_fit_err_array[logN_HDJ2_fit_array_idx2]
	else:
		logN_HDJ2_fit_array_max = logN_HDJ2_fit_array[logN_HDJ2_fit_array_idx2] + 1e-5

	
	N_HI_total = ((logN_HI_fit_min + logN_HI_fit_max)/2)
	N_HI_total_err = np.abs( (logN_HI_fit_max) - ((logN_HI_fit_min + logN_HI_fit_max)/2))

	N_H2J0_total = ((logN_H2J0_fit_array_min + logN_H2J0_fit_array_max)/2)
	N_H2J0_total_err = np.abs( (logN_H2J0_fit_array_max) - ((logN_H2J0_fit_array_min + logN_H2J0_fit_array_max)/2))

	N_H2J1_total = ((logN_H2J1_fit_array_min + logN_H2J1_fit_array_max)/2)
	N_H2J1_total_err = np.abs( (logN_H2J1_fit_array_max) - ((logN_H2J1_fit_array_min + logN_H2J1_fit_array_max)/2))

	N_H2J2_total = ((logN_H2J2_fit_array_min + logN_H2J2_fit_array_max)/2)
	N_H2J2_total_err = np.abs( (logN_H2J2_fit_array_max) - ((logN_H2J2_fit_array_min + logN_H2J2_fit_array_max)/2))

	N_H2J3_total = ((logN_H2J3_fit_array_min + logN_H2J3_fit_array_max)/2)
	N_H2J3_total_err = np.abs( (logN_H2J3_fit_array_max) - ((logN_H2J3_fit_array_min + logN_H2J3_fit_array_max)/2))

	N_H2J4_total = ((logN_H2J4_fit_array_min + logN_H2J4_fit_array_max)/2)
	N_H2J4_total_err = np.abs( (logN_H2J4_fit_array_max) - ((logN_H2J4_fit_array_min + logN_H2J4_fit_array_max)/2))

	N_H2J5_total = ((logN_H2J5_fit_array_min + logN_H2J5_fit_array_max)/2)
	N_H2J5_total_err = np.abs( (logN_H2J5_fit_array_max) - ((logN_H2J5_fit_array_min + logN_H2J5_fit_array_max)/2))
		
	N_H2J6_total = ((logN_H2J6_fit_array_min + logN_H2J6_fit_array_max)/2)
	N_H2J6_total_err = np.abs( (logN_H2J6_fit_array_max) - ((logN_H2J6_fit_array_min + logN_H2J6_fit_array_max)/2))

	N_H2J7_total = ((logN_H2J7_fit_array_min + logN_H2J7_fit_array_max)/2)
	N_H2J7_total_err = np.abs( (logN_H2J7_fit_array_max) - ((logN_H2J7_fit_array_min + logN_H2J7_fit_array_max)/2))

	N_HDJ0_total = ((logN_HDJ0_fit_array_min + logN_HDJ0_fit_array_max)/2)
	N_HDJ0_total_err = np.abs( (logN_HDJ0_fit_array_max) - ((logN_HDJ0_fit_array_min + logN_HDJ0_fit_array_max)/2))
	
	N_HDJ1_total = ((logN_HDJ1_fit_array_min + logN_HDJ1_fit_array_max)/2)
	N_HDJ1_total_err = np.abs( (logN_HDJ1_fit_array_max) - ((logN_HDJ1_fit_array_min + logN_HDJ1_fit_array_max)/2))

	N_HDJ2_total = ((logN_HDJ2_fit_array_min + logN_HDJ2_fit_array_max)/2)
	N_HDJ2_total_err = np.abs( (logN_HDJ2_fit_array_max) - ((logN_HDJ2_fit_array_min + logN_HDJ2_fit_array_max)/2))


	N_HI_total_log = np.log10(N_HI_total)
	N_HI_total_log_err = real_to_log_err(N_HI_total, N_HI_total_err)
	N_H2J0_total_log = np.log10(N_H2J0_total)
	N_H2J0_total_log_err = real_to_log_err(N_H2J0_total, N_H2J0_total_err)
	N_H2J1_total_log = np.log10(N_H2J1_total)
	N_H2J1_total_log_err = real_to_log_err(N_H2J1_total, N_H2J1_total_err)
	N_H2J2_total_log = np.log10(N_H2J2_total)
	N_H2J2_total_log_err = real_to_log_err(N_H2J2_total, N_H2J2_total_err)
	N_H2J3_total_log = np.log10(N_H2J3_total)
	N_H2J3_total_log_err = real_to_log_err(N_H2J3_total, N_H2J3_total_err)
	N_H2J4_total_log = np.log10(N_H2J4_total)
	N_H2J4_total_log_err = real_to_log_err(N_H2J4_total, N_H2J4_total_err)
	N_H2J5_total_log = np.log10(N_H2J5_total)
	N_H2J5_total_log_err = real_to_log_err(N_H2J5_total, N_H2J5_total_err)
	N_H2J6_total_log = np.log10(N_H2J6_total)
	N_H2J6_total_log_err = real_to_log_err(N_H2J6_total, N_H2J6_total_err)
	N_H2J7_total_log = np.log10(N_H2J7_total)
	N_H2J7_total_log_err = real_to_log_err(N_H2J7_total, N_H2J7_total_err)
	N_HDJ0_total_log = np.log10(N_HDJ0_total)
	N_HDJ0_total_log_err = real_to_log_err(N_HDJ0_total, N_HDJ0_total_err)
	N_HDJ1_total_log = np.log10(N_HDJ1_total)
	N_HDJ1_total_log_err = real_to_log_err(N_HDJ1_total, N_HDJ1_total_err)
	N_HDJ2_total_log = np.log10(N_HDJ2_total)
	N_HDJ2_total_log_err = real_to_log_err(N_HDJ2_total, N_HDJ2_total_err)



	b_val_new_array[10,:] = np.array([N_HI_total_log, N_H2J0_total_log, N_H2J1_total_log, N_H2J2_total_log, N_H2J3_total_log, N_H2J4_total_log, N_H2J5_total_log, N_H2J6_total_log, N_H2J7_total_log, N_HDJ0_total_log, N_HDJ1_total_log, N_HDJ2_total_log, b_val_array[0]])

	b_val_new_array_err[10,:] = np.array([N_HI_total_log_err, N_H2J0_total_log_err, N_H2J1_total_log_err, N_H2J2_total_log_err, N_H2J3_total_log_err, N_H2J4_total_log_err, N_H2J5_total_log_err, N_H2J6_total_log_err, N_H2J7_total_log_err, N_HDJ0_total_log_err, N_HDJ1_total_log_err, N_HDJ2_total_log_err, b_val_array_err[0]])

	print (b_val_new_array)
	print (b_val_new_array_err)
	print ('Success...')


	file_name_b_val = file_name1[:-9] + '_b_val_details.txt'
	file_name_b_val_err = file_name1[:-9] + '_b_val_details_err.txt'

	save_prompt3 = raw_input("Save file (y/n): ")
	if (save_prompt3.upper()=='Y'):
		np.savetxt(file_name_b_val, b_val_new_array)
		np.savetxt(file_name_b_val_err, b_val_new_array_err)
		print ('Saved')
	else:
		print ('Not Saved')


		
button16.on_clicked(fit_b_val_ax)






#Saving the above mentioned b-value fit 

save_b_val_ax = plt.axes([0.005, 0.400, 0.07, 0.02])
button17 = Button(save_b_val_ax, 'save_b_val')

def save_b_val_ax(event):
	print ('Test2....')

	file_name_b_val = file_name1[:-9] + '_b_val_details.txt'
	file_name_b_val_err = file_name1[:-9] + '_b_val_details_err.txt'

	np.savetxt(file_name_b_val, b_val_new_array)
	np.savetxt(file_name_b_val_err, b_val_new_array_err)

	print ('saved')




button17.on_clicked(save_b_val_ax)





#Loading the above mentioned b-value fit, if previously fitted 

get_b_val_ax = plt.axes([0.005, 0.375, 0.07, 0.02])
button18 = Button(get_b_val_ax, 'get_b_val')

def get_b_val_ax(event):
	print ('Test3....')

	file_name_b_val = file_name1[:-9] + '_b_val_details.txt'
	file_name_b_val_err = file_name1[:-9] + '_b_val_details_err.txt'

	global b_val_new_array, b_val_new_array_err
	b_val_new_array = np.loadtxt(file_name_b_val)
	b_val_new_array_err = np.loadtxt(file_name_b_val_err)
	print ('Loaded')


button18.on_clicked(get_b_val_ax)



#Making figure of the above mentioned b-value fit 

b_val_fig_ax = plt.axes([0.005, 0.350, 0.07, 0.02])
button19 = Button(b_val_fig_ax, 'relevant_figs')

def b_val_fig_ax(event):
	print ('Test4....')

	file_name_b_val = file_name1[:-9] + '_b_val_details.txt'
	file_name_b_val_err = file_name1[:-9] + '_b_val_details_err.txt'


	if (os.path.isfile(file_name_b_val)):

		str_prog_code3 = 'python ' + 'plotting_files/plot_b_val_vs_col_den.py' + ' ' + str(file_name_b_val) + ' ' + str(file_name_b_val_err)
		os.system(str_prog_code3)
		print ('Plotted b-val vs column-density....')

		b_val_new_array = np.loadtxt(file_name_b_val)
		b_val_new_array_err = np.loadtxt(file_name_b_val_err)


		b_val_new_array = np.nan_to_num(b_val_new_array)
		b_val_new_array_err = np.nan_to_num(b_val_new_array_err)

		print (b_val_new_array)
		print (b_val_new_array_err)
		


		N_HI_total_log2, N_H2J0_total_log2, N_H2J1_total_log2, N_H2J2_total_log2, N_H2J3_total_log2, N_H2J4_total_log2, N_H2J5_total_log2, N_H2J6_total_log2, N_H2J7_total_log2, N_HDJ0_total_log2, N_HDJ1_total_log2, N_HDJ2_total_log2, b_val_array2 = b_val_new_array[10,:]

		N_HI_total_log_err2, N_H2J0_total_log_err2, N_H2J1_total_log_err2, N_H2J2_total_log_err2, N_H2J3_total_log_err2, N_H2J4_total_log_err2, N_H2J5_total_log_err2, N_H2J6_total_log_err2, N_H2J7_total_log_err2, N_HDJ0_total_log_err2, N_HDJ1_total_log_err2, N_HDJ2_total_log_err2, b_val_array_err2 = b_val_new_array_err[10,:]
		


		file_name_b_val_excitation = file_name1[:-9] + '_b_val_excitation_diagram' + '.pdf'

		str_prog_code2 = 'python ' + 'plotting_files/excitation_diagram2.py' + ' ' + str(N_H2J0_total_log2) + ' ' + str(N_H2J0_total_log_err2) + ' ' + str(N_H2J1_total_log2) + ' ' + str(N_H2J1_total_log_err2) + ' ' + str(N_H2J2_total_log2) + ' ' + str(N_H2J2_total_log_err2) + ' ' + str(N_H2J3_total_log2) + ' ' + str(N_H2J3_total_log_err2) + ' ' + str(N_H2J4_total_log2) + ' ' + str(N_H2J4_total_log_err2) + ' ' + str(N_H2J5_total_log2) + ' ' + str(N_H2J5_total_log_err2) + ' ' + str(N_H2J6_total_log2) + ' ' + str(N_H2J6_total_log_err2) + ' ' + str(N_H2J7_total_log2) + ' ' + str(N_H2J7_total_log_err2) + ' ' + str(file_name_b_val_excitation)
		print (str_prog_code2)
		os.system(str_prog_code2)
		

		print ('Plotted excitation diagram....')
		print ('Success')

	else:
		print ('File not found....')


button19.on_clicked(b_val_fig_ax)




#Making table of the above mentioned b-value fit 

b_val_table_ax = plt.axes([0.005, 0.325, 0.07, 0.02])
button20 = Button(b_val_table_ax, 'b_val_table')

def b_val_table_ax(event):
	print ('Test5....')


	file_name_b_val = file_name1[:-9] + '_b_val_details.txt'
	file_name_b_val_err = file_name1[:-9] + '_b_val_details_err.txt'


	if (os.path.isfile(file_name_b_val)):

		str_prog_code3 = 'python ' + 'plotting_files/make_table_b_val_vs_col_den.py' + ' ' + str(file_name_b_val) + ' ' + str(file_name_b_val_err)
		os.system(str_prog_code3)
		print ('Table for b-val vs column-density created....')
		print ('Success')
	else:
		print ('File not found....')

	

button20.on_clicked(b_val_table_ax)







ax.set_xlim((xdata.min()-(0.02*(xdata.mean()))), (xdata.max()+(0.02*(xdata.mean()))))
ax.set_ylim((ydata.min()-(0.1*(ydata.mean()))), (ydata.max()+(1.6*(ydata.mean()))))

plt.show(1)
