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



##########################################FUNCTION_DEFINITIONS#########################################


##################################CHEBYSHEV_FUNCTIONS##################################

def chebyshev_order(wave, cont, stopping_number):
    wave_new = np.linspace(-1, 1, len(wave))
    i=1
    while True:
        roots = numpy.polynomial.chebyshev.chebfit(wave_new, cont, i, rcond=None, full=False, w=None)
        poly = numpy.polynomial.chebyshev.chebval(wave_new, roots, tensor=True)
        chi_sq = (poly - cont) ** 2
        chi_sq_sum = (np.sum(chi_sq))
        i+=1
        if not chi_sq_sum>(float(stopping_number)):
            break
    return (i, roots)

def chebyshev_fit(wave, cont, order):
    wave_new = np.linspace(-1, 1, len(wave))
    
    roots = numpy.polynomial.chebyshev.chebfit(wave_new, cont, order_test, rcond=None, full=False, w=None)
    poly_new = numpy.polynomial.chebyshev.chebval(wave_new, roots, tensor=True)
    return(roots, poly_new)

def chebyshev_disp(wave, coeff):
    wave_new = np.linspace(-1, 1, len(wave))
    poly_new = numpy.polynomial.chebyshev.chebval(wave_new, coeff, tensor=True)
    return(poly_new)

##################################CHEBYSHEV_FUNCTIONS##################################



############################CONVOLUTION#####################################

def convolved_prof3(wave, profile, res):
	wave_short = np.logspace(np.log10(wave.min()), np.log10(wave.max()), len(wave))
	center = np.searchsorted(wave_short, np.log10(np.median(wave_short)))
    	deltalam = wave_short[center + 1] - wave_short[center]
    	sigma = (wave_short[center]) / (res * (2 * np.sqrt(2 * np.log(2))) * deltalam)
    	gauss = scipy.signal.gaussian(len(profile), sigma, sym=True)
    	gauss = gauss / np.sum(gauss)
    	prof_new = signal.fftconvolve(profile, gauss, mode='same')
    	return (prof_new)

############################CONVOLUTION#####################################



######################VELOCITY PROFILE##################################

def vel_prof(x, centre):
	xnew = (c*1e-5) * ((x-centre)/x)	
	return (xnew)

######################VELOCITY PROFILE##################################


######################FIND_THE_NEAREST_VALUE##################################


def find_nearest(array,value):
	idx = (np.abs(array-value)).argmin()
	return array[idx]

######################FIND_THE_NEAREST_VALUE##################################


######################ARRAY_SMOOTHING_FUNCTION##################################

def smooth(y, box_pts):
	box = np.ones(box_pts)/box_pts
	y_smooth = np.convolve(y, box, mode='same')
	return y_smooth

######################ARRAY_SMOOTHING_FUNCTION##################################


##################################EXTRACTING_ATOMIC_DATA_LIST_FROM_ELEMENT##################################

initial_guesses = np.genfromtxt('initial_guess.dat', dtype=[('mystring','S30')], comments='#')
atomic_database_file = str(initial_guesses[35][0])


def atomic_list_extractor(atom_name):
	if (isinstance(atom_name, basestring)):
		atom_name = atom_name.split(',')
	readdata = csv.reader(open(str(atomic_database_file), 'r'))
	data = []
	for row in readdata:
		data.append(row)
	Header = data[0]
	del data[0]
	species_name = np.chararray([], itemsize=10)
	species_identification = np.chararray([], itemsize=10)
	species_wavelength = np.zeros([])
	species_oscillator_strength = np.zeros([])
	species_tau_value = np.zeros([])
	for i in range(len(data)):
		for j in range(len(atom_name)):
			if (atom_name[j] in (str(data[i][0].split()[0]) + str(int(float(str(data[i][0].split()[1])))))):
				species_name = np.append(species_name, (str(data[i][0].split()[0]) + str(int(float(str(data[i][0].split()[1]))))))
				species_wavelength = np.append(species_wavelength, float(str(data[i][0].split()[1])))
				species_oscillator_strength = np.append(species_oscillator_strength, float(str(data[i][0].split()[2])))
				species_tau_value = np.append(species_tau_value, float(str(data[i][0].split()[3])))
	species_name = np.delete(species_name, 0)
	species_wavelength = np.delete(species_wavelength, 0)
	species_oscillator_strength = np.delete(species_oscillator_strength, 0)
	species_tau_value = np.delete(species_tau_value, 0)
	return (species_name, species_wavelength, species_oscillator_strength, species_tau_value)


##################################EXTRACTING_ATOMIC_DATA_LIST_FROM_ELEMENT##################################


############################VOIGT_PROFILE_NEW#####################################
#A HUGE THANKS TO JK FOR THIS MODULE

def H(a, x):
	P = x**2
	H0 = np.exp(-x**2)
	Q = 1.5/x**2
	return H0 - a/np.sqrt(np.pi)/P * (H0*H0*(4.*P*P + 7.*P + 4. + Q) - Q - 1)

def Voigt3_H2(l, l0, f, N, b, gam, z, resolution):
	"""Calculate the Voigt profile of transition with
	rest frame transition wavelength: 'l0'
	oscillator strength: 'f'
	column density: N  cm^-2
	velocity width: b  cm/s
	"""
	# ==== PARAMETERS ==================
	c = 2.99792e10        # cm/s
	m_e = 9.1095e-28       # g
	e = 4.8032e-10        # cgs units
	# ==================================
	# Calculate Profile
	C_a = np.sqrt(np.pi)*e**2*f*l0*1.e-8/m_e/c/b
	a = l0*1.e-8*gam/(4.*np.pi*b)
	dl_D = b/c*l0
	l = l/(z+1.)
	x = (l - l0)/dl_D+0.0001
	tau = np.float64(C_a) * N * H(a, x)
	return (tau)



def group_voigt2_H2(wave, osc, tau, x_axis, logN, b, z, resolution):
	profile = 0.
	for i in range(len(wave)):
		if ( (((1+z[i])*wave[i])>x_axis.min()) and (((1+z[i])*wave[i])<x_axis.max()) ):
			profile += Voigt3_H2(x_axis, float(wave[i]), float(osc[i]), 10**logN[i], b[i]*1e5, float(tau[i]), z[i], resolution)
	profile = (np.exp(-profile))
	res = resolution
	profile_conv = convolved_prof3(x_axis, profile, res)
	return (profile_conv)


############################VOIGT_PROFILE_NEW#####################################





