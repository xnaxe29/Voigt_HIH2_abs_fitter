#python plotting_HI.py J0017+1307_uvb_fitted_combined.spec J0017+1307_uvb_selected_combined.spec 2.32415496741 2.32687235267 subplot_grid_custom


import numpy as np
import matplotlib.pyplot as plt
import sys
import itertools
import os.path
import os
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.gridspec as gridspec
plt.rcParams["font.family"] = "Times New Roman"
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker

styles = mpatches.ArrowStyle.get_styles()
c = 299792.458 # Speed in Light in Km/s


str1 = sys.argv[1]
str2 = sys.argv[2]
z_HI = float(sys.argv[3])
z_H2 = float(sys.argv[4])
subplot_grid = str(sys.argv[5])
name_DLA = str1[15:10]


#Setting y-scale limit and size of font
ylim_init = -0.9
ylim_end = 10
size_of_font = 17



######################GET THE WIDTH/HEIGHT OF A SUBPLOT########################

def get_ax_size(ax):
    bbox = ax.get_window_extent().transformed(f.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= f.dpi
    height *= f.dpi
    return width, height

######################GET THE WIDTH/HEIGHT OF A SUBPLOT########################




######################CAN THE NUMBER YIELD AN INTEGER SQUARE-ROOT?########################

def is_square(integer):
    root = math.sqrt(integer)
    if int(root + 0.5) ** 2 == integer: 
        return True
    else:
        return False

######################CAN THE NUMBER YIELD AN INTEGER SQUARE-ROOT?########################





######################IS THE NUMBER/STRING REALLY A FLOAT?########################

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

######################IS THE NUMBER/STRING REALLY A FLOAT?########################




######################ADDING SUBPLOTS#########################################

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3] 
    subax = fig.add_axes([x,y,width,height],facecolor=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

######################ADDING SUBPLOTS#########################################


######################EXTRACTOR#########################################
def extractor(namer):

	# Open file
	f5 = open('plotting_files/atomic_database_file.dat', 'r')

	# Loop over lines and extract variables of interest
	name = np.chararray([])
	namenew = np.chararray([])
	wavelength = np.array([])
	osc_str = np.array([])
	tau_value = np.array([])
	for line in f5:
		line = line.strip()
		if not line.startswith("#") and not line.startswith("!"):
			columns = line.split()
			name = np.append(name, columns[0])
			namenew = np.append(namenew, str(str(columns[0])+str(int(float(columns[1])))))
			wavelength = np.append(wavelength, columns[1])
			osc_str = np.append(osc_str, columns[2])
			tau_value = np.append(tau_value, columns[3])


	name = np.delete(name, 0)
	namenew = np.delete(namenew, 0)
	
	# Make a list of the required molecules with wavelength, oscillator strength and tau value
	name_list = np.array([])
	wave_list = np.array([])
	ocs_str_list = np.array([])
	tau_val_list = np.array([])

	
	for i in range(0, len(namer)):
		for j in range(0, len(namenew)):
			if (namer[i]==namenew[j]):
				name_list = np.append(name_list, namenew[j])
				wave_list = np.append(wave_list, wavelength[j].astype(np.float))
				ocs_str_list = np.append(ocs_str_list, osc_str[j].astype(np.float))
				tau_val_list = np.append(tau_val_list, tau_value[j].astype(np.float))
	
	


	f5.close()
	return (name_list, wave_list, ocs_str_list, tau_val_list)


######################EXTRACTOR#########################################


######################VELOCITY PROFILE##################################

def vel_prof(x, centre):
	xnew = c * ((x-centre)/x)	
	return (xnew)

######################VELOCITY PROFILE##################################




######################FIND THE NEAREST NUMBER###########################

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

######################FIND THE NEAREST NUMBER###########################




'''
#####OPEN THE CUSTOM TXT FILE CONTAINING NAMES OF THE METALS TO BE PLOTTED#####
b = open('list_H2J0_H2J1_custom.txt', 'r').readlines()
name = np.chararray([len(b)], itemsize=10)


for i in range(len(b)):
	name[i] = str(b[i].split()[0])


namer = name
name_list, wave_list, ocs_str_list, tau_val_list = extractor(namer)
waver2=wave_list
#####OPEN THE CUSTOM TXT FILE CONTAINING NAMES OF THE METALS TO BE PLOTTED#####
'''

#Loading the data
H2_wave, H2_flux, H2_flux_err, H2_cont, H2_prof = np.loadtxt(str1, unpack=True)
wave_selected, flux_selected, err_selected, cont_selected = np.loadtxt(str2, unpack=True)





wave = H2_wave
data = H2_flux
err = H2_flux_err
fit = H2_prof
zr = z_H2


y = data
err = err
z = fit
fit_area = wave_selected
fit_area_new = flux_selected
fit_area_new_err = err_selected






#Converting wavelength to velocity
wave_new2 = wave/(1+z_HI)
test_x = vel_prof(wave_new2, 1215.6701)
wave_new3 = fit_area/(1+z_HI)
test_x3 = vel_prof(wave_new3, 1215.6701)




#Plotting

fig, ax_new = plt.subplots(figsize=(25, 16), dpi=300)
ax_new.errorbar(test_x, y, yerr=err, fmt='k-', alpha=0.2, drawstyle='steps-mid', ecolor='k', elinewidth=0.5, capsize=1.0, capthick=1, label='data', zorder=1)
#plt.errorbar(test_x3, fit_area_new, yerr=fit_area_new_err, fmt='k.', alpha=0.8, drawstyle='steps-mid', ecolor='k', elinewidth=0.5, capsize=1.0, capthick=1, label='selected data', zorder=2)
ax_new.plot(test_x3, fit_area_new, 'k.', alpha=0.8, label='selected data', zorder=2)
ax_new.plot(test_x,z, 'r-', alpha=0.8, label='fit', zorder=3)
ax_new.plot(test_x, H2_cont, alpha=0.8, label='continuum', zorder=4)




#Setting x/y scale limits and tick params
 
#plt.locator_params(nbins=4, axis='x')
#plt.tick_params(axis='both', which='major', labelsize=1.5*size_of_font)
#plt.tick_params(axis='both', which='minor', labelsize=1.5*size_of_font)		
plt.xlim(-20000, 20000)
plt.ylim(ylim_init,ylim_end)
plt.axvline(0, color='green', linestyle='dashed', alpha=0.5, lw=2)
plt.axhline(0, color='green', linestyle='dashed', alpha=0.5, lw=2)
plt.axhline(1, color='green', linestyle='dashed', alpha=0.5, lw=2)		
#plt.tick_params(axis='both', which='major', direction='in', length=10, labelsize=size_of_font)
#plt.tick_params(axis='both', which='minor', direction='in', length=5, labelsize=size_of_font)
ax_new.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax_new.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
ax_new.xaxis.set_major_locator(ticker.MultipleLocator(10000))
ax_new.xaxis.set_minor_locator(ticker.MultipleLocator(1000))
ax_new.tick_params(axis = 'both', which = 'major', direction='in', length=15, width=2, colors='k')
ax_new.tick_params(axis = 'both', which = 'minor', direction='in', length=7, width=1, colors='k')


#Setting miscellaneous plot params

plt.text(15000, 0.25, 'HI 1215', fontsize=size_of_font*3, color='black', alpha=1.0, ha='right')
plt.title(name_DLA, fontsize=size_of_font*3)
plt.legend(loc=3, fontsize = size_of_font*2)
plt.xlabel(r'Relative velocity (km s$^{-1}$)', fontsize=size_of_font*3)
plt.ylabel(r'Flux (10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)', fontsize=size_of_font*3)
plt.tick_params(axis='both', labelsize=size_of_font*3)
plt.margins(y=0.2, x=0.2)



#Saving the plot
saved_file_name = str1[:-25] + '_HI_1215.pdf'
plt.savefig(saved_file_name)
#plt.show()
