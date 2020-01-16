#python plotting_H2J0_H2J1.py J0017+1307_uvb_fitted_combined.spec J0017+1307_uvb_selected_combined.spec 2.32415496741 2.32687235267 subplot_grid_custom


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
name_DLA = str1[:10]



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
	f5 = open('atomic_database_file.dat', 'r')

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





#####OPEN THE CUSTOM TXT FILE CONTAINING NAMES OF THE METALS TO BE PLOTTED#####
b = open('plotting_files/list_H2J0_H2J1_custom.txt', 'r').readlines()
name = np.chararray([len(b)], itemsize=10)


for i in range(len(b)):
	name[i] = str(b[i].split()[0])


namer = name
name_list, wave_list, ocs_str_list, tau_val_list = extractor(namer)
waver2=wave_list
#####OPEN THE CUSTOM TXT FILE CONTAINING NAMES OF THE METALS TO BE PLOTTED#####


#Loading data file for the fit and the selected values

H2_wave, H2_flux, H2_flux_err, H2_cont, H2_prof = np.loadtxt(str1, unpack=True)
wave_selected, flux_selected, err_selected, cont_selected = np.loadtxt(str2, unpack=True)



#Removing divide by zero dependancies which occur commonly while normalising

for i in range(len(wave_selected)):
	if (cont_selected[i]>0.0):
		flux_selected[i] = flux_selected[i]/cont_selected[i]
		err_selected[i] = err_selected[i]/cont_selected[i]
	else:
		flux_selected[i] = np.nan
		err_selected[i] = np.nan



for i in range(len(H2_cont)):
	if H2_cont[i]>0.0:
		H2_flux[i] = H2_flux[i]/H2_cont[i]
		H2_flux_err[i] = H2_flux_err[i]/H2_cont[i]
		H2_prof[i] = H2_prof[i]/H2_cont[i]
	else:
		H2_flux[i] = np.nan
		H2_flux_err[i] = np.nan
		H2_prof[i] = np.nan




#Creating array to display the residuals of the fit

residual=np.zeros([(len(H2_wave))])
residual[:] = np.nan
for i in range(len(wave_selected)):
	if H2_flux_err[i]!=0:
		residual[i] = ((H2_flux[i] - H2_prof[i])/H2_flux_err[i])


residual_new=np.zeros([(len(H2_wave))])
residual_new[:] = np.nan
for i in range(len(wave_selected)):
	if H2_flux_err[i]!=0 or H2_flux_err[i]!=np.nan:
		residual_new[i] = ((flux_selected[i] - H2_prof[i])/err_selected[i])



residual_rev=np.zeros([(len(H2_wave))])
residual_rev[:] = np.nan
for i in range(len(H2_wave)):
	if H2_flux_err[i]!=0 or H2_flux_err[i]!=np.nan:
		residual_rev[i] = ((H2_flux[i] - H2_prof[i])/H2_flux_err[i])








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
residual_filtered = residual_new/12 + 1.5
residual_non_filtered = residual/12 + 1.5
residual_filtered = residual_rev




residual_filtered_squeezed = np.zeros([len(residual_filtered)])
residual_filtered_squeezed[:] = np.nan
for i in range(len(residual_filtered)):
	if (-2 < residual_filtered[i] < 2):
		#residual_filtered_squeezed[i] = (residual_filtered[i]/12)+1.6
		residual_filtered_squeezed[i] = (residual_filtered[i]/12.5)+1.6





#Setting y-scale limit and size of font
ylim_init=-0.1
ylim_end=1.5
size_of_font = 20


#####PLOT THE DATA IN A CUSTOM SUBPLOTTING ROUTINE CREATING FOR MAKING SQUARE SUBPLOTS#####

if (subplot_grid=='subplot_grid_auto'):
	slots = len(waver2)
	if (is_square(slots)):
		x_slot = np.int(np.sqrt(slots))
		y_slot = np.int(np.sqrt(slots))
	else:
		x_slot = np.int(np.ceil(np.sqrt(slots)))
		y_slot = np.int(np.ceil(np.sqrt(slots)))
elif (subplot_grid=='subplot_grid_custom'):
	x_slot = 8
	y_slot = 8
	

#####PLOT THE DATA IN A CUSTOM SUBPLOTTING ROUTINE CREATING FOR MAKING SQUARE SUBPLOTS#####


#Plotting

f, axarr = plt.subplots(x_slot, y_slot, sharex=True, sharey=True, figsize=((x_slot*4), (y_slot*4)), dpi=300)
#plt.locator_params(axis='x', nticks=6)
#f.tight_layout()
f.subplots_adjust(wspace=0.1, hspace=0.1)


#Setting tick params
for axi in axarr.flat:
	axi.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
	axi.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

	axi.xaxis.set_major_locator(ticker.MultipleLocator(150))
	axi.xaxis.set_minor_locator(ticker.MultipleLocator(30))

	axi.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')

	axi.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')



#Creating each individual subplot window

for u in range(x_slot):
	for v in range(y_slot):
		if (len(waver2) > ((x_slot*v)+u)):

			wave_new2 = wave/(1+zr)
			test_x = vel_prof(wave_new2, waver2[(x_slot*v)+u])
			
			wave_new3 = fit_area/(1+zr)
			test_x3 = vel_prof(wave_new3, waver2[(x_slot*v)+u])



	
			axarr[u, v].errorbar(test_x, y, yerr=err, fmt='k-', alpha=0.2, drawstyle='steps-mid', ecolor='k', elinewidth=0.5, capsize=1.0, capthick=1, label='data', zorder=1)
			axarr[u, v].errorbar(test_x3, fit_area_new, yerr=fit_area_new_err, fmt='k.', alpha=0.8, drawstyle='steps-mid', ecolor='k', elinewidth=0.5, capsize=1.0, capthick=1, label='data', zorder=2)
			axarr[u, v].plot(test_x,z, 'r-', alpha=0.8, label='fit', zorder=3)
	
			axarr[u, v].plot(test_x,residual_filtered_squeezed, 'r.', markersize=3.0)
			axarr[u, v].axhline(1.6, color='blue', linestyle='dashed', alpha=1.0, lw=1)
			axarr[u, v].axhline((1.6-0.16), color='green', linestyle='dashed', alpha=1.0, lw=1)
			axarr[u, v].text(180, (1.6-0.16), r'$\rm -2\sigma$', fontsize=size_of_font/2, color='green')
			axarr[u, v].text(180, (1.6+0.06), r'$\rm +2\sigma$', fontsize=size_of_font/2, color='green')
			axarr[u, v].axhline((1.6+0.16), color='green', linestyle='dashed', alpha=1.0, lw=1)
			#axarr[u, v].set_ylim(-2.0/12, 2.0/12)	
		

			axarr[u, v].set_xlim(-260, 260)
			axarr[u, v].set_ylim(ylim_init,ylim_end)
			axarr[u, v].axvline(0, color='green', linestyle='dashed', alpha=0.5, lw=1)
			axarr[u, v].axhline(0, color='green', linestyle='dashed', alpha=0.5, lw=1)
			axarr[u, v].axhline(1, color='green', linestyle='dashed', alpha=0.5, lw=1)
			
			






			axarr[u, v].tick_params(axis='both', which='major', direction='in', length=10, labelsize=size_of_font)
			axarr[u, v].tick_params(axis='both', which='minor', direction='in', length=5, labelsize=size_of_font)


			axarr[u, v].text(250, 1.25, name_list[(x_slot*v)+u], fontsize=size_of_font/1.1, color='black', alpha=1.0, ha='right')
			

		else:
			axarr[u, v].set_visible(False)
		

	for i in range(x_slot-1):
		plt.setp([a.get_xticklabels() for a in axarr[i, :]], visible=False)
		plt.setp([a.get_yticklabels() for a in axarr[:, (i+1)]], visible=False)




#Miscellaneous text

f.text(0.5, 0.05, r'Relative velocity (km s$\rm^{-1}$)', ha='center', va='center', fontsize=1.5*size_of_font)
f.text(0.05, 0.5, 'Normalized flux', ha='center', va='center', rotation='vertical', fontsize=1.5*size_of_font)
f.text(0.5, 0.92, name_DLA, ha='center', va='center', fontsize=size_of_font*1.5)
#plt.title(name_DLA, fontsize=size_of_font*1.5)




#Save the file
saved_file_name = str1[:-25] + '_H2J0_H2J1.pdf'
plt.savefig(saved_file_name)
#plt.show()
