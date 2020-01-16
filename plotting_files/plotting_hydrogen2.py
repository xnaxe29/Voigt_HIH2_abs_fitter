import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.ticker as ticker

str1 = sys.argv[1]
str2 = sys.argv[2]
z_abs = float(sys.argv[3])
size_of_font = 30
name_DLA = str1[:10]


#Arrow characteristics. Please change values to change the size of the arrow pointing for each rotatational level of H2.
top_lyman = 1.4
bottom_lyman = 1.1
top_werner = 1.8+0.1
werner_arrow = 1.7+0.1
werner_text = 1.85+0.1
even_lyman_arrow = 1.35
even_lyman_text = 1.45
odd_lyman_arrow = 1.05
odd_lyman_text = 1.15


#Custom values to define the window (in wavelength) for each subplot
idx1 = int(3195/(1+z_abs))
idx2 = int(3400/(1+z_abs))
idx3 = int(3580/(1+z_abs))
idx4 = int(3770/(1+z_abs))
idx5 = int(4250/(1+z_abs))

diff1 = ((3440-3195)/10)



# Loading the fitted data file
HI_wave, HI_flux, HI_flux_err, HI_cont, HI_prof = np.loadtxt(str1, unpack=True)
H2_wave, H2_flux, H2_flux_err, H2_cont, H2_prof = np.loadtxt(str1, unpack=True)
wave, flux, err, cont, fit = np.loadtxt(str1, unpack=True)

#Loading the file used for indicating selected data points for the fit
wave_selected, flux_selected, err_selected, cont_selected = np.loadtxt(str2, unpack=True)



#Removing the divide by zero issue that occurs commonly while normalising the data
for i in range(len(wave_selected)):
	if (cont_selected[i]>0.0):
		flux_selected[i] = flux_selected[i]/cont_selected[i]
		err_selected[i] = err_selected[i]/cont_selected[i]
	else:
		flux_selected[i] = np.nan
		err_selected[i] = np.nan


for i in range(len(HI_cont)):
	if HI_cont[i]>0.0:
		HI_flux[i] = HI_flux[i]/HI_cont[i]
		HI_flux_err[i] = HI_flux_err[i]/HI_cont[i]
		HI_prof[i] = HI_prof[i]/HI_cont[i]
	else:
		HI_flux[i] = 0.0
		HI_flux_err[i] = 0.0
		HI_prof[i] = 0.0
		

for i in range(len(H2_cont)):
	if H2_cont[i]>0.0:
		H2_flux[i] = H2_flux[i]/H2_cont[i]
		H2_flux_err[i] = H2_flux_err[i]/H2_cont[i]
		H2_prof[i] = H2_prof[i]/H2_cont[i]
	else:
		H2_flux[i] = 0.0
		H2_flux_err[i] = 0.0
		H2_prof[i] = 0.0
		

for i in range(len(cont)):
	if cont[i]>0.0:
		flux[i] = flux[i]/cont[i]
		err[i] = err[i]/cont[i]
		fit[i] = fit[i]/cont[i]
	else:
		flux[i] = 0.0
		err[i] = 0.0
		fit[i] = 0.0
		





#Plotting the data
f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=False, sharey=False, figsize=(20, 20), dpi=300)


#Setting x-scale limits
ax1.set_xlim(idx1, idx2)
ax2.set_xlim(idx2, idx3)
ax3.set_xlim(idx3, idx4)
ax4.set_xlim(idx4, idx5)


HI_wave2 = HI_wave/(1+z_abs)
H2_wave2 = H2_wave/(1+z_abs)
wave_selected2 = wave_selected/(1+z_abs)


ax1.errorbar(HI_wave2, HI_flux, yerr=HI_flux_err, drawstyle='steps-mid', alpha=0.1, color='k', zorder=6, rasterized=False)
ax1.plot(HI_wave2, HI_prof, 'r-', zorder=7, rasterized=False)
ax1.errorbar(H2_wave2, H2_flux, yerr=H2_flux_err, drawstyle='steps-mid', alpha=0.1, color='k', zorder=6, rasterized=False)
ax1.plot(H2_wave2, H2_prof, 'r-', zorder=7, rasterized=False)
ax1.plot(wave_selected2, flux_selected, 'g.', rasterized=False)


#Setting y-scale limits
ylim_up = 2.2
ax1.set_ylim(-0.2, 2.2)


#ax2.plot(wave, flux, 'b.', alpha=0.8)
ax2.errorbar(HI_wave2, HI_flux, yerr=HI_flux_err, drawstyle='steps-mid', alpha=0.1, color='k', zorder=6, rasterized=False)
#ax2.plot(HI_wave, HI_cont, 'g-')
ax2.plot(HI_wave2, HI_prof, 'r-', zorder=7, rasterized=False)
ax2.errorbar(H2_wave2, H2_flux, yerr=H2_flux_err, drawstyle='steps-mid', alpha=0.1, color='k', zorder=6, rasterized=False)
#ax2.plot(H2_wave, H2_cont, 'g-', rasterized=False)
ax2.plot(H2_wave2, H2_prof, 'r-', zorder=7, rasterized=False)
ax2.plot(wave_selected2, flux_selected, 'g.', rasterized=False)


#Setting y-scale limits
ax2.set_ylim(-0.2, 1.7)


#ax3.plot(wave, flux, 'b.', alpha=0.8)
ax3.errorbar(HI_wave2, HI_flux, yerr=HI_flux_err, drawstyle='steps-mid', alpha=0.1, color='k', zorder=6, rasterized=False)
#ax3.plot(HI_wave, HI_cont, 'g-', zorder=7)
ax3.plot(HI_wave2, HI_prof, 'r-', rasterized=False)
ax3.errorbar(H2_wave2, H2_flux, yerr=H2_flux_err, drawstyle='steps-mid', alpha=0.1, color='k', zorder=6, rasterized=False)
#ax3.plot(H2_wave, H2_cont, 'g-')
ax3.plot(H2_wave2, H2_prof, 'r-', zorder=7, rasterized=False)
ax3.plot(wave_selected2, flux_selected, 'g.', rasterized=False)


#Setting y-scale limits
ax3.set_ylim(-0.2, 1.7)

#'''
ax4.errorbar(HI_wave2, HI_flux, yerr=HI_flux_err, drawstyle='steps-mid', alpha=0.1, color='k', zorder=6, rasterized=False)
ax4.plot(HI_wave, HI_prof, 'r-', zorder=7, rasterized=False)
ax4.errorbar(H2_wave2, H2_flux, yerr=H2_flux_err, drawstyle='steps-mid', alpha=0.1, color='k', zorder=6, rasterized=False)
ax4.plot(H2_wave2, H2_prof, 'r-', zorder=7, rasterized=False)
ax4.plot(wave_selected, flux_selected, 'g.', rasterized=False)


#Setting y-scale limits
ax4.set_ylim(-0.2, 1.2)


ax1.axhline(1.0, color='black', linestyle='dashed', zorder=5, rasterized=False)
ax2.axhline(1.0, color='black', linestyle='dashed', zorder=5, rasterized=False)
ax3.axhline(1.0, color='black', linestyle='dashed', zorder=5, rasterized=False)
ax4.axhline(1.0, color='black', linestyle='dashed', zorder=5, rasterized=False)

ax1.axhline(0.0, color='black', linestyle='dashed', zorder=5, rasterized=False)
ax2.axhline(0.0, color='black', linestyle='dashed', zorder=5, rasterized=False)
ax3.axhline(0.0, color='black', linestyle='dashed', zorder=5, rasterized=False)
ax4.axhline(0.0, color='black', linestyle='dashed', zorder=5, rasterized=False)


#f.subplots_adjust(hspace=0.15)
#f.subplots_adjust(wspace=0)

#plt.tight_layout()
f.subplots_adjust(hspace=0.15)



















#Loading the H2 atomic database file. 
with open ('plotting_files/H2_atomic_file5.csv', "r") as myfile:
	data=myfile.readlines()

name = np.chararray([], itemsize=10)
at_wave = np.chararray([], itemsize=10)
osc = np.chararray([], itemsize=10)
tau = np.chararray([], itemsize=20)
mass_amu = np.chararray([], itemsize=10)
comments1 = np.chararray([], itemsize=30)
band_name = np.chararray([], itemsize=10)
comments2 = np.chararray([], itemsize=10)




for i in range(len(data)):
	name = np.append(name, str(data[i].split()[0]))
	at_wave = np.append(at_wave, str(data[i].split()[1]))
	osc = np.append(osc, str(data[i].split()[2]))
	tau = np.append(tau, str(data[i].split()[3]))
	mass_amu = np.append(mass_amu, str(data[i].split()[4]))
	comments1 = np.append(comments1, str(data[i].split()[5]))
	band_name = np.append(band_name, str(data[i].split()[6]))
	comments2 = np.append(comments2, str(data[i].split()[7]))





name = np.delete(name, 0)
at_wave = np.delete(at_wave, 0)
osc = np.delete(osc, 0)
tau = np.delete(tau, 0)
mass_amu = np.delete(mass_amu, 0)
comments1 = np.delete(comments1, 0)
band_name = np.delete(band_name, 0)
comments2 = np.delete(comments2, 0)







#Creating data for plotting

l0_waves = []
l0_names = []

l1_waves = []
l1_names = []

l2_waves = []
l2_names = []

l3_waves = []
l3_names = []

l4_waves = []
l4_names = []

l5_waves = []
l5_names = []

l6_waves = []
l6_names = []

l7_waves = []
l7_names = []

l8_waves = []
l8_names = []

l9_waves = []
l9_names = []

l10_waves = []
l10_names = []

l11_waves = []
l11_names = []

l12_waves = []
l12_names = []

l13_waves = []
l13_names = []

l14_waves = []
l14_names = []

l15_waves = []
l15_names = []

l16_waves = []
l16_names = []

l17_waves = []
l17_names = []

l18_waves = []
l18_names = []

l19_waves = []
l19_names = []

w0_waves = []
w0_names = []

w1_waves = []
w1_names = []

w2_waves = []
w2_names = []

w3_waves = []
w3_names = []

w4_waves = []
w4_names = []

w5_waves = []
w5_names = []





for i in range(len(band_name)):
	if band_name[i]=='L0':
		l0_names = np.append(l0_names, (name[i][-2:]))
		l0_waves = np.append(l0_waves, (float(at_wave[i])*(1+0.)))

	elif band_name[i]=='L1':
		l1_names = np.append(l1_names, (name[i][-2:]))
		l1_waves = np.append(l1_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L2':
		l2_names = np.append(l2_names, (name[i][-2:]))
		l2_waves = np.append(l2_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L3':
		l3_names = np.append(l3_names, (name[i][-2:]))
		l3_waves = np.append(l3_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L4':
		l4_names = np.append(l4_names, (name[i][-2:]))
		l4_waves = np.append(l4_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L5':
		l5_names = np.append(l5_names, (name[i][-2:]))
		l5_waves = np.append(l5_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L6':
		l6_names = np.append(l6_names, (name[i][-2:]))
		l6_waves = np.append(l6_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L7':
		l7_names = np.append(l7_names, (name[i][-2:]))
		l7_waves = np.append(l7_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L8':
		l8_names = np.append(l8_names, (name[i][-2:]))
		l8_waves = np.append(l8_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L9':
		l9_names = np.append(l9_names, (name[i][-2:]))
		l9_waves = np.append(l9_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L10':
		l10_names = np.append(l10_names, (name[i][-2:]))
		l10_waves = np.append(l10_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L11':
		l11_names = np.append(l11_names, (name[i][-2:]))
		l11_waves = np.append(l11_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L12':
		l12_names = np.append(l12_names, (name[i][-2:]))
		l12_waves = np.append(l12_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L13':
		l13_names = np.append(l13_names, (name[i][-2:]))
		l13_waves = np.append(l13_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L14':
		l14_names = np.append(l14_names, (name[i][-2:]))
		l14_waves = np.append(l14_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L15':
		l15_names = np.append(l15_names, (name[i][-2:]))
		l15_waves = np.append(l15_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L16':
		l16_names = np.append(l16_names, (name[i][-2:]))
		l16_waves = np.append(l16_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L17':
		l17_names = np.append(l17_names, (name[i][-2:]))
		l17_waves = np.append(l17_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L18':
		l18_names = np.append(l18_names, (name[i][-2:]))
		l18_waves = np.append(l18_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='L19':
		l19_names = np.append(l19_names, (name[i][-2:]))
		l19_waves = np.append(l19_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='W0':
		w0_names = np.append(w0_names, (name[i][-2:]))
		w0_waves = np.append(w0_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='W1':
		w1_names = np.append(w1_names, (name[i][-2:]))
		w1_waves = np.append(w1_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='W2':
		w2_names = np.append(w2_names, (name[i][-2:]))
		w2_waves = np.append(w2_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='W3':
		w3_names = np.append(w3_names, (name[i][-2:]))
		w3_waves = np.append(w3_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='W4':
		w4_names = np.append(w4_names, (name[i][-2:]))
		w4_waves = np.append(w4_waves, (float(at_wave[i])*(1+0.)))


	elif band_name[i]=='W5':
		w5_names = np.append(w5_names, (name[i][-2:]))
		w5_waves = np.append(w5_waves, (float(at_wave[i])*(1+0.)))





#Plotting the data


ax1.hlines(y=top_lyman, xmin=l0_waves.min(), xmax=l0_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_lyman, xmin=l0_waves.min(), xmax=l0_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_lyman, xmin=l0_waves.min(), xmax=l0_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_lyman, xmin=l0_waves.min(), xmax=l0_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l0_waves.mean(), top_lyman+0.15, 'L0', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l0_waves.mean(), top_lyman+0.15, 'L0', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l0_waves.mean(), top_lyman+0.15, 'L0', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l0_waves.mean(), top_lyman+0.15, 'L0', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)



ax1.hlines(y=bottom_lyman, xmin=l1_waves.min(), xmax=l1_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=bottom_lyman, xmin=l1_waves.min(), xmax=l1_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=bottom_lyman, xmin=l1_waves.min(), xmax=l1_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=bottom_lyman, xmin=l1_waves.min(), xmax=l1_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l1_waves.mean(), bottom_lyman+0.15, 'L1', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l1_waves.mean(), bottom_lyman+0.15, 'L1', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l1_waves.mean(), bottom_lyman+0.15, 'L1', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l1_waves.mean(), bottom_lyman+0.15, 'L1', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)



ax1.hlines(y=top_lyman, xmin=l2_waves.min(), xmax=l2_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_lyman, xmin=l2_waves.min(), xmax=l2_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_lyman, xmin=l2_waves.min(), xmax=l2_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_lyman, xmin=l2_waves.min(), xmax=l2_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l2_waves.mean(), top_lyman+0.15, 'L2', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l2_waves.mean(), top_lyman+0.15, 'L2', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l2_waves.mean(), top_lyman+0.15, 'L2', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l2_waves.mean(), top_lyman+0.15, 'L2', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=bottom_lyman, xmin=l3_waves.min(), xmax=l3_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=bottom_lyman, xmin=l3_waves.min(), xmax=l3_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=bottom_lyman, xmin=l3_waves.min(), xmax=l3_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=bottom_lyman, xmin=l3_waves.min(), xmax=l3_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l3_waves.mean(), bottom_lyman+0.15, 'L3', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l3_waves.mean(), bottom_lyman+0.15, 'L3', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l3_waves.mean(), bottom_lyman+0.15, 'L3', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l3_waves.mean(), bottom_lyman+0.15, 'L3', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=top_lyman, xmin=l4_waves.min(), xmax=l4_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_lyman, xmin=l4_waves.min(), xmax=l4_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_lyman, xmin=l4_waves.min(), xmax=l4_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_lyman, xmin=l4_waves.min(), xmax=l4_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l4_waves.mean(), top_lyman+0.15, 'L4', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l4_waves.mean(), top_lyman+0.15, 'L4', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l4_waves.mean(), top_lyman+0.15, 'L4', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l4_waves.mean(), top_lyman+0.15, 'L4', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=bottom_lyman, xmin=l5_waves.min(), xmax=l5_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=bottom_lyman, xmin=l5_waves.min(), xmax=l5_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=bottom_lyman, xmin=l5_waves.min(), xmax=l5_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=bottom_lyman, xmin=l5_waves.min(), xmax=l5_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l5_waves.mean(), bottom_lyman+0.15, 'L5', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l5_waves.mean(), bottom_lyman+0.15, 'L5', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l5_waves.mean(), bottom_lyman+0.15, 'L5', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l5_waves.mean(), bottom_lyman+0.15, 'L5', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=top_lyman, xmin=l6_waves.min(), xmax=l6_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_lyman, xmin=l6_waves.min(), xmax=l6_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_lyman, xmin=l6_waves.min(), xmax=l6_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_lyman, xmin=l6_waves.min(), xmax=l6_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l6_waves.mean(), top_lyman+0.15, 'L6', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l6_waves.mean(), top_lyman+0.15, 'L6', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l6_waves.mean(), top_lyman+0.15, 'L6', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l6_waves.mean(), top_lyman+0.15, 'L6', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=bottom_lyman, xmin=l7_waves.min(), xmax=l7_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=bottom_lyman, xmin=l7_waves.min(), xmax=l7_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=bottom_lyman, xmin=l7_waves.min(), xmax=l7_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=bottom_lyman, xmin=l7_waves.min(), xmax=l7_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l7_waves.mean(), bottom_lyman+0.15, 'L7', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l7_waves.mean(), bottom_lyman+0.15, 'L7', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l7_waves.mean(), bottom_lyman+0.15, 'L7', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l7_waves.mean(), bottom_lyman+0.15, 'L7', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=top_lyman, xmin=l8_waves.min(), xmax=l8_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_lyman, xmin=l8_waves.min(), xmax=l8_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_lyman, xmin=l8_waves.min(), xmax=l8_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_lyman, xmin=l8_waves.min(), xmax=l8_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l8_waves.mean(), top_lyman+0.15, 'L8', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l8_waves.mean(), top_lyman+0.15, 'L8', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l8_waves.mean(), top_lyman+0.15, 'L8', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l8_waves.mean(), top_lyman+0.15, 'L8', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=bottom_lyman, xmin=l9_waves.min(), xmax=l9_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=bottom_lyman, xmin=l9_waves.min(), xmax=l9_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=bottom_lyman, xmin=l9_waves.min(), xmax=l9_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=bottom_lyman, xmin=l9_waves.min(), xmax=l9_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l9_waves.mean(), bottom_lyman+0.15, 'L9', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l9_waves.mean(), bottom_lyman+0.15, 'L9', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l9_waves.mean(), bottom_lyman+0.15, 'L9', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l9_waves.mean(), bottom_lyman+0.15, 'L9', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=top_lyman, xmin=l10_waves.min(), xmax=l10_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_lyman, xmin=l10_waves.min(), xmax=l10_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_lyman, xmin=l10_waves.min(), xmax=l10_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_lyman, xmin=l10_waves.min(), xmax=l10_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l10_waves.mean(), top_lyman+0.15, 'L10', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l10_waves.mean(), top_lyman+0.15, 'L10', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l10_waves.mean(), top_lyman+0.15, 'L10', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l10_waves.mean(), top_lyman+0.15, 'L10', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=bottom_lyman, xmin=l11_waves.min(), xmax=l11_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=bottom_lyman, xmin=l11_waves.min(), xmax=l11_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=bottom_lyman, xmin=l11_waves.min(), xmax=l11_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=bottom_lyman, xmin=l11_waves.min(), xmax=l11_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l11_waves.mean(), bottom_lyman+0.15, 'L11', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l11_waves.mean(), bottom_lyman+0.15, 'L11', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l11_waves.mean(), bottom_lyman+0.15, 'L11', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l11_waves.mean(), bottom_lyman+0.15, 'L11', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=top_lyman, xmin=l12_waves.min(), xmax=l12_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_lyman, xmin=l12_waves.min(), xmax=l12_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_lyman, xmin=l12_waves.min(), xmax=l12_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_lyman, xmin=l12_waves.min(), xmax=l12_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l12_waves.mean(), top_lyman+0.15, 'L12', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l12_waves.mean(), top_lyman+0.15, 'L12', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l12_waves.mean(), top_lyman+0.15, 'L12', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l12_waves.mean(), top_lyman+0.15, 'L12', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)



ax1.hlines(y=bottom_lyman, xmin=l13_waves.min(), xmax=l13_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=bottom_lyman, xmin=l13_waves.min(), xmax=l13_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=bottom_lyman, xmin=l13_waves.min(), xmax=l13_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=bottom_lyman, xmin=l13_waves.min(), xmax=l13_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l13_waves.mean(), bottom_lyman+0.15, 'L13', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l13_waves.mean(), bottom_lyman+0.15, 'L13', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l13_waves.mean(), bottom_lyman+0.15, 'L13', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l13_waves.mean(), bottom_lyman+0.15, 'L13', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=top_lyman, xmin=l14_waves.min(), xmax=l14_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_lyman, xmin=l14_waves.min(), xmax=l14_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_lyman, xmin=l14_waves.min(), xmax=l14_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_lyman, xmin=l14_waves.min(), xmax=l14_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l14_waves.mean(), top_lyman+0.15, 'L14', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l14_waves.mean(), top_lyman+0.15, 'L14', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l14_waves.mean(), top_lyman+0.15, 'L14', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l14_waves.mean(), top_lyman+0.15, 'L14', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=bottom_lyman, xmin=l15_waves.min(), xmax=l15_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=bottom_lyman, xmin=l15_waves.min(), xmax=l15_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=bottom_lyman, xmin=l15_waves.min(), xmax=l15_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=bottom_lyman, xmin=l15_waves.min(), xmax=l15_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l15_waves.mean(), bottom_lyman+0.15, 'L15', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l15_waves.mean(), bottom_lyman+0.15, 'L15', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l15_waves.mean(), bottom_lyman+0.15, 'L15', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l15_waves.mean(), bottom_lyman+0.15, 'L15', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=top_lyman, xmin=l16_waves.min(), xmax=l16_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_lyman, xmin=l16_waves.min(), xmax=l16_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_lyman, xmin=l16_waves.min(), xmax=l16_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_lyman, xmin=l16_waves.min(), xmax=l16_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l16_waves.mean(), top_lyman+0.15, 'L16', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l16_waves.mean(), top_lyman+0.15, 'L16', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l16_waves.mean(), top_lyman+0.15, 'L16', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l16_waves.mean(), top_lyman+0.15, 'L16', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=bottom_lyman, xmin=l17_waves.min(), xmax=l17_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=bottom_lyman, xmin=l17_waves.min(), xmax=l17_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=bottom_lyman, xmin=l17_waves.min(), xmax=l17_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=bottom_lyman, xmin=l17_waves.min(), xmax=l17_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l17_waves.mean(), bottom_lyman+0.15, 'L17', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l17_waves.mean(), bottom_lyman+0.15, 'L17', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l17_waves.mean(), bottom_lyman+0.15, 'L17', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l17_waves.mean(), bottom_lyman+0.15, 'L17', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=top_lyman, xmin=l18_waves.min(), xmax=l18_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_lyman, xmin=l18_waves.min(), xmax=l18_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_lyman, xmin=l18_waves.min(), xmax=l18_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_lyman, xmin=l18_waves.min(), xmax=l18_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l18_waves.mean(), top_lyman+0.15, 'L18', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l18_waves.mean(), top_lyman+0.15, 'L18', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l18_waves.mean(), top_lyman+0.15, 'L18', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l18_waves.mean(), top_lyman+0.15, 'L18', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=bottom_lyman, xmin=l19_waves.min(), xmax=l19_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=bottom_lyman, xmin=l19_waves.min(), xmax=l19_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=bottom_lyman, xmin=l19_waves.min(), xmax=l19_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=bottom_lyman, xmin=l19_waves.min(), xmax=l19_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(l19_waves.mean(), bottom_lyman+0.15, 'L19', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(l19_waves.mean(), bottom_lyman+0.15, 'L19', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(l19_waves.mean(), bottom_lyman+0.15, 'L19', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(l19_waves.mean(), bottom_lyman+0.15, 'L19', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)




ax1.hlines(y=top_werner, xmin=w0_waves.min(), xmax=w0_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_werner, xmin=w0_waves.min(), xmax=w0_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_werner, xmin=w0_waves.min(), xmax=w0_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_werner, xmin=w0_waves.min(), xmax=w0_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(w0_waves.mean(), top_werner+0.15, 'W0', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(w0_waves.mean(), top_werner+0.15, 'W0', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(w0_waves.mean(), top_werner+0.15, 'W0', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(w0_waves.mean(), top_werner+0.15, 'W0', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)



ax1.hlines(y=top_werner, xmin=w1_waves.min(), xmax=w1_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_werner, xmin=w1_waves.min(), xmax=w1_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_werner, xmin=w1_waves.min(), xmax=w1_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_werner, xmin=w1_waves.min(), xmax=w1_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(w1_waves.mean(), top_werner+0.15, 'W1', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(w1_waves.mean(), top_werner+0.15, 'W1', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(w1_waves.mean(), top_werner+0.15, 'W1', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(w1_waves.mean(), top_werner+0.15, 'W1', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=top_werner, xmin=w2_waves.min(), xmax=w2_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_werner, xmin=w2_waves.min(), xmax=w2_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_werner, xmin=w2_waves.min(), xmax=w2_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_werner, xmin=w2_waves.min(), xmax=w2_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(w2_waves.mean(), top_werner+0.15, 'W2', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(w2_waves.mean(), top_werner+0.15, 'W2', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(w2_waves.mean(), top_werner+0.15, 'W2', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(w2_waves.mean(), top_werner+0.15, 'W2', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=top_werner, xmin=w3_waves.min(), xmax=w3_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_werner, xmin=w3_waves.min(), xmax=w3_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_werner, xmin=w3_waves.min(), xmax=w3_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_werner, xmin=w3_waves.min(), xmax=w3_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(w3_waves.mean(), top_werner+0.15, 'W3', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(w3_waves.mean(), top_werner+0.15, 'W3', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(w3_waves.mean(), top_werner+0.15, 'W3', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(w3_waves.mean(), top_werner+0.15, 'W3', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)



ax1.hlines(y=top_werner, xmin=w4_waves.min(), xmax=w4_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_werner, xmin=w4_waves.min(), xmax=w4_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_werner, xmin=w4_waves.min(), xmax=w4_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_werner, xmin=w4_waves.min(), xmax=w4_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(w4_waves.mean(), top_werner+0.15, 'W4', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(w4_waves.mean(), top_werner+0.15, 'W4', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(w4_waves.mean(), top_werner+0.15, 'W4', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(w4_waves.mean(), top_werner+0.15, 'W4', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)


ax1.hlines(y=top_werner, xmin=w5_waves.min(), xmax=w5_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax2.hlines(y=top_werner, xmin=w5_waves.min(), xmax=w5_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax3.hlines(y=top_werner, xmin=w5_waves.min(), xmax=w5_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)
ax4.hlines(y=top_werner, xmin=w5_waves.min(), xmax=w5_waves.max(), linewidth=2, color='b', zorder=10, rasterized=False)

ax1.text(w5_waves.mean(), top_werner+0.15, 'W5', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax2.text(w5_waves.mean(), top_werner+0.15, 'W5', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax3.text(w5_waves.mean(), top_werner+0.15, 'W5', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)
ax4.text(w5_waves.mean(), top_werner+0.15, 'W5', color='blue', fontsize=0.5*size_of_font,  zorder=10, clip_on=True, rasterized=False)



for i in range(len(l0_waves)):
	ax1.plot(l0_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l0_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l0_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l0_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l0_waves[i]-1, even_lyman_text, l0_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l0_waves[i]-1, even_lyman_text, l0_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)

	if (i%2==0.):	
		ax3.text(l0_waves[i]-1, even_lyman_text, l0_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	else:	
		ax3.text(l0_waves[i]+1, even_lyman_text, l0_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)

	ax4.text(l0_waves[i]-1, even_lyman_text, l0_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


for i in range(len(l1_waves)):
	ax1.plot(l1_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l1_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l1_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l1_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l1_waves[i]-1, odd_lyman_text, l1_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l1_waves[i]-1, odd_lyman_text, l1_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	#ax3.text(l1_waves[i]-1, odd_lyman_text, l1_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	if (i%2==0.):	
		ax3.text(l1_waves[i]-1, odd_lyman_text, l1_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	else:	
		ax3.text(l1_waves[i]+1, odd_lyman_text, l1_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


	ax4.text(l1_waves[i]-1, odd_lyman_text, l1_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


for i in range(len(l2_waves)):
	ax1.plot(l2_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l2_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l2_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l2_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l2_waves[i]-1, even_lyman_text, l2_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l2_waves[i]-1, even_lyman_text, l2_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l2_waves[i]-1, even_lyman_text, l2_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l2_waves[i]-1, even_lyman_text, l2_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


for i in range(len(l3_waves)):
	ax1.plot(l3_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l3_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l3_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l3_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l3_waves[i]-1, odd_lyman_text, l3_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l3_waves[i]-1, odd_lyman_text, l3_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l3_waves[i]-1, odd_lyman_text, l3_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l3_waves[i]-1, odd_lyman_text, l3_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


for i in range(len(l4_waves)):
	ax1.plot(l4_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l4_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l4_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l4_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l4_waves[i]-1, even_lyman_text, l4_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l4_waves[i]-1, even_lyman_text, l4_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l4_waves[i]-1, even_lyman_text, l4_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l4_waves[i]-1, even_lyman_text, l4_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


for i in range(len(l5_waves)):
	ax1.plot(l5_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l5_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l5_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l5_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l5_waves[i]-1, odd_lyman_text, l5_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l5_waves[i]-1, odd_lyman_text, l5_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l5_waves[i]-1, odd_lyman_text, l5_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l5_waves[i]-1, odd_lyman_text, l5_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


for i in range(len(l6_waves)):
	ax1.plot(l6_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l6_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l6_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l6_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l6_waves[i]-1, even_lyman_text, l6_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l6_waves[i]-1, even_lyman_text, l6_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l6_waves[i]-1, even_lyman_text, l6_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l6_waves[i]-1, even_lyman_text, l6_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


for i in range(len(l7_waves)):
	ax1.plot(l7_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l7_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l7_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l7_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l7_waves[i]-1, odd_lyman_text, l7_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l7_waves[i]-1, odd_lyman_text, l7_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l7_waves[i]-1, odd_lyman_text, l7_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l7_waves[i]-1, odd_lyman_text, l7_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)



for i in range(len(l8_waves)):
	ax1.plot(l8_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l8_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l8_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l8_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l8_waves[i]-1, even_lyman_text, l8_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l8_waves[i]-1, even_lyman_text, l8_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l8_waves[i]-1, even_lyman_text, l8_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l8_waves[i]-1, even_lyman_text, l8_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


for i in range(len(l9_waves)):
	ax1.plot(l9_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l9_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l9_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l9_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l9_waves[i]-1, odd_lyman_text, l9_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l9_waves[i]-1, odd_lyman_text, l9_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l9_waves[i]-1, odd_lyman_text, l9_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l9_waves[i]-1, odd_lyman_text, l9_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)



for i in range(len(l10_waves)):
	ax1.plot(l10_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l10_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l10_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l10_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l10_waves[i]-1, even_lyman_text, l10_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l10_waves[i]-1, even_lyman_text, l10_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l10_waves[i]-1, even_lyman_text, l10_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l10_waves[i]-1, even_lyman_text, l10_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)



for i in range(len(l11_waves)):
	ax1.plot(l11_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l11_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l11_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l11_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l11_waves[i]-1, odd_lyman_text, l11_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l11_waves[i]-1, odd_lyman_text, l11_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l11_waves[i]-1, odd_lyman_text, l11_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l11_waves[i]-1, odd_lyman_text, l11_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)



for i in range(len(l12_waves)):
	ax1.plot(l12_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l12_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l12_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l12_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l12_waves[i]-1, even_lyman_text, l12_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l12_waves[i]-1, even_lyman_text, l12_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l12_waves[i]-1, even_lyman_text, l12_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l12_waves[i]-1, even_lyman_text, l12_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)



for i in range(len(l13_waves)):
	ax1.plot(l13_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l13_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l13_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l13_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l13_waves[i]-1, odd_lyman_text, l13_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l13_waves[i]-1, odd_lyman_text, l13_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l13_waves[i]-1, odd_lyman_text, l13_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l13_waves[i]-1, odd_lyman_text, l13_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)



for i in range(len(l14_waves)):
	ax1.plot(l14_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l14_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l14_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l14_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l14_waves[i]-1, even_lyman_text, l14_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l14_waves[i]-1, even_lyman_text, l14_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l14_waves[i]-1, even_lyman_text, l14_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l14_waves[i]-1, even_lyman_text, l14_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)



for i in range(len(l15_waves)):
	ax1.plot(l15_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l15_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l15_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l15_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l15_waves[i]-1, odd_lyman_text, l15_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l15_waves[i]-1, odd_lyman_text, l15_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l15_waves[i]-1, odd_lyman_text, l15_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l15_waves[i]-1, odd_lyman_text, l15_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)



for i in range(len(l16_waves)):
	ax1.plot(l16_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l16_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l16_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l16_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l16_waves[i]-1, even_lyman_text, l16_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l16_waves[i]-1, even_lyman_text, l16_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l16_waves[i]-1, even_lyman_text, l16_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l16_waves[i]-1, even_lyman_text, l16_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)



for i in range(len(l17_waves)):
	ax1.plot(l17_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l17_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l17_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l17_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l17_waves[i]-1, odd_lyman_text, l17_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l17_waves[i]-1, odd_lyman_text, l17_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l17_waves[i]-1, odd_lyman_text, l17_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l17_waves[i]-1, odd_lyman_text, l17_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)




for i in range(len(l18_waves)):
	ax1.plot(l18_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l18_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l18_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l18_waves[i], even_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l18_waves[i]-1, even_lyman_text, l18_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l18_waves[i]-1, even_lyman_text, l18_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l18_waves[i]-1, even_lyman_text, l18_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l18_waves[i]-1, even_lyman_text, l18_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)



for i in range(len(l19_waves)):
	ax1.plot(l19_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(l19_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(l19_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(l19_waves[i], odd_lyman_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(l19_waves[i]-1, odd_lyman_text, l19_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(l19_waves[i]-1, odd_lyman_text, l19_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(l19_waves[i]-1, odd_lyman_text, l19_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(l19_waves[i]-1, odd_lyman_text, l19_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


#'''
for i in range(len(w0_waves)):
	ax1.plot(w0_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(w0_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(w0_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(w0_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(w0_waves[i]-1, werner_text, w0_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(w0_waves[i]-1, werner_text, w0_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(w0_waves[i]-1, werner_text, w0_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(w0_waves[i]-1, werner_text, w0_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


for i in range(len(w1_waves)):
	ax1.plot(w1_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(w1_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(w1_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(w1_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(w1_waves[i]-1, werner_text, w1_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(w1_waves[i]-1, werner_text, w1_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(w1_waves[i]-1, werner_text, w1_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(w1_waves[i]-1, werner_text, w1_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


for i in range(len(w2_waves)):
	ax1.plot(w2_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(w2_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(w2_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(w2_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(w2_waves[i]-1, werner_text, w2_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(w2_waves[i]-1, werner_text, w2_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(w2_waves[i]-1, werner_text, w2_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(w2_waves[i]-1, werner_text, w2_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)



for i in range(len(w3_waves)):
	ax1.plot(w3_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(w3_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(w3_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(w3_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(w3_waves[i]-1, werner_text, w3_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(w3_waves[i]-1, werner_text, w3_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(w3_waves[i]-1, werner_text, w3_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(w3_waves[i]-1, werner_text, w3_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


for i in range(len(w4_waves)):
	ax1.plot(w4_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(w4_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(w4_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(w4_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(w4_waves[i]-1, werner_text, w4_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(w4_waves[i]-1, werner_text, w4_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(w4_waves[i]-1, werner_text, w4_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(w4_waves[i]-1, werner_text, w4_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)



for i in range(len(w5_waves)):
	ax1.plot(w5_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax2.plot(w5_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax3.plot(w5_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax4.plot(w5_waves[i], werner_arrow, linestyle='none', marker=r'$\downarrow$', markersize=10, rasterized=False)
	ax1.text(w5_waves[i]-1, werner_text, w5_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax2.text(w5_waves[i]-1, werner_text, w5_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax3.text(w5_waves[i]-1, werner_text, w5_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)
	ax4.text(w5_waves[i]-1, werner_text, w5_names[i][1:], fontsize=0.5*size_of_font, clip_on=True, rasterized=False)


#'''





#Setting tick paramters

ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

ax3.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax3.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

ax4.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax4.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))


ax1.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(2))

ax2.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(2))

ax3.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax3.xaxis.set_minor_locator(ticker.MultipleLocator(2))

ax4.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax4.xaxis.set_minor_locator(ticker.MultipleLocator(10))

ax1.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')
ax2.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')
ax3.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')
ax4.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')

ax1.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')
ax2.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')
ax3.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')
ax4.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')






#Miscellneous text for the plot

f.text(0.5, 0.92, (str('QSO SDSS ') + str(name_DLA)), ha='center', va='center', fontsize=1.5*size_of_font, rasterized=False)
#plt.xlabel('Observed wavelength ($\AA$)', fontsize=size_of_font)
#plt.ylabel('Normalised flux', fontsize=size_of_font)
f.text(0.5, 0.04, r'Rest wavelength ($\AA$)', ha='center', fontsize=1.5*size_of_font, rasterized=False)
f.text(0.04, 0.5, 'Normalised flux', va='center', rotation='vertical', fontsize=1.5*size_of_font, rasterized=False)

ax1.tick_params(axis='both', labelsize=size_of_font)
ax2.tick_params(axis='both', labelsize=size_of_font)
ax3.tick_params(axis='both', labelsize=size_of_font)
ax4.tick_params(axis='both', labelsize=size_of_font)





#Saving the plot in a file
file_name = name_DLA + '_HI_H2_rest_wave.pdf'

plt.savefig(file_name)
#plt.title(name_DLA, fontsize=size_of_font)
#plt.show()



