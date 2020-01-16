import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.ticker as ticker

b_val_plot = np.array([1., 2.,3.,4.,5.,6.,7.,8.,9.,10.])


size_of_font = 20

file_name_b_val = str(sys.argv[1])
file_name_b_val_err = str(sys.argv[2])



#print (file_name_b_val)
#print (file_name_b_val_err)
#quit()

b_val_new_array = np.loadtxt(file_name_b_val)
b_val_new_array_err = np.loadtxt(file_name_b_val_err)


b_val_new_array = np.nan_to_num(b_val_new_array)
b_val_new_array_err = np.nan_to_num(b_val_new_array_err)

#print (b_val_new_array)
#print (b_val_new_array_err)
#quit()

logN_HI_fit_array3 = b_val_new_array[:-1,0]
logN_HI_fit_array_err3 = b_val_new_array_err[:-1,0]
logN_H2J0_fit_array3 = b_val_new_array[:-1,1]
logN_H2J0_fit_err_array3 = b_val_new_array_err[:-1,1]
logN_H2J1_fit_array3 = b_val_new_array[:-1,2]
logN_H2J1_fit_err_array3 = b_val_new_array_err[:-1,2]
logN_H2J2_fit_array3 = b_val_new_array[:-1,3]
logN_H2J2_fit_err_array3 = b_val_new_array_err[:-1,3]
logN_H2J3_fit_array3 = b_val_new_array[:-1,4]
logN_H2J3_fit_err_array3 = b_val_new_array_err[:-1,4]
logN_H2J4_fit_array3 = b_val_new_array[:-1,5]
logN_H2J4_fit_err_array3 = b_val_new_array_err[:-1,5]
logN_H2J5_fit_array3 = b_val_new_array[:-1,6]
logN_H2J5_fit_err_array3 = b_val_new_array_err[:-1,6]
logN_H2J6_fit_array3 = b_val_new_array[:-1,7]
logN_H2J6_fit_err_array3 = b_val_new_array_err[:-1,7]
logN_H2J7_fit_array3 = b_val_new_array[:-1,8]
logN_H2J7_fit_err_array3 = b_val_new_array_err[:-1,8]
logN_HDJ0_fit_array3 = b_val_new_array[:-1,9]
logN_HDJ0_fit_err_array3 = b_val_new_array_err[:-1,9]
logN_HDJ1_fit_array3 = b_val_new_array[:-1,10]
logN_HDJ1_fit_err_array3 = b_val_new_array_err[:-1,10]
logN_HDJ2_fit_array3 = b_val_new_array[:-1,11]
logN_HDJ2_fit_err_array3 = b_val_new_array_err[:-1,11]
b_fit_array3 = b_val_new_array[:-1,12]
b_fit_err_array3 = b_val_new_array_err[:-1,12]


f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, sharex='col', sharey='row', figsize=(16, 8), dpi=300)
#f, axarr = plt.subplots(5, 5, sharex='col', sharey='row')

#print (logN_H2J0_fit_array3)

ax1.errorbar(b_val_plot, logN_H2J0_fit_array3, yerr=logN_H2J0_fit_err_array3)
ax1.set_title('H2J0', fontsize=1.5*size_of_font)
ax1.set_ylabel('logN', fontsize=size_of_font)

ax1.set_ylim(logN_H2J0_fit_array3.min()-5, logN_H2J0_fit_array3.max()+5)
#print (logN_H2J0_fit_array3)

ax2.errorbar(b_val_plot, logN_H2J1_fit_array3, yerr=logN_H2J1_fit_err_array3)
ax2.set_title('H2J1', fontsize=1.5*size_of_font)

ax3.errorbar(b_val_plot, logN_H2J2_fit_array3, yerr=logN_H2J2_fit_err_array3)
ax3.set_title('H2J2', fontsize=1.5*size_of_font)

ax4.errorbar(b_val_plot, logN_H2J3_fit_array3, yerr=logN_H2J3_fit_err_array3)
ax4.set_title('H2J3', fontsize=1.5*size_of_font)

ax5.errorbar(b_val_plot, logN_H2J4_fit_array3, yerr=logN_H2J4_fit_err_array3)
ax5.set_title('H2J4', fontsize=1.5*size_of_font)
ax5.set_ylabel('logN', fontsize=size_of_font)


ax5.set_ylim(logN_H2J4_fit_array3.min()-5, logN_H2J4_fit_array3.max()+5)

ax6.errorbar(b_val_plot, logN_H2J5_fit_array3, yerr=logN_H2J5_fit_err_array3)
ax6.set_title('H2J5', fontsize=1.5*size_of_font)

ax7.errorbar(b_val_plot, logN_H2J6_fit_array3, yerr=logN_H2J6_fit_err_array3)
ax7.set_title('H2J6', fontsize=1.5*size_of_font)
ax5.set_xlabel(r'b value (Km s$^{-1}$)', fontsize=size_of_font)
ax6.set_xlabel(r'b value (Km s$^{-1}$)', fontsize=size_of_font)
ax7.set_xlabel(r'b value (Km s$^{-1}$)', fontsize=size_of_font)
ax8.set_xlabel(r'b value (Km s$^{-1}$)', fontsize=size_of_font)



saved_filename = file_name_b_val[:-11] + 'vs_col_den.pdf'


ax1.margins(y=.1, x=.1)
ax2.margins(y=.1, x=.1)
ax3.margins(y=.1, x=.1)
ax4.margins(y=.1, x=.1)
ax5.margins(y=.1, x=.1)
ax6.margins(y=.1, x=.1)
ax7.margins(y=.1, x=.1)
ax8.margins(y=.1, x=.1)

'''

ax1.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

ax2.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
ax2.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

ax3.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax3.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
ax3.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax3.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

ax4.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax4.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
ax4.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax4.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

ax5.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax5.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax5.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax5.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

ax6.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax6.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax6.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax6.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

ax7.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax7.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax7.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax7.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

ax8.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax8.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax8.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax8.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))

'''

ax1.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')
ax2.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')
ax3.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')
ax4.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')
ax5.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')
ax6.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')
ax7.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')
ax8.tick_params(axis = 'both', which = 'major', direction='in', length=10, width=2, colors='k')

ax1.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')
ax2.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')
ax3.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')
ax4.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')
ax5.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')
ax6.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')
ax7.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')
ax8.tick_params(axis = 'both', which = 'minor', direction='in', length=5, width=1, colors='k')




plt.savefig(saved_filename)
print ('saved')
quit()

