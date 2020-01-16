import matplotlib.pyplot as plt
import numpy as np
import sys
from tabulate import tabulate


b_val_plot = np.array([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.])


file_name_b_val = str(sys.argv[1])
file_name_b_val_err = str(sys.argv[2])



#print (file_name_b_val)
#print (file_name_b_val_err)
#quit()

b_val_new_array = np.loadtxt(file_name_b_val)
b_val_new_array_err = np.loadtxt(file_name_b_val_err)

table_array = np.chararray([12,13], itemsize=35)
table_array[1:,:] = np.round(b_val_new_array[:], 2)
table_array[0,:] = ['N(HI)', 'N(H2J0)', 'N(H2J1)', 'N(H2J2)', 'N(H2J3)', 'N(H2J4)', 'N(H2J5)', 'N(H2J6)', 'N(H2J7)', 'N(HDJ0)', 'N(HDJ1)', 'N(HDJ2)', 'b_val']
table_array_new = np.transpose(table_array)

print ('\n\n\n')
print(tabulate(table_array_new, tablefmt="latex", floatfmt="2.2f"))
print ('\n\n\n')


table_array_err = np.chararray([12,13], itemsize=35)
table_array_err[1:,:] = np.round(b_val_new_array_err[:], 2)
table_array_err[0,:] = ['N(HI)_err', 'N(H2J0)_err', 'N(H2J1)_err', 'N(H2J2)_err', 'N(H2J3)_err', 'N(H2J4)_err', 'N(H2J5)_err', 'N(H2J6)_err', 'N(H2J7)_err', 'N(HDJ0)_err', 'N(HDJ1)_err', 'N(HDJ2)_err', 'b_val_err']
table_array_new_err = np.transpose(table_array_err)

print ('\\vspace{20pt}')
print ('\n\n\n')
print(tabulate(table_array_new_err, tablefmt="latex", floatfmt="2.2f"))
print ('\n\n\n')



quit()











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


f, ((ax1, ax2, ax3, ax4), (ax6, ax7, ax8, ax9), (ax10, ax11, ax12, ax13)) = plt.subplots(3, 4, sharex='col', sharey='row', figsize=(9, 12), dpi=300)
#f, axarr = plt.subplots(5, 5, sharex='col', sharey='row')

ax1.errorbar(b_val_plot, logN_H2J0_fit_array3, yerr=logN_H2J0_fit_err_array3)
ax1.set_title('H2J0')
ax1.set_ylabel('logN')


ax2.errorbar(b_val_plot, logN_H2J1_fit_array3, yerr=logN_H2J1_fit_err_array3)
ax2.set_title('H2J1')

ax3.errorbar(b_val_plot, logN_H2J2_fit_array3, yerr=logN_H2J2_fit_err_array3)
ax3.set_title('H2J2')

ax4.errorbar(b_val_plot, logN_H2J3_fit_array3, yerr=logN_H2J3_fit_err_array3)
ax4.set_title('H2J3')

ax6.errorbar(b_val_plot, logN_H2J4_fit_array3, yerr=logN_H2J4_fit_err_array3)
ax6.set_title('H2J4')
ax6.set_ylabel('logN')

ax7.errorbar(b_val_plot, logN_H2J5_fit_array3, yerr=logN_H2J5_fit_err_array3)
ax7.set_title('H2J5')

ax8.errorbar(b_val_plot, logN_H2J6_fit_array3, yerr=logN_H2J6_fit_err_array3)
ax8.set_title('H2J6')

ax9.errorbar(b_val_plot, logN_H2J7_fit_array3, yerr=logN_H2J7_fit_err_array3)
ax9.set_title('H2J7')

ax10.errorbar(b_val_plot, logN_HDJ0_fit_array3, yerr=logN_HDJ0_fit_err_array3)
ax10.set_title('HDJ0')
ax10.set_ylabel('logN')
ax10.set_xlabel('b value (Km/s)')

ax11.errorbar(b_val_plot, logN_HDJ1_fit_array3, yerr=logN_HDJ1_fit_err_array3)
ax11.set_title('HDJ1')
ax11.set_xlabel('b value (Km/s)')
ax12.set_xlabel('b value (Km/s)')
ax13.set_xlabel('b value (Km/s)')

#ax12.errorbar(b_val_plot, logN_HDJ2_fit_array3, yerr=logN_HDJ2_fit_err_array3)
#ax12.set_title('HDJ2')

#ax13.errorbar(b_val_plot, logN_HDJ2_fit_array3, yerr=logN_HDJ2_fit_err_array3)
#ax13.set_title('HDJ2')

#print (logN_HDJ2_fit_err_array3)
#print (logN_HDJ2_fit_err_array3)


#plt.show()

saved_filename = file_name_b_val[:-11] + 'vs_col_den.pdf'
#print (saved_filename)
#quit()

ax1.margins(y=.1, x=.1)
ax2.margins(y=.1, x=.1)
ax3.margins(y=.1, x=.1)
ax4.margins(y=.1, x=.1)
ax6.margins(y=.1, x=.1)
ax7.margins(y=.1, x=.1)
ax8.margins(y=.1, x=.1)
ax9.margins(y=.1, x=.1)
ax10.margins(y=.1, x=.1)
ax11.margins(y=.1, x=.1)
ax12.margins(y=.1, x=.1)
ax13.margins(y=.1, x=.1)

plt.savefig(saved_filename)
print ('saved')


