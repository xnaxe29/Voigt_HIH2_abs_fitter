# Voigt_HIH2_abs_fitter
GUI code to fit Lyman-series absorption bands of HI and Lyman-Werner absorption bands of H2 and HD

This is a GUI code that takes in a spectra from a source file and estimate any absorption feature from HI, H2 and other ions as well. However, the GUI features are currently available only for HI and H2 rest frame UV absorption. The code further lets you draw a continuum by hand for an initial guess on the column densitities of HI and rotational levels of H2/HD. 

Required dependancies - 
This is a simple python code (version - 2.7). Due to issues in running matplotlib widges in version 3, unfortuantely the code might not work in python 3. However, experts are welcome to change the code and try. Here are the required python modules for your reference - 


numpy
scipy
matplotlib
pyastronomy
sys
from pathlib import Path
os
pyfits
csv
itertools
tabulate

Apart from this, the primary GUI code also imports a custom file - 'relevant_functions.py', which is given in the main folder of this repository. Other files are useful for making plots, but are not essential for running the code itself.



Other essential non python files include - 

1. atomic_database_file.dat - This is the file from which the atomic data of the relevant ions will be picked up from. Any other file with the same format can be used to replace this.

2. initial_parameters.dat - This is the primary parameter file and is required to give initial values for many paramters used by the GUI. The file is self descriptive. Please have a look to change the default.

3. J2140-0321_uvb.spec - This is an example of an acceptable format of the source file. It has to be a binary with wavelength, flux and flux error in that order. Any comments should start with '#'.

4. quasar_template_selsing_2015.dat - This is an example of the optical quasar template spectra that I have used for comparing and making intial guesses on the emission line transitions as well as continuum bluer than lyman-alpha emission. Any other similar file can be used. The name of the template file should be updated in the GUI code, if changed.

5. relevant_functions.py - This is a file with all relevant custom functions defined which is utilised in the GUI code.



Non essential files

6. relevant_functions.pyc - This is the compiled bytecode of 'relevant_functions.py' and can be deleted. This will be generated everytime the main GUI code (main_GUI.py) will be run.

7. Files inside 'plotting_files' folder - Inside this folder, there are plotting routines (.py codes) and files essential to run those routines. 

8. Files inside 'plotting_files' folder - This is a dedicated folder for plotting routine (.py code) and relevant files, to plot the H2 excitation diagram for the system. 


Some common warnings/errors one might encounter - 

1. 'x0' is infeasible - One might get this error if they are trying to fit the HI/H2 lines with estimated column density of the rotational levels, logN[atoms cm^-2] < 9.0. Please try to keep the values of all column densities, logN[atoms cm^-2] > 9.0

2







Here is a short description of interactive functions that one can perform and all the buttons that you see in the GUI  -   

1. Drawing a continnum by hand - Hover your mouse over the plot and press 'a' to add a continuum points. For modifying the same point, hover the mouse in close proximity to that point and press 'm'. For removing a point, press 'r'. I recommend smoothing the data a bit before drawing the continuum. The smoothed data helps a lot in  



