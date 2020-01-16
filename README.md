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

Other import files include - 

1. atomic_database_file.dat - This is the file from which the atomic data of the relevant ions will be picked up from. 


initial_guess.dat




Here is a short description of all the buttons that you see in the GUI and other interactive functions that you can perform with keyboard and mouse.  
