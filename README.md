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

pathlib

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

2. "File "main_GUI.py", line 232, in on_key
    x1, y1 = zip(*points)
ValueError: need more than 0 values to unpack" - 
This is a common warning indicating that there are not enough spline points to create a continuum. Please ignore this warning and keep adding continuum points. Once 3-4 points are added, this warning will vanish and a continuum will be drawn.

3. 






Here is a short description of interactive functions that one can perform and all the buttons that you see in the GUI  -   

1. Drawing a continnum by hand - Hover your mouse over the plot and press 'a' to add a continuum points. For modifying the same point, hover the mouse in close proximity to that point and press 'm'. For removing a point, press 'r'. I recommend smoothing the data a bit before drawing the continuum. The smoothed data helps a lot in identifying the correct continuum.

2. Selecting data points for fitting purpose - The absorption lines of HI (Lyman-series) and H2 (Lyman-Werner) are heavily contaminated by the lyman-alpha forest lines. Apart from this, there is also absorption from metal lines from the same redshift as that of a strong absorber (say, e.g. a Damped Lyman-alpha system). For the purpose of fitting the HI/H2 lines, all these can be treated as contaminants. Hence, for fitting, I have created a point selection method. One has to specify data points that are free of absorption line contaminants before one starts the fit. To do the same, one has to press the buttons 't' -> 'y' -> 'h' in that sequence. This activates the selection function and now one can draw a rectangle around data points to select them for fitting. In case you would like to remove some data points that have been selected by mistake, you can activate the remove function by pressing 'j' (or 't' -> 'y' -> 'j', if the function is inactive). After this, the points inside any selected rectangle will disappear.

3. All the matplotlib basic plotting buttons and keys work as their defaults.

All custom button functionalities (clockwise) - 

4. Bar: 'Flux_reduction' - The quasar template has a specific continuum flux that may or may not be on the same scale as that of the quasar being analysed. This bar is used to scale the continuum up and down to match the flux of the quasar. We have used the template of a quasar for this code. However, the template can easily be changed to other type of bright source.

5. Bar: 'Resolution' - This bar can be used to vary the spectral resolution of the observation in real time. It can be used when the precise resolution of the observation is not known. One can play around with this to see the changes in absorption lines to estimate the correct resolution.  

6. Radio buttons: 'Vary/Fix HI Cont' - In all cases of strong absorption (also seen in the example), the quasar continuum is completely absorbed by the HI line at the absorber rest frame. In such cases, it is best to let the code fit the continuum along with the absorption profile. This code uses chebyshev polynomials of certain order (which is decide again on different parameters) to model a continuum over HI 1215 absorption line. Please note that for this to work properly, you also have to specify the region of HI fitting in the 'initial_parameters.dat' file. The continuum fitting will only occur on that specified region.

7. Radio buttons: 'Vary/Fix H2 Cont' - Similar to above, this is used to specify if you want the code to fit a continuum for two H2 regions. The area of these H2 regions also needs to be specified in the 'initial_parameters.dat' file.

8. Radio buttons: 'Vary/Fix H2 Redshift' - This button is used to specify if you would like to keep the H2 redshift fixed or would like the fit the redshift along with the column densities. HI redshift is always decided by the fit.

9. Radio buttons: 'Vary/Fix H2J0/H2J1' - This button is used to specify whether you would like to vary the column densities of intial rotational levels of H2 (J=0,1, and 2) during the fit. 

10. Radio buttons: 'Vary/Fix H2J2+' - This button is used to specify whether you would like to vary the column densities of other rotational levels of H2 (J=3 to J=7) during the fit.

11. Radio buttons: 'Vary/Fix HD' - This button is used to specify whether you would like to vary the column densities of all rotational levels of HD (J=0,1 and 2) during the fit.

12. Button: 'Fit' - The button is for executing the fitting function.

13. Button: 'save_fit' - The button is used for saving the fitted parameters.

14. Button: 'save_params' - The button is used for all the current parameters displayed in the GUI regardless of whether fitting has occured.

15. Button: 'get_params' - Load the last saved parameters.

16. Button: 'save_cont' - Save the drawn or fitted continuum

17. Button: 'load_cont' - Load the drawn or fitted continuum
  
18. Button: 'save_data' - Save the regions that have been selected to perform a fit.

19. Button: 'load_data' - Load the regions that have been selected to perform a fit.

20. Bar: 'HI' - This is the bar that can be used to vary the HI column density.

21. Bars: 'H2J0' to 'H2J7' - These are 8 bars that can be used to vary the column densities of individual rotational levels of H2. 

21. Bars: 'HDJ0' to 'HDJ2' - These are 3 bars that can be used to vary the column densities of individual rotational levels of HD.

22. Button: 'fit_b_val' - This is a special function that is created to fit the column densities of H2 with varying doppler parameter from 1km/s to 10km/s (This range can be changed in the code). This (and other four buttons, detailed from 23 to 26) are dedicated to find the fit in case where the spectral resolution is not good enough to estimate the correct doppler parameter for H2.  

23. Button: 'save_b_val' - Saves the fit performed above.

24. Button: 'get_b_val' - Loads the fit performed above.

25. Button: 'relevant_figs' - Creates the H2 excitation diagram and H2 rotational level velocity plots for the mean column density obtained from the fit (with errors displaying the variance in column density). 

26. Button: 'b_val_table' - Creates a latex script indicating the details of the obtained column density for the fit.

27. Button: 'plot_total_fit' - Plots the H2 and HI lines in 4 subplots indication each rotational level of H2 in different Lyman and Werner Bands.

28. Button: 'plot_H2J0_H2J1' - Plots in velocity different rotational levels of H2. The transitions for which the plots will be created can be given in the file - 'plotting_files/list_H2J0_H2J1_custom.txt'.

29. Button: 'plot_HI' - Plots the HI 1215 line in velocity space. The code can also be altered to plot all Lyman series lines in velocity.

30. Button: 'latex_table' - Prints out a latex script will all the details of all parameters and values currently in the GUI.

31. Button: 'excitation' - Creates an H2 rotational level excitation diagram that can be further used to estimate the gas temperature. Such diagrams are common used in absorption line literature (see e.g. Ranjan+18, 2018A&A...618A.184R)

32. Button: 'save_binary_table' - Saves the table with all parameters and values in a binary format.

33. Checkbox: 'model' - Checkbox to indicate whether the quasar template is visible in the plot.

34. Checkbox: 'data' - Checkbox to indicate whether the spectral data to be analysed is visible in the plot.

35. Checkbox: 'fit' - Checkbox to indicate whether the fit to the spectral data is visible in the plot.

36. Checkbox: 'cont' - Checkbox to indicate whether the drawn/estimated continuum is visible in the plot.

37. Bar: 'DLA-H2-Redshift' - Bar that can be used to vary the H2/HD absorption line redshift

38. Bar: 'DLA-HI-Redshift' - Bar that can be used to vary the HI absorption line redshift

39. Bar: 'QSO-Redshift' - Bar that can be used to vary the QSO emission redshift


The code is free to use for all. Thank you for using the code and please feel free to contact me at - 'ranjan_adarsh@yahoo.com' for any comments, suggestions and bugs.
