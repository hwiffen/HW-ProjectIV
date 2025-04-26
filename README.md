# HW-ProjectIV

**Code repository by Harry Wiffen for Project IV: *Stokes and the Brinkmanlets;* for the degree of Master of Mathematics at Durham University.**

This repository contains all code used in ***Stokes and the Brinkmanlets***. All code was created using Python or R. A breakdown of the contents can be seen below.

Throughout, consistent colour maps are used to distinguish different set of simulations. From the matplotlib package, the following colour maps are used:

-   Stokes flows: 'cool'

-   Brinkman flows: 'winter'

-   Magnetised flows: 'autumn'

## Folders:

### npy_output_files

Simulation results from all models. Can be imported into visualisation scripts.

**Note:** Many longer simulations are excluded, due to file size limitations.

### Simulation_scripts 

**optimised_stokes_save_to_npy.py** - After selecting simulation parameters, this script simulates a falling cloud in a Stokes flow, saving the results as a .npy file.

**optimised_brinkman_save_to_npy.py** - After selecting simulation parameters, this script simulates a falling cloud in a Brinkman flow, saving the results as a .npy file.

**optimised_magnetised_save_to_npy.py** - After selecting simulation parameters, this script simulates a falling cloud in a magnetised Brinkman flow, saving the results as a .npy file.

**articulated_body_save_to_npy.py** - After selecting simulation parameters, this script simulates a falling 3D articulated body in a Brinkman flow, saving the results as a .npy file.

### Main_visualisation_scripts

In both of these scripts, simulation parameters must be specified. These inputted values must align with a completed simulation. The .npy file from this simulation is then read in.

**interactive_online_plots.py** - An interactive, online 3D plot is produced using plotly. A HTML file is also produced, containing the interactive plot.

**multi_frame_graphs.py** - Specific time-steps must be specified, along with the desired z_limits at each time. The simulation results from each time-step are then superposed into one holistic visual. There is also an option to specify a number of zoomed-in frames. Plots of these specific time-steps are then returned in addition to the holistic plot.

### Spread_analysis

These scripts all relate to the calculation and creation of the particle dispersion plots. The raw data is calculated in Python, before being saved as a .csv file and subsequently imported into R, for plotting.

**spread_analysis_save_to_csv.py** - After specifying simulation parameters, the corresponding .npy file is imported. A .csv file is created, with the average particle width in each vertical window.

**Raw_spread_data_vs_LLRE.R** - Imports a specified .csv file. Defines a function to evaluate the local linear regression estimator (LLRE). Plots the raw data alongside the LLRE of the data, for visual comparison.

**Stokes_spread_analysis.R** - Imports the specified .csv files. Defines a function to evaluate the local linear regression estimator (LLRE). Plots the LLRE of each .csv, with corresponding labels, for Stokes flow simulations.

**Brinkman_spread_analysis.R** - Imports the specified .csv files. Defines a function to evaluate the local linear regression estimator (LLRE). Plots the LLRE of each .csv, with corresponding labels, for Brinkman flow simulations.

**Magnetised_spread_analysis.R** - Imports the specified .csv files. Defines a function to evaluate the local linear regression estimator (LLRE). Plots the LLRE of each .csv, with corresponding labels, for magnetised flow simulations.

### Misc_visualisation_scripts

**scalar_blobs_visualisation.R** - Plots multiple scalar blob functions, as Gaussian curves, with variance $\delta^2$. Included a visualisation of the Dirac delta measure for comparison.

**regularised_soln_blob_heatmaps.R** - Plots heatmaps of the precise 3D blob function found for the regularised Brinkman equations, for different values of $\delta$.

**streamlines_fundamental_vs_regularised_brinkman.py** - Plots streamlines of fundamental and regularised Brinkman flow, with circles representing the concentrations of the forces.

