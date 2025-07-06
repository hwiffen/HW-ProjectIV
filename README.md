![GitHub repo size](https://img.shields.io/github/repo-size/hwiffen/HW-ProjectIV)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Made with Python](https://img.shields.io/badge/Made%20with-Python-blue)
![Made with R](https://img.shields.io/badge/Made%20with-R-1f425f)

# HW-ProjectIV

**Repository by Harry Wiffen for Project IV: *Stokes and the Brinkmanlets*, submitted for the degree of Master of Mathematics at Durham University.**

This repository contains all code used in ***Stokes and the Brinkmanlets***. Code is written in **`Python`** and **`R`**.

The complete dissertation report can be found at [`Stokes_and_the_Brinkmanlets.pdf`](Stokes_and_the_Brinkmanlets.pdf)

---

Consistent colour maps are used throughout to distinguish simulation types (from `matplotlib`):

- **Stokes flows**: `'cool'`
- **Brinkman flows**: `'winter'`
- **Magnetised flows**: `'autumn'`

---

## Folder Structure


### [`HTML_interactive_plots/`](HTML_interactive_plots)

3D interactive simulation results, saved as `.HTML` files. Files were written using `write_html` from the `plotly` package.

> **Note:** Only one plot of each simulation type is included, due to GitHub file size limits. Any other plots can be produced directly using **[`interactive_online_plots.py`](Main_visualisation_scripts/interactive_online_plots.py)**.

---

### [`npy_output_files/`](npy_output_files)

Raw simulation results from all models, saved as `.npy` files. Can be imported into visualisation scripts.

> **Note:** Many longer simulations are excluded, due to GitHub file size limits.

---

### [`Simulation_scripts/`](Simulation_scripts)

Scripts for running simulations and saving outputs as `.npy` files. Simulation parameters must be specified:

- **[`optimised_stokes_save_to_npy.py`](Simulation_scripts/optimised_stokes_save_to_npy.py)** — Simulates a falling particle cloud in Stokes flow.
- **[`optimised_brinkman_save_to_npy.py`](Simulation_scripts/optimised_brinkman_save_to_npy.py)** — Simulates a falling cloud in Brinkman flow.
- **[`optimised_magnetised_save_to_npy.py`](Simulation_scripts/optimised_magnetised_save_to_npy.py)** — Simulates a falling cloud in a magnetised Brinkman flow.
- **[`articulated_body_save_to_npy.py`](Simulation_scripts/articulated_body_save_to_npy.py)** — Simulates a falling 3D articulated body in Brinkman flow.

---

### [`Main_visualisation_scripts/`](Main_visualisation_scripts)

Scripts for visualising simulation results. Simulation parameters must be specified, and the corresponding `.npy` file is then imported:

- **[`interactive_online_plots.py`](Main_visualisation_scripts/interactive_online_plots.py)** — Produces an online, interactive 3D plot with `plotly` and saves it as an `.HTML` file.
- **[`multi_frame_graphs.py`](Main_visualisation_scripts/multi_frame_graphs.py)** — Creates holistic multi-frame visualisations. Option to also return separate zoomed-in frames.

---

### [`Spread_analysis/`](Spread_analysis)

Scripts for analysing and visualising particle dispersion. Raw data is calculated using **`Python`** and saved as a `.csv` file, to be imported into **`R`** for visualisation:

- **[`spread_analysis_save_to_csv.py`](Spread_analysis/spread_analysis_save_to_csv.py)** — Processes simulation data and saves particle spread data as `.csv`.
- **[`Raw_spread_data_vs_LLRE.R`](Spread_analysis/Raw_spread_data_vs_LLRE.R)** — Plots raw dispersion data with the local linear regression estimator (LLRE), for visual comparison.
- **[`Stokes_spread_analysis.R`](Spread_analysis/Stokes_spread_analysis.R)** — LLRE spread plots for Stokes flow simulations.
- **[`Brinkman_spread_analysis.R`](Spread_analysis/Brinkman_spread_analysis.R)** — LLRE spread plots for Brinkman flow simulations.
- **[`Magnetised_spread_analysis.R`](Spread_analysis/Magnetised_spread_analysis.R)** — LLRE spread plots for magnetised flow simulations.

---

### [`Misc_visualisation_scripts/`](Misc_visualisation_scripts)

Additional visualisation scripts:

- **[`scalar_blobs_visualisation.R`](Misc_visualisation_scripts/scalar_blobs_visualisation.R)** — Plots multiple Gaussian scalar blob functions, with the Dirac delta measure for comparison.
- **[`regularised_soln_blob_heatmaps.R`](Misc_visualisation_scripts/regularised_soln_blob_heatmaps.R)** — Heatmaps of regularised Brinkman blob solutions for different width parameters.
- **[`streamlines_fundamental_vs_regularised_brinkman.py`](Misc_visualisation_scripts/streamlines_fundamental_vs_regularised_brinkman.py)** — Compares streamlines of fundamental and regularised Brinkman solutions, with a visual representation of the force concentrations.

---

## License

This project is licensed under the [MIT License](LICENSE).

---
