import numpy as np
import pandas as pd

## Set parameters to align with a completed simulation:

N0 = 500  # Number of particles
R = 1  # Initial sphere radius
alpha = 5  # Permeability parameter (Brinkman)
beta = 5  # Strength of magnetic field
t_span = (0,300)  # Simulation duration
Brinkman = False # Set to True for Brinkman flow
Magnetic = True # Set to True for magntised flow

t = t_span[1]

# Load simulation results from .npy files:
if Brinkman:
    reshaped = np.load(f"npy_output_files/BM_{alpha}_{N0}_{t}_output.npy")
elif Magnetic:
    reshaped = np.load(f"npy_output_files/Magnetic_{beta}_{alpha}_{N0}_{t}_output.npy")
else:
    reshaped = np.load(f"npy_output_files/{N0}_{t}_output.npy")

## Filter out outlier particles:
velocities = np.linalg.norm(reshaped[:, :, -1] - reshaped[:, :, 0], axis=1)
final_positions = reshaped[:, :, -1]
center_of_mass = np.mean(final_positions, axis=0)
distances = np.linalg.norm(final_positions - center_of_mass, axis=1)

velocity_threshold = np.percentile(velocities, 95)
distance_threshold = np.percentile(distances, 95)
mask = (velocities < velocity_threshold) & (distances < distance_threshold)

filtered_positions = reshaped[mask, :, :t]
filtered_N = filtered_positions.shape[0]

## Define function to compute widths and save to .csv file:
def compute_widths_vs_height(filtered_positions, n_slices=200, z_min=None, z_max=None):
    N, _, T = filtered_positions.shape
    all_z = filtered_positions[:, 2, :].flatten()

    if z_min is None:
        z_min = np.percentile(all_z, 1)
    if z_max is None:
        z_max = np.percentile(all_z, 99)

    z_levels = np.linspace(z_min, z_max, n_slices)
    results = []

    for z in z_levels:
        x_vals_at_z = []
        y_vals_at_z = []

        for i in range(N):
            z_traj = filtered_positions[i, 2, :]
            below = np.where(z_traj < z)[0]
            if below.size > 0:
                idx = below[0]
                x = filtered_positions[i, 0, idx]
                y = filtered_positions[i, 1, idx]
                x_vals_at_z.append(x)
                y_vals_at_z.append(y)

        if len(x_vals_at_z) >= 5:
            x_vals_at_z = np.array(x_vals_at_z)
            y_vals_at_z = np.array(y_vals_at_z)

            x_trimmed = x_vals_at_z[
                (x_vals_at_z >= np.percentile(x_vals_at_z, 5)) &
                (x_vals_at_z <= np.percentile(x_vals_at_z, 95))
            ]
            y_trimmed = y_vals_at_z[
                (y_vals_at_z >= np.percentile(y_vals_at_z, 5)) &
                (y_vals_at_z <= np.percentile(y_vals_at_z, 95))
            ]

            if len(x_trimmed) >= 5 and len(y_trimmed) >= 5:
                x_width = np.max(x_trimmed) - np.min(x_trimmed)
                y_width = np.max(y_trimmed) - np.min(y_trimmed)
                results.append((z, x_width, y_width, len(x_trimmed)))

    if Brinkman:
        filename = f"BM_{alpha}_{N0}_{t}_particle_widths.csv"
    elif Magnetic:
        filename = f"Magnetic_{beta}_{alpha}_{N0}_{t}_particle_widths.csv"
    else:
        filename = f"Stokes_{N0}_{t}_particle_widths.csv"

    df = pd.DataFrame(results, columns=["z", "x_width", "y_width", "n_particles"])
    df.to_csv(filename, index=False)
    print(f"Saved data to {filename}")

## Run the function on specified simualtion data:
compute_widths_vs_height(filtered_positions, n_slices=300)
