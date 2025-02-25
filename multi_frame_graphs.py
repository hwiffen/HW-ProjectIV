import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as ticker

####################################################################################
"""Parameters to be adjusted to customise plot appearance:"""

## Simulation data to be loaded in for plotting:
Brinkman = True # Set True to plot a Brinkman simulation, False for Stokes
Zoomed_frames = True # Set True to also plot zoomed in versions of selected frames
alpha = 10 
N0 = 2000
T = 500

## Manually select time indices for freeze frames:
timesteps = [20, 56, 114, 226, 438] 

selected_timesteps = [56, 226, 438]  # Frames to plot zoomed-in, if Zoomed_frames = True

## For each desired freezeframe, manually select z-limits to appear:
z_limits = {
    20: (-7, 1),
    56: (-10, 1),
    114: (-20, 1),
    226: (-25, -10),
    438: (-25, -15)
}

## Define global x and y limits:
x_lim = (-6, 6)
y_lim = (-6, 6)

## Define custom axis limits for zoomed-in freeze frames:
zoomed_x_limits = {
    56: (-2, 2),
    226: (-2, 2),
    438: (-4, 4),
}  

zoomed_y_limits = {
    56: (-2, 2),
    226: (-2, 2),
    438: (-4, 4),
}  

zoomed_z_limits = {
    56: (-500, -1),
    226: (-500, -1),
    438: (-500, -1),
} 
####################################################################################

if Brinkman:
    reshaped = np.load(f"npy_output_files/BM_{alpha}_{N0}_{T}_output.npy")
else:
    reshaped = np.load(f"npy_output_files/{N0}_{T}_output.npy")

## Filter for <95th percentile velocity + final position particles for visual clarity:
velocities = np.linalg.norm(reshaped[:, :, -1] - reshaped[:, :, 0], axis=1)
final_positions = reshaped[:, :, -1]
center_of_mass = np.mean(final_positions, axis=0)
distances = np.linalg.norm(final_positions - center_of_mass, axis=1)

velocity_threshold = np.percentile(velocities, 95)
distance_threshold = np.percentile(distances, 95)

mask = (velocities < velocity_threshold) & (distances < distance_threshold)
filtered_positions = reshaped[mask, :, :]
filtered_N = filtered_positions.shape[0]


## Define colour map for plotting: cm.winter for Brinkman, cm.cool for Stokes
if Brinkman:
    time_colors = cm.winter(np.linspace(0, 1, len(timesteps)+1))
else:
    time_colors = cm.cool(np.linspace(0, 1, len(timesteps)+1))    

## Plot selected freeze frames:
fig = plt.figure(figsize=(6, 12))
ax = fig.add_subplot(111, projection='3d')

for i, t in enumerate(timesteps):
    x, y, z = filtered_positions[:, 0, t], filtered_positions[:, 1, t], filtered_positions[:, 2, t]

    z_min, z_max = z_limits.get(t, (-np.inf, np.inf)) 
    valid_mask = (x >= x_lim[0]) & (x <= x_lim[1]) & \
                 (y >= y_lim[0]) & (y <= y_lim[1]) & \
                 (z >= z_min) & (z <= z_max)

    ax.scatter(x[valid_mask], y[valid_mask], z[valid_mask], color=time_colors[-i-1], alpha=1, s=0.1)

## Adjust axis:
ax.set_box_aspect([1, 1, 4])  # 4x stretch in vertical dirction for visual clarity
ax.set_xlabel("X", fontsize=8)
ax.set_ylabel("Y", fontsize=8)
ax.set_zlabel("Z", fontsize=8)

ax.set_xlim(x_lim)
ax.set_ylim(y_lim)

ax.tick_params(labelsize=8)
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))  
ax.yaxis.set_major_locator(ticker.MultipleLocator(5))  

## Plot zoomed-in freeze frames if desired:
if Zoomed_frames:
    for t in selected_timesteps:
        fig_freeze = plt.figure(figsize=(5, 5))
        ax_freeze = fig_freeze.add_subplot(111, projection='3d')

        # Get index of timestep for colour mapping:
        color_idx = timesteps.index(t) 
        frame_color = time_colors[color_idx]

        # Extract positions at this timestep:
        x, y, z = filtered_positions[:, 0, t], filtered_positions[:, 1, t], filtered_positions[:, 2, t]

        # Use custom axis limits for the zoomed-in plot:
        z_min, z_max = zoomed_z_limits.get(t, (-np.inf, np.inf))
        x_min, x_max = zoomed_x_limits.get(t, (-np.inf, np.inf))
        y_min, y_max = zoomed_y_limits.get(t, (-np.inf, np.inf))
        valid_mask = (x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max) & (z >= z_min) & (z <= z_max)

        # Plot zoomed-in graph:
        ax_freeze.scatter(x[valid_mask], y[valid_mask], z[valid_mask], color=frame_color, alpha=1, s=1)

        ax_freeze.set_xlim(x_min, x_max)
        ax_freeze.set_ylim(y_min, y_max)
        ax_freeze.set_zlim(z_min, z_max)

        ax_freeze.set_xticks(np.linspace(x_min, x_max, 5))  
        ax_freeze.set_yticks(np.linspace(y_min, y_max, 5))  
        ax_freeze.set_zticks(np.linspace(z_min, z_max, 5))

        ax_freeze.set_xlabel("X", fontsize=8)
        ax_freeze.set_ylabel("Y", fontsize=8)
        ax_freeze.set_zlabel("Z", fontsize=8)
        ax_freeze.set_title(f"Timestep {t}", fontsize=10)

        plt.show()

plt.show()


