import numpy as np
import plotly.graph_objects as go
from matplotlib import cm

####################################################################################
"""Parameters to determine which simulation to plot:"""

N0 = 500  # Number of particles
R = 1  # Initial sphere radius
alpha = 5  # Permeability parameter (Brinkman)
beta = 5 # Strength of magnetic field
t_span = (0,500)  # Simulation duration
Brinkman = True # Set to True for Brinkman flow
Magnetic = False # Set to True for magntised flow

t = t_span[1]
max_frames = 500 # Set maximum number of frames to include in animation
####################################################################################

## Load simulation results:
if Brinkman:
    reshaped = np.load(f"npy_output_files/BM_{alpha}_{N0}_{t}_output.npy")
elif Magnetic:
    reshaped = np.load(f"npy_output_files/Magnetic_{beta}_{alpha}_{N0}_{t}_output.npy")
else:
    reshaped = np.load(f"npy_output_files/{N0}_{t}_output.npy")

N, _, T = reshaped.shape 

T = min(T, max_frames) 

## Filter out outlier particles:
velocities = np.linalg.norm(reshaped[:, :, -1] - reshaped[:, :, 0], axis=1)

final_positions = reshaped[:, :, -1]
center_of_mass = np.mean(final_positions, axis=0)

distances = np.linalg.norm(final_positions - center_of_mass, axis=1)

velocity_threshold = np.percentile(velocities, 78)
distance_threshold = np.percentile(distances, 85)

mask = (velocities < velocity_threshold) & (distances < distance_threshold)

filtered_positions = reshaped[mask, :, :T]
filtered_N = filtered_positions.shape[0]

## Generate random colour maps using consistent colour maps for each simulation type:
if Brinkman:
    rand_colors = cm.winter(np.random.rand(filtered_N))  
    hex_colors = ["rgba({},{},{},{})".format(
        int(c[0] * 255), int(c[1] * 255), int(c[2] * 255), c[3]) for c in rand_colors]
    
elif Magnetic:
    rand_colors = cm.autumn(np.random.rand(filtered_N)) 
    hex_colors = ["rgba({},{},{},{})".format(
        int(c[0] * 255), int(c[1] * 255), int(c[2] * 255), c[3]) for c in rand_colors]

else:
    rand_colors = cm.cool(np.random.rand(filtered_N)) 
    hex_colors = ["rgba({},{},{},{})".format(
        int(c[0] * 255), int(c[1] * 255), int(c[2] * 255), c[3]) for c in rand_colors]

## Create list of particle positions for each frame:
frames = []
for frame in range(T):
    frame_data = go.Scatter3d(
        x=filtered_positions[:, 0, frame],
        y=filtered_positions[:, 1, frame],
        z=filtered_positions[:, 2, frame],
        mode="markers",
        marker=dict(size=3, opacity=0.7, color=hex_colors),
    )
    frames.append(go.Frame(data=[frame_data], name=str(frame)))

## Define custom title based on simulation type and parameters:
if Brinkman:
    custom_title=str(f"Brinkman flow simulation for α={alpha}, N0={N0}")

elif Magnetic:
    custom_title=str(f"Brinkman flow in magnetic field simulation for β={beta}, α={alpha}, N0={N0}")

else:
    custom_title=str(f"Stokes flow simulation for N0={N0}")

## Create interactive 3D plot:
fig = go.Figure(
    data=[go.Scatter3d(
        x=filtered_positions[:, 0, 0],
        y=filtered_positions[:, 1, 0],
        z=filtered_positions[:, 2, 0],
        mode="markers",
        marker=dict(size=3, opacity=0.7, color=hex_colors),
    )],
    layout=go.Layout(
        title=custom_title,      
        updatemenus=[{
            "buttons": [
                {
                    "args": [None, {"frame": {"duration": 50, "redraw": True}, "fromcurrent": True}],
                    "label": "▶ Play",
                    "method": "animate",
                },
                {
                    "args": [[None], {"frame": {"duration": 0, "redraw": False}, "mode": "immediate", "transition": {"duration": 0}}],
                    "label": "❚❚ Pause",
                    "method": "animate",
                },
            ],
            "direction": "left",
            "pad": {"r": 10, "t": 87},
            "showactive": False,
            "type": "buttons",
            "x": 0.1,
            "xanchor": "right",
            "y": 0,
            "yanchor": "top",
        }],
        sliders=[{
            "active": 0,
            "steps": [
                {"args": [[str(k)], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}], "label": str(k), "method": "animate"}
                for k in range(T)
            ],
        }],
    ),
    frames=frames,
)

## Show interactive plot:
fig.show()

## Save interactive plot as a HTML file:
if Brinkman:
    html_filename = f"BM_{alpha}_{N0}_{t}_particle_motion.html"
elif Magnetic:
    html_filename = f"Magnetic_{beta}_{alpha}_{N0}_{t}_particle_motion.html"
else:
    html_filename = f"Stokes_{N0}_{t}_particle_motion.html"

fig.write_html(html_filename)
print(f"Interactive plot saved as {html_filename}")