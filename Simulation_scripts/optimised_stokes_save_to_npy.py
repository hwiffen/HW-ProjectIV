import numpy as np
import time
from scipy.integrate import solve_ivp

## Assign parameters:
N0 = 300  # Number of particles
R = 1  # Radius of initial sphere
t_span = (0, 1000)  # Simulation time span
progress_intervals = np.linspace(t_span[0], t_span[1], 101)[1:] # Intervals for progress printing

## Define initial sphere function:
def initial_sphere(R, N):
    points = []
    while len(points) < N:
        point = np.random.uniform(-R, R, 3)
        if np.linalg.norm(point) <= R:
            points.append(point)
    return np.array(points)

## Vectorised derivative function with progress updates:
def deriv_func(t, r):

    global progress_index, start_time

    if progress_index < len(progress_intervals) and t >= progress_intervals[progress_index]:
        elapsed_time = time.time() - start_time
        print(f"Progress: {int((t/t_span[1]) * 100)}% complete. Time elapsed: {elapsed_time:.2f} seconds")
        progress_index += 1
    r = r.reshape(N0, 3) 
    velocities = np.zeros_like(r)

    r_i = r[:, np.newaxis, :]  
    r_j = r[np.newaxis, :, :]
    
    r_vec = r_i - r_j  
    r_abs = np.linalg.norm(r_vec, axis=2)  
    
    mask = (r_abs > 1e-6)  # Avoid r -> infinity error
    r_abs[~mask] = 1  

    r_abs3 = r_abs ** 3
    sum_i = np.zeros_like(r)

    sum_i[:, 0] = np.sum(mask * (-r_vec[:, :, 0] * r_vec[:, :, 2] / r_abs3), axis=1)
    sum_i[:, 1] = np.sum(mask * (-r_vec[:, :, 1] * r_vec[:, :, 2] / r_abs3), axis=1)
    sum_i[:, 2] = np.sum(mask * ((-1 / r_abs) - (r_vec[:, :, 2] ** 2) / r_abs3), axis=1)

    velocities = (5 / (8 * N0)) * sum_i  
    return velocities.flatten()

## Initialise positions:
init_positions = initial_sphere(R, N0)

## Flatten initial positions for ODE solver:
init_pos_flat = init_positions.flatten()

## Initialise progress tracking:
progress_index = 0
start_time = time.time()

## Solve ODE system:
solution = solve_ivp(deriv_func, t_span, init_pos_flat, method="RK45")

## Reshape output for saving:
reshaped = solution.y.reshape(N0, 3, len(solution.t))

## Save to .npy file:
output_filename = f"{N0}_{t_span[1]}_output.npy"
np.save(output_filename, reshaped)
print(f"3D array saved to {output_filename}")
