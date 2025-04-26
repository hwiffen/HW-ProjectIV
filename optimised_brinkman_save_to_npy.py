import numpy as np
import time
from scipy.integrate import solve_ivp


## Assign parameters:
N0 = 500  # Number of particles
R = 1  # Radius of initial sphere
alpha = 5  # Permeability parameter
t_span = (0, 500)  # Time span for simulation
progress_intervals = np.linspace(t_span[0], t_span[1], 101)[1:] # Intervals for progress printing

## Define initial sphere function:
def initial_sphere(R, N):
    points = []
    while len(points) < N:
        point = np.random.uniform(-R, R, 3)
        if np.linalg.norm(point) <= R:
            points.append(point)
    return np.array(points)

## Define vectorised H1 and H2 functions:
def H1(r):
    mask = r > 1e-6  # Avoid r -> infinity error
    result = np.zeros_like(r)
    exp_term = np.exp(-alpha * r[mask])
    result[mask] = (exp_term / (4 * np.pi * r[mask])) * (1 + 1/(alpha*r[mask]) + 1/((alpha**2) * (r[mask]**2))) - 1/(4 * np.pi * (alpha**2) * (r[mask]**3))
    return result

def H2(r):
    mask = r > 1e-6 # Avoid r -> infinity error
    result = np.zeros_like(r)
    exp_term = np.exp(-alpha * r[mask])
    result[mask] = (-exp_term / (4 * np.pi * r[mask]**3)) * (1 + 3/(alpha*r[mask]) + 3/((alpha**2) * (r[mask]**2))) + 3 / (4 * np.pi * (alpha**2) * (r[mask]**5))
    return result

## Vectorised derivative function with progress updates:
def deriv_func(t, r):
    
    global progress_index, start_time

    if progress_index < len(progress_intervals) and t >= progress_intervals[progress_index]:
        elapsed_time = time.time() - start_time
        print(f"Progress: {int((t/t_span[1]) * 100)}% complete. Time elapsed: {elapsed_time:.2f} seconds")
        progress_index += 1

    positions = r.reshape(N0, 3)  
    velocities = np.zeros_like(positions)  

    r_ij = positions[:, None, :] - positions[None, :, :] 
    r_abs = np.linalg.norm(r_ij, axis=2) 
 
    H1_vals = H1(r_abs)
    H2_vals = H2(r_abs)

    vel_x = -r_ij[:, :, 0] * r_ij[:, :, 2] * H2_vals
    vel_y = -r_ij[:, :, 1] * r_ij[:, :, 2] * H2_vals
    vel_z = -H1_vals - (r_ij[:, :, 2] ** 2) * H2_vals

    velocities[:, 0] = np.sum(vel_x * (r_abs > 1e-6), axis=1)
    velocities[:, 1] = np.sum(vel_y * (r_abs > 1e-6), axis=1)
    velocities[:, 2] = np.sum(vel_z * (r_abs > 1e-6), axis=1)

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

## Save to file:
np.save(f"BM_{alpha}_{N0}_{t_span[1]}_output.npy", reshaped)
print(f"3D array saved to BM_{alpha}_{N0}_{t_span[1]}_output.npy")

