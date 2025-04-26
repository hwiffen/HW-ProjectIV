import numpy as np
import time
from scipy.integrate import solve_ivp

## Assign parameters:
N0 = 500  # Number of particles
R = 1  # Radius of initial sphere
alpha = 3  # Permeability parameter
beta = 5 # Magnetic field parameter
t_span = (0, 300)  # Time span for simulation
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
    result[mask] = (-exp_term / (4 * np.pi * r[mask]**3)) * (1 + 3/(alpha*r[mask]) + 3/((alpha**2) * (r[mask]**2))) + 3 / (4 * np.pi * (alpha**2)) * (r[mask]**-5)
    return result


## Define function to compute resultant magnetic force:
def resultant_force(positions, velocities, beta):
    N0 = positions.shape[0] 
    
    r_ij = positions[:, None, :] - positions[None, :, :]
    r_abs = np.linalg.norm(r_ij, axis=2)
    
    mask = r_abs > 1e-2
    r_abs_safe = np.where(mask, r_abs, 1.0)

    B_all = np.zeros_like(r_ij)
    
    valid_i, valid_j = np.where(mask)  
    B_all[valid_i, valid_j] = (1/(N0**2)) * beta * np.cross(velocities[valid_j], r_ij[valid_i, valid_j]) / r_abs_safe[valid_i, valid_j, None]**2

    B_totals = np.sum(B_all, axis=1)

    B_totals[B_totals > 2] = 2 
    B_totals[B_totals < -2] = -2

    F_all = np.cross(velocities, B_totals)

    F_all[:, 2] += -1  # add gravity

    return F_all


## Define function to compute mobility matrix:
def create_mobility_matrix(positions):
 
    r_ij = positions[:, None, :] - positions[None, :, :]
    r_abs = np.linalg.norm(r_ij, axis=2)

    H1_vals = H1(r_abs)
    H2_vals = H2(r_abs)

    r_x, r_y, r_z = r_ij[..., 0], r_ij[..., 1], r_ij[..., 2]

    x_sum = np.stack([H1_vals + r_x**2 * H2_vals, 
                      r_x * r_y * H2_vals, 
                      r_x * r_z * H2_vals], axis=0)
    y_sum = np.stack([r_y * r_x * H2_vals, 
                      H1_vals + r_y**2 * H2_vals, 
                      r_y * r_z * H2_vals], axis=0)
    z_sum = np.stack([r_z * r_x * H2_vals, 
                      r_z * r_y * H2_vals, 
                      H1_vals + r_z**2 * H2_vals], axis=0)

    mobility_matrix = np.stack([x_sum, y_sum, z_sum], axis=0)

    return mobility_matrix

## Vectorised derivative function with progress updates:
def deriv_func(t, r):
    
    global progress_index, start_time, previous_pos, previous_t, velocities

    if progress_index < len(progress_intervals) and t >= progress_intervals[progress_index]:
        elapsed_time = time.time() - start_time
        print(f"Progress: {int((t/t_span[1]) * 100)}% complete. Time elapsed: {elapsed_time:.2f} seconds")
        progress_index += 1

    positions = r.reshape(N0, 3)  

    F_all = resultant_force(positions, velocities, beta)

    mobility_matrix = create_mobility_matrix(positions)

    velocities = np.einsum('uhij,ih->iu', mobility_matrix, F_all) 

    print(f"velocities: {velocities.min()}, {velocities.max()}")

    return velocities.flatten()

## Initialise progress tracking:
progress_index = 0
start_time = time.time()

## Initialise global variables:
init_positions = initial_sphere(R, N0)
previous_pos = init_positions.copy()
previous_t = 0
init_pos_flat = init_positions.flatten()
velocities = np.zeros((N0, 3)) 

## Solve ODE system:
solution = solve_ivp(deriv_func, t_span, init_pos_flat, method="RK45")

## Reshape output for saving:
new_positions = solution.y.reshape(N0, 3, len(solution.t))

## Save to file:
np.save(f"Magnetic_{beta}_{alpha}_{N0}_{t_span[1]}_output.npy", new_positions)
print(f"Magnetic_{beta}_{alpha}_{N0}_{t_span[1]}_output.npy saved" )

