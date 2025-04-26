import numpy as np
import time
from scipy.integrate import solve_ivp

## Assign parameters:
alpha = 1  # Permeability parameter
delta = 0.1 # Regularisation parameter
K = 1  # Spring constant
L = 1  # Rest length of springs
t_span = (0, 500)  # Time span for simulation
progress_intervals = np.linspace(t_span[0], t_span[1], 101)[1:] # Intervals for progress printing

## Establish initial stickman position:
"""
Stickman made up of 32 points, in the x-z plane.
Connected points are ~ 1 unit apart (with small, random perturbations to allow 3D movement).
"""

# Head:
head_x = np.array([-0.9238795325,  -0.9238795325 - 0.3826834324,  -0.9238795325,  0, 0.9238795325, 0.9238795325 + 0.3826834324, 
                   0.9238795325,  0])
head_y = np.array([-1+0.3826834324,  -1+0.3826834324+0.9238795325,  -1+0.3826834324+2*0.9238795325,  -1+2*0.3826834324+2*0.9238795325,
                   -1+0.3826834324+2*0.9238795325,  -1+0.3826834324+0.9238795325,  -1+0.3826834324,  -1])

# Body:
body_x = np.array([0, 0, 0, 0, 0, 0])
body_y = np.array([-2, -3, -4, -5, -6, -7])

# Arms:
arms_x = np.array([-4, -3, -2, -1, 1, 2, 3, 4])
arms_y = np.array([-3, -3, -3, -3, -3, -3, -3, -3])

# Legs:
legs_x = np.array([-3/5, -6/5, -9/5, -12/5, -3, 3/5, 6/5, 9/5, 12/5, 3])
legs_y = np.array([-39/5, -43/5, -47/5, -51/5, -11, -39/5, -43/5, -47/5, -51/5, -11])

stickman_positions = np.zeros((32,3))
stickman_positions[:,0] = np.concatenate([head_x, body_x, arms_x, legs_x])
stickman_positions[:,2] = np.concatenate([head_y, body_y, arms_y, legs_y])

# Add small random pertubations to y-coordinates:
stickman_positions[:,1] = np.random.uniform(-0.1, 0.1, 32)

## Define spring connections between stickman points:
connections = [(i, i+1) for i in range(31)]

connections.remove((13, 14))
connections.remove((17, 18))
connections.remove((21, 22))
connections.remove((26, 27))

connections.append((0,7))
connections.append((9,17))
connections.append((9,18))
connections.append((13,22))
connections.append((13,27))

## Define regularised Brinkmanlet interaction terms:
def H2(r):
    mask = r > 1e-6 # Avoid r -> infinity error
    R = np.sqrt(r[mask]**2 + delta**2)
    result = np.zeros_like(r)
    exp_term = np.exp(-alpha * R)
    result[mask] = (-exp_term / (4 * np.pi * R**3)) * (1 + 3/(alpha*R) + 3/((alpha**2) * (R**2))) + 3 / (4 * np.pi * (alpha**2)) * (R**-5)
    return result

def H1(r):
    mask = r > 1e-6  # Avoid r -> infinity error
    R = np.sqrt(r[mask]**2 + delta**2)
    result = np.zeros_like(r)
    exp_term = np.exp(-alpha * R)
    H2_term = H2(r[mask])
    result[mask] = (exp_term / (4 * np.pi * R)) * (1 + 1/(alpha*R) + 1/((alpha**2) * (R**2))) - 1/(4 * np.pi * (alpha**2) * (R**3)) + (delta**2 * H2_term)
    return result


## Define function to compute resultant spring force:
def resultant_force(positions, connections, L, K):
    
    F_total = np.zeros_like(positions)

    for i, j in connections:
        r_ij = positions[i] - positions[j]
        r_abs = np.linalg.norm(r_ij)
        if r_abs > 1e-6: 
            r_hat = r_ij / r_abs
            F_ij = -K * (r_abs - L) * r_hat  
            F_total[i] += F_ij
            F_total[j] -= F_ij

    F_total[:, 2] += -1 # add gravity

    return F_total

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
    
    global progress_index, start_time

    if progress_index < len(progress_intervals) and t >= progress_intervals[progress_index]:
        elapsed_time = time.time() - start_time
        print(f"Progress: {int((t/t_span[1]) * 100)}% complete. Time elapsed: {elapsed_time:.2f} seconds")
        progress_index += 1

    positions = r.reshape(N0, 3)  

    F_all = resultant_force(positions, connections, L, K)

    mobility_matrix = create_mobility_matrix(positions)

    velocities = np.einsum('uhij,ih->iu', mobility_matrix, F_all) 

    return velocities.flatten()

## Initialise global variables:
progress_index = 0
start_time = time.time()
N0 = stickman_positions.shape[0]
init_pos_flat = stickman_positions.flatten()

## Solve ODE:
solution = solve_ivp(deriv_func, t_span, init_pos_flat, method="RK45")

## Reshape output for saving:
new_positions = solution.y.reshape(N0, 3, len(solution.t))

## Save to file:
np.save(f"Stickman_{alpha}_{delta}_{K}_{L}_{t_span[1]}_output.npy", new_positions)
print(f"Stickman_{alpha}_{delta}_{K}_{L}_{t_span[1]}_output.npy saved" )
