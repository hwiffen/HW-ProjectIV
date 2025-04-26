import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

## Define model parameters:
alpha = 1.0
delta = 0.1
F = np.array([1, 1])  

## Define plot parameters:
x = np.linspace(-0.3, 0.3, 1000)
y = np.linspace(-0.3, 0.3, 1000)
X, Y = np.meshgrid(x, y)
r = np.sqrt(X**2 + Y**2)
r[r == 0] = 1e-10 

## Define fundamental functions:
def H1_fundamental(r, alpha):
    return (np.exp(-alpha * r) / (4 * np.pi * r)) * (1 + 1/(alpha*r) + 1/(alpha**2 * r**2)) - 1/(4 * np.pi * alpha**2 * r**3)

def H2_fundamental(r, alpha):
    return -np.exp(-alpha * r) / (4 * np.pi * r**3) * (1 + 3/(alpha*r) + 3/(alpha**2 * r**2)) + 3/(4 * np.pi * alpha**2 * r**5)

## Define regularised functions:
def H1_regularised(r, alpha, delta):
    R_blob = np.sqrt(r**2 + delta**2)
    return (np.exp(-alpha * R_blob) / (4 * np.pi * R_blob)) * (1 + 1/(alpha*R_blob) + 1/(alpha**2 * R_blob**2)) - 1/(4 * np.pi * alpha**2 * R_blob**3) + delta**2 * H2_regularised(r, alpha, delta)

def H2_regularised(r, alpha, delta):
    R_blob = np.sqrt(r**2 + delta**2)
    return -np.exp(-alpha * R_blob) / (4 * np.pi * R_blob**3) * (1 + 3/(alpha*R_blob) + 3/(alpha**2 * R_blob**2)) + 3/(4 * np.pi * alpha**2 * R_blob**5)

## Calculate values for the velocity field:
r_hat_x = X / r
r_hat_y = Y / r
F_dot_r = F[0]*X + F[1]*Y

## Define velocity field function:
def velocity_field(H1, H2):
    u_x = F[0] * H1 + F_dot_r * X * H2
    u_y = F[1] * H1 + F_dot_r * Y * H2
    return u_x, u_y

## Calculate solutions:
H1_f = H1_fundamental(r, alpha)
H2_f = H2_fundamental(r, alpha)
u_fx, u_fy = velocity_field(H1_f, H2_f)

H1_r = H1_regularised(r, alpha, delta)
H2_r = H2_regularised(r, alpha, delta)
u_rx, u_ry = velocity_field(H1_r, H2_r)

## Plot using Brinkman colour scheme: 
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

cmap = cm.get_cmap('winter', 10)
cmap = [cmap(i) for i in range(cmap.N)]

axs[0].streamplot(x, y, u_fx, u_fy, color=cmap[3], linewidth=1)
axs[0].set_title("Fundamental Brinkman flow", fontsize=20)
axs[0].set_aspect('equal')

axs[1].streamplot(x, y, u_rx, u_ry, color=cmap[7], linewidth=1)
axs[1].set_title("Regularised Brinkman flow", fontsize=20)
axs[1].set_aspect('equal')

for ax in axs:
    ax.set_xlim([-0.3, 0.3])
    ax.set_ylim([-0.3, 0.3])
    ax.set_xlabel('x')
    ax.set_ylabel('y')

## Add markers at the origin to represent concentrations of the forces:
axs[0].plot(0, 0, 'ro', markersize=5)

circle = plt.Circle((0, 0), 0.06, edgecolor='red', facecolor='none', linestyle='--', linewidth=1.5)
axs[1].add_patch(circle)

plt.tight_layout()
plt.show()
