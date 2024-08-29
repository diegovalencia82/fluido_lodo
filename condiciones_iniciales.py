import numpy as np
import matplotlib.pyplot as plt

# Parámetros para el espacio de partículas
particle_spacing = 0.01  # Espaciado entre partículas del fluido
delta = particle_spacing  # Pequeño delta para no sobreponerse con las paredes

# Parámetros de la caja y del bloque de agua
box_length = 1.61
box_height = 0.6

# Parámetros para el bloque de agua
water_length = 0.6
water_height = 0.3

# Constantes físicas
g = 9.82
gamma = 7.0
beta0 = 25.0  # Ajuste según sea necesario
beta = beta0 * 1.0
ht = water_height
c = beta * np.sqrt(2. * g * ht)
rho0 = 1000.
b = rho0 * c * c / gamma

# Crear partículas para el bloque de agua
x_water = np.arange(1.5 * delta, water_length, particle_spacing)
z_water = np.arange(1.5 * delta, water_height, particle_spacing)
xx_water, zz_water = np.meshgrid(x_water, z_water)
water_particles = np.vstack([xx_water.ravel(), zz_water.ravel()]).T

# Crear partículas para las paredes de la caja
x_box = np.arange(0, box_length + particle_spacing, particle_spacing)
z_box = np.arange(particle_spacing, box_height + particle_spacing, particle_spacing)

wall_particles = []
# Paredes laterales
for x in [0, box_length]:
    for z in z_box:
        wall_particles.append([x, z])
# Piso
for x in x_box:
    wall_particles.append([x, 0])

wall_particles = np.array(wall_particles)

# Combinar partículas de la caja y del agua
num_fluid_particles = water_particles.shape[0]
num_box_particles = wall_particles.shape[0]
num_particles = num_fluid_particles + num_box_particles

# Inicializar IDs y velocidades
ids_fluid = np.arange(1, num_fluid_particles + 1)
ids_wall = np.arange(num_fluid_particles + 1, num_particles + 1)
velocities = np.zeros((num_particles, 2))

# Densidad, presión y masa de las partículas
dx = particle_spacing
rho = np.zeros(num_particles)
pressure = np.zeros(num_particles)
mass = np.zeros(num_particles)
for i in range(num_fluid_particles):
    rho[i] = rho0 * (1 + (rho0 * g * (ht - abs(water_particles[i, 1])) / b)) ** (1.0 / gamma)
    pressure[i] = b * ((rho[i] / rho0) ** gamma - 1.0)
    mass[i] = dx * dx * rho[i]

# Para las partículas de pared, mantener rho0
rho[num_fluid_particles:] = rho0
pressure[num_fluid_particles:] = 0.0
mass[num_fluid_particles:] = dx * dx * rho0

# Presión y energía interna
internal_energy = np.full(num_particles, 357.1)

# Tipo de partícula y hsml
itype = np.full(num_particles, -1)
itype[:num_fluid_particles] = 1  # Partículas del fluido
hsml = np.full(num_particles, dx)

# Guardar los datos en un archivo
with open('snapshot_000', 'w') as f:
    f.write(f"0   0.0   {num_particles}   {num_fluid_particles}   {num_box_particles}\n")

    # Escribir partículas del fluido
    for i in range(num_fluid_particles):
        f.write(f"{ids_fluid[i]}   {water_particles[i, 0]:.6f}   {water_particles[i, 1]:.6f}   "
                f"{velocities[i, 0]:.6f}   {velocities[i, 1]:.6f}   {mass[i]:.6f}   {rho[i]:.6f}   "
                f"{pressure[i]:.6f}   {internal_energy[i]:.6f}   {itype[i]}   {hsml[i]:.6f}\n")

    # Escribir partículas de las paredes
    for i in range(num_box_particles):
        j = num_fluid_particles + i
        f.write(f"{ids_wall[i]}   {wall_particles[i, 0]:.6f}   {wall_particles[i, 1]:.6f}   "
                f"{velocities[j, 0]:.6f}   {velocities[j, 1]:.6f}   {mass[j]:.6f}   {rho[j]:.6f}   "
                f"{pressure[j]:.6f}   {internal_energy[j]:.6f}   {itype[j]}   {hsml[j]:.6f}\n")

# Graficar las partículas en 2D
fig, ax = plt.subplots()

# Paredes de la caja (puntos más pequeños y transparentes)
ax.scatter(wall_particles[:, 0], wall_particles[:, 1], c='r', s=10, alpha=0.1, label='Box')

# Bloque de agua (puntos más grandes)
ax.scatter(water_particles[:, 0], water_particles[:, 1], c='b', s=10, label='Water')

ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.legend()
ax.set_aspect('equal')  # Asegurar que los ejes tengan la misma escala
plt.show()

