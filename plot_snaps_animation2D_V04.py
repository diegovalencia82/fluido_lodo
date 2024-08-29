import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

def read_snapshot(filepath):
    # Leer datos del archivo snapshot
    with open(filepath, 'r') as file:
        lines = file.readlines()
        # Extraer el tiempo del segundo dato de la primera línea que contiene metadatos
        time = float(lines[0].split()[1])
        # Saltar la primera línea que contiene metadatos
        lines = lines[1:]
        # Extraer las posiciones x e y, el tipo y las velocidades de cada partícula de cada línea
        data = np.array([[float(val) for val in line.split()[1:]] for line in lines])
    return time, data

def calculate_velocity(data):
    # Calcular la magnitud de la velocidad para cada partícula
    velocity = np.sqrt(data[:, 2]**2 + data[:, 3]**2)  # Considerando las columnas 4 y 5 como las velocidades
    return velocity

def update(frame):
    ax.clear()
    ax.set_xlim(-0.1, 1.7)  # Ajustar límites x según tus necesidades
    ax.set_ylim(-0.1, 0.7)  # Ajustar límites y según tus necesidades
    
    # Calcular la velocidad de las partículas en el frame actual
    velocity = calculate_velocity(frame[1])
    
    # Crear una escala de colores que vaya de azul a rojo
    norm = Normalize(vmin=0, vmax=np.sqrt(2*9.8*0.6))
    cmap = plt.get_cmap('coolwarm')
    sm = ScalarMappable(norm=norm, cmap=cmap)
    
    # Obtener el color de cada partícula según su velocidad
    colors = sm.to_rgba(velocity)
    
    # Dibujar las partículas en el frame actual con el color correspondiente
    ax.scatter(frame[1][:, 0], frame[1][:, 1], color=colors, s=20)

    # Seleccionar las partículas con itype = -1 y pintarlas de color gris
    gray_particles = frame[1][frame[1][:, 8] == -1]
    ax.scatter(gray_particles[:, 0], gray_particles[:, 1], color='gray', s=30)
    
    # Mostrar el tiempo en la gráfica
    time_factor = 1.#2.6 / 3.75
    adjusted_time = frame[0] * time_factor
    ax.text(0.05, 0.95, f't = {adjusted_time:.2f}', transform=ax.transAxes, ha='left', va='top', fontsize=12)

# Definir la carpeta donde se encuentran los archivos
folder = './'

# Crear una figura y un eje para la animación
fig, ax = plt.subplots()

# Lista para almacenar los datos de todas las instantáneas
data_frames = []

# Iterar sobre los nombres de archivo del tipo 'snapshot_0001', 'snapshot_0002', ..., 'snapshot_2000'
for i in range(1, 800):
    # Formatear el nombre de archivo con el número de snapshot
    filename = f'snapshot_{i:04d}'
    # Crear la ruta completa al archivo
    filepath = os.path.join(folder, filename)
    # Verificar si el archivo existe
    if os.path.exists(filepath):
        # Leer datos del archivo snapshot y almacenarlos en la lista
        time, data = read_snapshot(filepath)
        data_frames.append((time, data))
    else:
        print(f'El archivo {filename} no existe en la carpeta especificada.')

# Crear la animación
anim = FuncAnimation(fig, update, frames=data_frames, interval=10)  # Intervalo de 5 milisegundos entre cada frame

# Añadir la barra de colores en la parte superior izquierda
cbar_ax = fig.add_axes([0.15, 0.92, 0.7, 0.02])  # [left, bottom, width, height]
norm = Normalize(vmin=0, vmax=np.sqrt(2*9.8*0.6))
cmap = plt.get_cmap('coolwarm')
sm = ScalarMappable(norm=norm, cmap=cmap)
cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')
cbar_ax.text(1.0, 1.3, 'Velocity magnitude', transform=cbar_ax.transAxes, ha='right', va='bottom', fontsize=12)

# Guardar la animación como un archivo .mp4
anim.save('animacion.mp4', writer='ffmpeg')

plt.tight_layout(rect=[0, 0, 1, 0.9])  # Ajustar el layout para dejar espacio a la barra de colores
plt.show()

