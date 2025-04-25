import numpy as np
import matplotlib.pyplot as plt

# Parámetros de la viga
q = 20  # Carga distribuida en kN/m
L = 5   # Longitud de la viga en metros
w = q   # Carga por unidad de longitud (kN/m)

# Fuerzas cortantes en los nodos
V1 = 10  # Reacción en el primer apoyo (kN)
V2 = 10  # Reacción en el segundo apoyo (kN)

# Momentos en los nodos
M1 = 0  # Momento en el primer apoyo (kN.m)
M2 = 0 # Momento en el segundo apoyo (kN.m)

# Suponemos que la fuerza axial es constante a lo largo de la viga
A1 = 5  # Fuerza axial en el primer nodo (kN)
A2 = 5  # Fuerza axial en el segundo nodo (kN)

# Crear una malla de puntos a lo largo de la viga
x = np.linspace(0, L, 100)

# Interpolación de la fuerza cortante
V = V1 - (x * (V1 - V2) / L)  # La fuerza cortante varía linealmente a lo largo de la viga
V_scaled = V * 0.1  # Escala para el diagrama

# Calcular el momento integrando la fuerza cortante
M = -np.cumsum(V) * (x[1] - x[0])  # Aproximación de la integral
M -= M[0]  # Ajustar para que M(x=0) = 0

# Sumamos M1 y M2 al momento en los extremos
M += M1  # Añadir el momento en el primer nodo
M_scaled = M * 0.1  # Aplicar un factor de escala para el diagrama de momento

# Calcular la fuerza axial (suponiendo que es constante a lo largo de la viga)
A = A1  # Asumiendo fuerza axial constante
A_scaled = np.full_like(x, A * 0.1)  # Escalar la fuerza axial para la visualización

# Graficar los tres diagramas
fig, ax = plt.subplots(3, 1, figsize=(8, 9))  # Ajustar el tamaño de la figura

# Diagrama de esfuerzo cortante
ax[0].plot(x, V_scaled, label="Fuerza Cortante", color="red")
ax[0].set_title("Diagrama de Fuerza Cortante")
ax[0].set_xlabel("Posición a lo largo de la viga (m)")
ax[0].set_ylabel("Esfuerzo Cortante (kN)")
ax[0].grid(True)

# Diagrama de momento con escala aplicada
ax[1].plot(x, M_scaled, label="Momento Flector (Escalado)", color="blue")
ax[1].set_title("Diagrama de Momento Flector")
ax[1].set_xlabel("Posición a lo largo de la viga (m)")
ax[1].set_ylabel("Momento Flector (kN.m) Escalado")
ax[1].grid(True)

# Diagrama de fuerza axial
ax[2].plot(x, A_scaled, label="Fuerza Axial", color="green")
ax[2].set_title("Diagrama de Fuerza Axial")
ax[2].set_xlabel("Posición a lo largo de la viga (m)")
ax[2].set_ylabel("Fuerza Axial (kN)")
ax[2].grid(True)

plt.tight_layout()  # Ajusta el espaciado entre los gráficos
plt.show()
