import numpy as np


coords = np.array([
    
    [0.0, 0.0],  #Nodo 0   
    [0.5, 0.0],  # Nodo 1
    [0.0, 0.5],  # Nodo 2
    [0.0, 1.0]   # Nodo 3
])


xi = 0
eta = 0


dN_dxi = 0.25 * np.array([
    -(1 - eta), (1 - eta), (1 + eta), -(1 + eta)
])

dN_deta = 0.25 * np.array([
    -(1 - xi), -(1 + xi), (1 + xi), (1 - xi)
])

# Jacobiano
J = np.zeros((2, 2))
for i in range(4):
    J[0, 0] += dN_dxi[i] * coords[i, 0]   
    J[0, 1] += dN_deta[i] * coords[i, 0]  
    J[1, 0] += dN_dxi[i] * coords[i, 1]   
    J[1, 1] += dN_deta[i] * coords[i, 1]  

# Imprimir Jacobiano y su determinante
print("Jacobiano:\n", J)
print("Determinante:", np.linalg.det(J))

