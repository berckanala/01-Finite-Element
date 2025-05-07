import numpy as np
from scipy.optimize import fsolve

# Coordenadas del punto físico
P = np.array([4.0, 4.0])

# Coordenadas de los nodos del Q4 en el espacio físico (u_1)
# Formato: [x1, y1, x2, y2, ..., x4, y4]
u_1 = [3, 0, 9, 0, 0, 9, 0, 3]
coords = np.array(u_1).reshape(4, 2)

# Función que define el sistema de ecuaciones
def sistema(vars):
    xi, eta = vars

    # Funciones de forma para Q4
    N1 = 0.25 * (1 - xi) * (1 - eta)
    N2 = 0.25 * (1 + xi) * (1 - eta)
    N3 = 0.25 * (1 + xi) * (1 + eta)
    N4 = 0.25 * (1 - xi) * (1 + eta)

    # Interpolación de la posición física x(ξ, η)
    N = np.array([N1, N2, N3, N4])
    x_eta = N @ coords  # x_eta es el punto interpolado

    return x_eta - P  # Queremos que sea igual a P

# Suposición inicial
xi_eta_0 = [0.0, 0.0]

# Resolver
sol = fsolve(sistema, xi_eta_0)
print("Question A:")
print(f"Solución: ξ = {sol[0]:.4f}, η = {sol[1]:.4f}")
print("")

# Función de interpolación isoparamétrica
def interpolar_xy(xi, eta, coords):
    N1 = 0.25 * (1 - xi) * (1 - eta)
    N2 = 0.25 * (1 + xi) * (1 - eta)
    N3 = 0.25 * (1 + xi) * (1 + eta)
    N4 = 0.25 * (1 - xi) * (1 + eta)
    
    N = np.array([N1, N2, N3, N4])
    xy = N @ coords
    return xy

# Ejemplo: valores dados de xi y eta
xi = 2*((7-1)/(30-1))-1
eta = 2*((11-1)/(11))-1

xy = interpolar_xy(xi, eta, coords)
print("Question B:")
print(f"Coordenadas físicas: x = {xy[0]:.4f}, y = {xy[1]:.4f}")
print("")
import numpy as np

# Coordenadas nodales del Quad4 (debes ajustar si las reales son distintas)
# Nodo 1, 2, 3, 4 en sentido antihorario
coords = np.array([
    [0.0, 3.0],   # Nodo 1
    [3.0, 0.0],   # Nodo 2
    [9.0, 0.0],   # Nodo 3
    [0.0, 9.0]    # Nodo 4
])
coords_desplazado = np.array([
    [0.0, 1.0],   # Nodo 1
    [5.0, 0.0],   # Nodo 2
    [12.0, 0.0],   # Nodo 3
    [0.0, 6.0]    # Nodo 4
])
# Desplazamientos globales u = [ux1, uy1, ux2, uy2, ux3, uy3, ux4, uy4]
u = np.array([3, 0, 9, 0, 0, 9, 0, 3])

# Punto de evaluación en coordenadas naturales
xi = -0.5862
eta = 0.8182

# Derivadas de las shape functions respecto a xi y eta
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

# Inversa del Jacobiano
J_inv = np.linalg.inv(J)

# Derivadas respecto a x, y
dN_dx = np.zeros(4)
dN_dy = np.zeros(4)

for i in range(4):
    dN = np.array([dN_dxi[i], dN_deta[i]])
    dN_global = J_inv @ dN
    dN_dx[i] = dN_global[0]
    dN_dy[i] = dN_global[1]

# Construcción de la matriz B (deformación plana: εx, εy, γxy)
B = np.zeros((3, 8))
for i in range(4):
    B[0, 2*i]     = dN_dx[i]
    B[1, 2*i + 1] = dN_dy[i]
    B[2, 2*i]     = dN_dy[i]
    B[2, 2*i + 1] = dN_dx[i]

# Deformación en el punto
strain = B @ u

print("Question C:")

print("Matriz B:")
print(B)
print("")

print("Question D:")
print(strain)
print("")


print("Question E:")
#Explain what Iso-Parameteric finite-elements are and why they are attractive
print("Iso-Parametric finite-elements are a type of finite element method (FEM) where the shape functions used to interpolate the geometry of the element are the same as those used to interpolate the field variables (like displacements).")