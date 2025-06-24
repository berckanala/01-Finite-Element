import numpy as np
import matplotlib.pyplot as plt
from elementclass import Element
from nodeclass import Node
from structuresolver import StructureSolver
globalParameters = {'nDoF': 3}  
# Graficar la estructura y los diagramas de corte y momento
def plot_shear_and_moment_diagram(elemento, desplazamientos):
    # Nodo 1 (inicial) y nodo 2 (final) del elemento
    xi, yi = elemento.node_i.coordenadas
    xf, yf = elemento.node_j.coordenadas
    L = elemento.length

    # Cálculo de las fuerzas internas (corte y momento)
    dofs = np.hstack([elemento.node_i.idx, elemento.node_j.idx])
    K_local = elemento.k_local
    desplazamientos_local = desplazamientos[dofs]
    
    # Corte en los nodos
    V1 = -K_local[1, 0] * desplazamientos_local[0] - K_local[1, 1] * desplazamientos_local[1]  # Corte en nodo i
    V2 = -K_local[4, 3] * desplazamientos_local[3] - K_local[4, 4] * desplazamientos_local[4]  # Corte en nodo j

    # Momento en los nodos
    M1 = -K_local[2, 0] * desplazamientos_local[0] - K_local[2, 1] * desplazamientos_local[1]  # Momento en nodo i
    M2 = -K_local[5, 3] * desplazamientos_local[3] - K_local[5, 4] * desplazamientos_local[4]  # Momento en nodo j

    # Crear un conjunto de puntos para la viga
    x_local = np.linspace(0, L, 100)  # Interpolación a lo largo de la viga
    V_vals = (1 - x_local / L) * V1 + (x_local / L) * V2  # Interpolación lineal de corte

    # Integrar la fuerza cortante para obtener el diagrama de momento
    M_vals = np.cumsum(V_vals) * (x_local[1] - x_local[0])  # Aproximación de la integral
    M_vals -= M_vals[0]  # Ajustar para que M(x=0) = 0, o usar M1 si es necesario

    # Graficar el diagrama de corte y momento
    fig, ax = plt.subplots(2, 1, figsize=(8, 6))

    # Diagrama de corte
    ax[0].plot(x_local, V_vals, 'b-', label='Corte (V)')
    ax[0].set_title('Diagrama de Corte')
    ax[0].set_xlabel('Posición en la viga (m)')
    ax[0].set_ylabel('Fuerza Cortante (kN)')
    ax[0].grid(True)

    # Diagrama de momento
    ax[1].plot(x_local, M_vals, 'r-', label='Momento (M)')
    ax[1].set_title('Diagrama de Momentos')
    ax[1].set_xlabel('Posición en la viga (m)')
    ax[1].set_ylabel('Momento (kN·m)')
    ax[1].grid(True)

    plt.tight_layout()
    plt.show()

# Función principal para ejecutar el análisis con un solo elemento (viga)
def main():
    # Crear un nodo y un elemento (viga)
    nodo1 = Node(name=1, coordenadas=[0, 0], nodalLoad=[0, 0, 0], restrain=['r', 'r', 'f'])
    nodo2 = Node(name=2, coordenadas=[5, 0], nodalLoad=[0, 0, 0], restrain=['r', 'r', 'f'])
    nodo3 = Node(name=3, coordenadas=[5, 4], nodalLoad=[0, 0, 0], restrain=['f', 'f', 'f'])
    nodo4 = Node(name=4, coordenadas=[0, 4], nodalLoad=[0, 0, 0], restrain=['f', 'f', 'f'])

    # Crear un solo elemento (viga)
    elemento = Element(node_i=nodo1, node_j=nodo3, E=200e3, I=1.5, A=0.1)  # Viga entre nodo 1 y nodo 3

    # Ensamblar la matriz de rigidez global
    nDoF_total = len([nodo1, nodo2, nodo3, nodo4]) * 3
    K_global = np.zeros((nDoF_total, nDoF_total))
    dofs = np.hstack([nodo1.idx, nodo3.idx])
    for i in range(6):
        for j in range(6):
            K_global[dofs[i], dofs[j]] += elemento.k_local[i, j]

    # Crear el solucionador
    solver = StructureSolver([nodo1, nodo2, nodo3, nodo4], K_global, globalParameters['nDoF'])

    # Obtener los desplazamientos
    desplazamientos = solver.get_displacements()

    # Graficar los diagramas de corte y momento
    plot_shear_and_moment_diagram(elemento, desplazamientos)

# Ejecutar la función principal
main()
