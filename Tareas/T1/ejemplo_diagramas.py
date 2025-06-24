import numpy as np
import matplotlib.pyplot as plt
from elementclass import Element
from nodeclass import Node
from structuresolver import StructureSolver

def plot_structure_with_shear_diagram(elements):
    """
    Esta función recibe una lista de elementos y genera una figura con la estructura y el diagrama de corte.
    Se distinguen los elementos horizontales y verticales.
    """
    fig, ax = plt.subplots(figsize=(12, 8))  # Crear una figura de mayor tamaño

    # Graficar la estructura (vigas y columnas)
    for elem in elements:
        xi, yi = elem.node_i.coordenadas
        xj, yj = elem.node_j.coordenadas
        ax.plot([xi, xj], [yi, yj], 'b-', linewidth=2)  # Elemento en azul

    # Para cada elemento, calcular y graficar el diagrama de corte
    for elem in elements:
        xi, yi = elem.node_i.coordenadas
        xf, yf = elem.node_j.coordenadas

        # Llamar a la función de cálculo de fuerzas
        V_scaled, M_scaled, A_scaled = elem.calculate_forces()

        if xi == xf:  # Elemento vertical
            # Ajustar el corte para los elementos verticales (trasponer los valores de coordenadas)
            ax.plot(V_scaled, np.linspace(yi, yf, 100), label=f'Elemento {elem.node_i.name}-{elem.node_j.name} - Corte Vertical', linewidth=1)
        else:  # Elemento horizontal
            # Ajustar el corte para los elementos horizontales (añadir yf para mover el diagrama a la coordenada final)
            ax.plot(np.linspace(xi, xf, 100), V_scaled + yf, label=f'Elemento {elem.node_i.name}-{elem.node_j.name} - Corte Horizontal', linewidth=1)

    # Añadir etiquetas y formato
    ax.set_title("Diagrama de Corte de la Estructura")
    ax.set_xlabel("Posición a lo largo de la viga (m)")
    ax.set_ylabel("Esfuerzo Cortante (kN)")
    ax.grid(True)

    # Mostrar la leyenda
    ax.legend(loc='best')

    # Mostrar la figura
    plt.tight_layout()
    plt.show()

# Ejemplo de uso con algunos elementos
# Crear nodos
nodos = [
    Node(name=1, coordenadas=[0, 0], nodalLoad=[0, 0, 0], restrain=['r', 'r', 'f']),
    Node(name=2, coordenadas=[0, 4], nodalLoad=[0, 10, 0], restrain=['f', 'f', 'f']),
    Node(name=3, coordenadas=[5, 4], nodalLoad=[0, 10,0], restrain=['f', 'f', 'f']),
    Node(name=4, coordenadas=[5, 0], nodalLoad=[0, 0, 0], restrain=['r', 'r', 'f'])
]

# Crear elementos (vigas)
elements = [
    Element(node_i=nodos[0], node_j=nodos[1], E=200000, I=1.5, A=1),
    Element(node_i=nodos[1], node_j=nodos[2], E=200000, I=1.5, A=1),
    Element(node_i=nodos[2], node_j=nodos[3], E=200000, I=1.5, A=1)
]

# Llamar a la función para graficar la estructura con el diagrama de corte
plot_structure_with_shear_diagram(elements)