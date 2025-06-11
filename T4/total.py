import numpy as np
import matplotlib.pyplot as plt
import gmsh
import os
import meshio
import math
import sys
from matplotlib.patches import Polygon
from matplotlib.tri import Triangulation
from LST import LST
from Node import Node
from solve import Pro_Solver
from functions import matrix_extract,  get_nodes_from_physical_id, get_line_load_global_vector, matrix_replace
input_file = "T4.geo"
output_file = "mallas/malla_LST.msh"  
titulo = "Caso LST"
def fix_LST_order(node_ids, nodes):
    """
    Corrige el orden de los nodos de un triángulo LST para asegurar sentido antihorario.
    Muestra advertencia si se aplica corrección.
    """
    # Obtener coordenadas de los 3 nodos vértice
    p1 = np.array([nodes[node_ids[0] - 1].x, nodes[node_ids[0] - 1].y])
    p2 = np.array([nodes[node_ids[1] - 1].x, nodes[node_ids[1] - 1].y])
    p3 = np.array([nodes[node_ids[2] - 1].x, nodes[node_ids[2] - 1].y])

    # Calcular área orientada
    area = 0.5 * ((p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]))

    if area > 0:
        return node_ids  # OK
    else:
        # Invertir orden vértices y medios
        print(f"⚠️ Reordenando nodos LST (ID vértices: {node_ids[:3]})")
        return [node_ids[0], node_ids[2], node_ids[1],  # vértices
                node_ids[5], node_ids[4], node_ids[3]]  # medios reordenados
import meshio

# Leer archivo de malla
mesh = meshio.read(output_file)

# Crear nodos desde las coordenadas del archivo (asumiendo 2D)
nodes = [Node(i + 1, x, y) for i, (x, y, _) in enumerate(mesh.points)]

# Crear elementos LST (triángulo de 6 nodos)
# Crear nodos desde las coordenadas del archivo (asumiendo 2D)
nodes = [Node(i + 1, x, y) for i, (x, y, _) in enumerate(mesh.points)]

# Crear elementos LST (triángulo de 6 nodos)
lst_elements = []
for cell_block in mesh.cells:
    if cell_block.type == "triangle6":
        for i, elem in enumerate(cell_block.data):
            node_ids = [int(id) + 1 for id in elem]  # base 1
            fixed_ids = fix_LST_order(node_ids, nodes)
            lst_elements.append(LST(i + 1, fixed_ids))
        break
  # Ya procesamos los triángulos, no hace falta seguir

# Identificar nodos en los bordes físicos (Physical Groups de Gmsh)
boundary_nodes = {1: set(), 2: set(), 3: set(), 4: set()}  # etiquetas de borde

for cell_block in mesh.cells:
    if cell_block.type == "line3":
        physical_ids = mesh.cell_data_dict["gmsh:physical"]["line3"]
        for line, phys_id in zip(cell_block.data, physical_ids):
            if phys_id in boundary_nodes:
                for node_id in line:
                    boundary_nodes[phys_id].add(int(node_id) + 1)  # Convertir a base 1

# Asignar etiquetas de borde a cada nodo
for node in nodes:
    node.boundary_label = []
    for label_id, node_set in boundary_nodes.items():
        if node.id in node_set:
            node.boundary_label.append("Dirichlet Boundary")

# Lista final de elementos
elements = lst_elements
import matplotlib.pyplot as plt

for elem in lst_elements:
    coords = np.array([[nodes[i - 1].x, nodes[i - 1].y] for i in elem.node_ids])
    x, y = coords[:, 0], coords[:, 1]
    plt.plot(np.append(x[:3], x[0]), np.append(y[:3], y[0]), 'k-')  # vértices
plt.axis("equal")
plt.show()



#Ahora debo calcular la solucion u
alpha = 3

for node in nodes:
    node.solve_u(alpha)


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

def plot_u_field(nodes, elements, use_attr="u", mode="colormap"):
    """
    Visualiza el campo u(x,y) o u_fem(x,y) usando distintos modos:
    - 'colormap'  → mapa de color interpolado (2D)
    """
    # Extraer coordenadas y valores de campo
    x = np.array([node.x for node in nodes])
    y = np.array([node.y for node in nodes])
    z = np.array([getattr(node, use_attr) for node in nodes])

    # Conectividad de elementos (base 1 → base 0)
    triangles = np.array([[nid - 1 for nid in elem.node_ids] for elem in elements])
    #tri = Triangulation(x, y, triangles)

    if mode == "colormap":
        # Crear malla regular
        xi = np.linspace(min(x), max(x), 100)
        yi = np.linspace(min(y), max(y), 100)
        X, Y = np.meshgrid(xi, yi)

        # Interpolación desde nodos hacia malla
        Z = griddata((x, y), z, (X, Y), method='cubic')

        # Graficar mapa de colores
        plt.figure(figsize=(8, 6))
        plt.contourf(X, Y, Z, levels=30, cmap='viridis')
        plt.colorbar(label=use_attr)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(f"Mapa interpolado de {use_attr}(x,y)")
        plt.axis("equal")
        plt.tight_layout()
        plt.savefig(f"CST_{use_attr}_{mode}_plot.png", dpi=300)
        plt.show()

    else:
        raise ValueError("Modo no válido. Usa solo 'colormap' en esta versión.")

# Uso: solución exacta u
plot_u_field(nodes, elements, use_attr="u", mode="colormap")
free_nodes = [node.id for node in nodes if not node.boundary_label]
print(f"Nodos libres: {len(free_nodes)} → {free_nodes}")
dirichlet_nodes = []
for node in nodes:
    if hasattr(node, "boundary_label") and any("Dirichlet" in label for label in node.boundary_label):
        dirichlet_nodes.append((node.id, node.x, node.y, node.u))

print(f"🧱 Total Dirichlet nodes: {len(dirichlet_nodes)}")
for nid, x, y, u in dirichlet_nodes:
    print(f"Nodo {nid} → ({x:.3f}, {y:.3f}) → u = {u}")
Estructure = Pro_Solver(nodes, elements, alpha)
Estructure.solve_system()

plot_u_field(nodes, elements, use_attr="u_fem", mode="colormap")