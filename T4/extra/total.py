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

def reorder_triangle6(node_ids, nodes):
    """
    Asegura que los nodos del triángulo de 6 nodos estén en orden:
    [v1, v2, v3, m12, m23, m31]
    según convención FEM (Gmsh puede usar [v1, v2, v3, m23, m31, m12])
    """

    # Coordenadas de los 6 nodos
    coords = np.array([[nodes[nid - 1].x, nodes[nid - 1].y] for nid in node_ids])

    # Los primeros tres son los vértices
    v1, v2, v3 = coords[:3]

    # Los últimos tres son los nodos de borde (midside)
    mids = coords[3:]

    # Función auxiliar: distancia al segmento
    def is_between(p, a, b):
        mid = 0.5 * (a + b)
        return np.linalg.norm(p - mid) < 1e-8  # tolerancia

    # Identifica cuál nodo es m12, m23, m31
    m12 = m23 = m31 = None
    for i, p in enumerate(mids):
        if is_between(p, v1, v2):
            m12 = node_ids[3 + i]
        elif is_between(p, v2, v3):
            m23 = node_ids[3 + i]
        elif is_between(p, v3, v1):
            m31 = node_ids[3 + i]

    if None in (m12, m23, m31):
        print(f"⚠️ Error reordenando nodos LST: {node_ids}")
        return node_ids  # fallback

    # Devuelve orden correcto
    return [node_ids[0], node_ids[1], node_ids[2], m12, m23, m31]


from LST import LST
from Node import Node
# 👇 AQUI MISMO pega esta función:

def fix_LST_order(node_ids, nodes):
    """
    Corrige el orden de los nodos de un triángulo LST para asegurar sentido antihorario.
    Muestra advertencia si se aplica corrección.
    """
    p1 = np.array([nodes[node_ids[0] - 1].x, nodes[node_ids[0] - 1].y])
    p2 = np.array([nodes[node_ids[1] - 1].x, nodes[node_ids[1] - 1].y])
    p3 = np.array([nodes[node_ids[2] - 1].x, nodes[node_ids[2] - 1].y])
    detJ = 0.5 * ((p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]))

    if detJ <= 0:
        print(f"⚠️ Elemento con nodos {node_ids} tiene det(J) = {detJ:.4e} (corrigiendo orden)")
        return [node_ids[0], node_ids[2], node_ids[1], node_ids[5], node_ids[4], node_ids[3]]
    return node_ids

mesh = meshio.read(output_file)
nodes = [Node(i + 1, x, y) for i, (x, y, _) in enumerate(mesh.points)]

lst_elements = []
for cell_block in mesh.cells:
    if cell_block.type == "triangle6":
        for i, node_ids in enumerate(cell_block.data):
            node_ids = [int(id) + 1 for id in node_ids]  # pasar a base 1
            ordered_ids = reorder_triangle6(node_ids, nodes)
            lst_elements.append(LST(i + 1, ordered_ids))

        break

boundary_nodes = {1: set(), 2: set(), 3: set(), 4: set()}
if "line3" in mesh.cell_data_dict["gmsh:physical"]:
    for cell_block in mesh.cells:
        if cell_block.type == "line3":
            physical_ids = mesh.cell_data_dict["gmsh:physical"]["line3"]
            for line, phys_id in zip(cell_block.data, physical_ids):
                if phys_id in boundary_nodes:
                    for node_id in line:
                        boundary_nodes[phys_id].add(int(node_id) + 1)

for node in nodes:
    node.boundary_label = []
    for label_id, node_set in boundary_nodes.items():
        if node.id in node_set:
            node.boundary_label.append("Dirichlet Boundary")

# Verificación visual
for elem in lst_elements:
    coords = np.array([[nodes[i - 1].x, nodes[i - 1].y] for i in elem.node_ids])
    x, y = coords[:, 0], coords[:, 1]
    plt.plot(np.append(x[:3], x[0]), np.append(y[:3], y[0]), 'k-')
plt.axis("equal")
plt.show()


alpha = 5
for node in nodes:
    node.solve_u(alpha)


from scipy.interpolate import griddata
def plot_u_field(nodes, elements, use_attr="u", mode="colormap"):
    x = np.array([node.x for node in nodes])
    y = np.array([node.y for node in nodes])
    z = np.array([getattr(node, use_attr) for node in nodes])
    triangles = np.array([[nid - 1 for nid in elem.node_ids] for elem in elements])
    if mode == "colormap":
        xi = np.linspace(min(x), max(x), 100)
        yi = np.linspace(min(y), max(y), 100)
        X, Y = np.meshgrid(xi, yi)
        Z = griddata((x, y), z, (X, Y), method='cubic')
        plt.figure(figsize=(8, 6))
        plt.contourf(X, Y, Z, levels=30, cmap='viridis')
        plt.colorbar(label=use_attr)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(f"Mapa interpolado de {use_attr}(x,y)")
        plt.axis("equal")
        plt.tight_layout()
        plt.savefig(f"LST_{use_attr}_{mode}_plot.png", dpi=300)
        plt.show()

plot_u_field(nodes, lst_elements, use_attr="u", mode="colormap")

# Resolución FEM
Estructure = Pro_Solver(nodes, lst_elements, alpha)
Estructure.solve_system()

# Gráfico FEM
plot_u_field(nodes, lst_elements, use_attr="u_fem", mode="colormap")