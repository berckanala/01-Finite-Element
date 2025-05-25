import numpy as np
import matplotlib.pyplot as plt
import gmsh
import os
import meshio
import math
import sys
from matplotlib.patches import Polygon

from Quad4 import Quad4
from Node import Node
from units import mm, cm, m, kgf, N, tf, kN, MPa, GPa
from fem import Material, Membrane
from functions import matrix_extract,  get_nodes_from_physical_id, get_line_load_global_vector, matrix_replace


geo_file = "T3.geo"
lc = 50

gmsh.initialize()
gmsh.open(geo_file)

gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.Algorithm", 8)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)

output_file = "malla_quad4.msh"
gmsh.write(output_file)
gmsh.fltk.run()
gmsh.finalize()

ASTM_36 = Material.Material(name="ASTM 36", E=200.0 * GPa, nu=0.30, rho=7850 * kgf / m**3)
Sup = Membrane.Membrane(name="Steel", thickness=2 * cm, material=ASTM_36)
section_dict = {"Steel": Sup}

mesh = meshio.read(output_file)
tag_to_name = {v[0]: k for k, v in mesh.field_data.items()}

mesh = meshio.read(output_file)
tag_to_name = {v[0]: k for k, v in mesh.field_data.items()}

# === CONSTRUCCIÓN DE GRUPOS DE NODOS ===
grupos = {}

def crear_nodo(node_id, coords, restrain=None):
    return Node(node_id, coords, restrain=restrain)

for cell_block, phys_tags in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    for cell, tag in zip(cell_block.data, phys_tags):
        nombre = tag_to_name.get(tag, str(tag))
        if nombre not in grupos:
            grupos[nombre] = []
        if cell_block.type == "quad":
            for node_id in cell:
                x, y = mesh.points[node_id][:2]
                grupos[nombre].append(crear_nodo(node_id, (x, y)))
        elif cell_block.type in ["line", "vertex"]:
            for node_id in cell:
                x, y = mesh.points[node_id][:2]
                restrain = ["r", "r"] if nombre in ["BC_R1", "BC_R2"] else ["f", "f"]
                grupos[nombre].append(crear_nodo(node_id, (x, y), restrain))



# Eliminar duplicados por grupo
for nombre in grupos:
    nodos_unicos = {}
    for n in grupos[nombre]:
        nodos_unicos[n.name] = n
    grupos[nombre] = list(nodos_unicos.values())
         


nodes_dict = {}
for group in grupos:
    for node in grupos[group]:
        nodes_dict[node.name] = node

def ordenar_nodos_quad(nodos):
    coords = np.array([n.coordenadas for n in nodos])
    cx, cy = np.mean(coords, axis=0)
    angulos = np.arctan2(coords[:, 1] - cy, coords[:, 0] - cx)
    orden = np.argsort(angulos)
    return [nodos[i] for i in orden]

quads = mesh.cells_dict["quad"]
tags = mesh.cell_data_dict["gmsh:physical"]["quad"]
elements = []
nodes = set()

for i in range(len(tags)):
    tag = tags[i]
    group_name = tag_to_name[tag]
    material = section_dict[group_name]
    node_ids = quads[i]

    nodo_a = nodes_dict[node_ids[0]]
    nodo_b = nodes_dict[node_ids[1]]
    nodo_c = nodes_dict[node_ids[2]]
    nodo_d = nodes_dict[node_ids[3]]

    for nodo in [nodo_a, nodo_b, nodo_c, nodo_d]:
        nodes.add(nodo)
    nodos = ordenar_nodos_quad([nodo_a, nodo_b, nodo_c, nodo_d])
    elem = Quad4(i + 1, nodos, section=material)
    elements.append(elem)

nodes = list(set(nodes) | set(grupos.get("BC_R1", [])) | set(grupos.get("BC_1", [])) | set(grupos.get("otros_grupos", [])))
nDoF = 2  # Por ejemplo, 2 grados de libertad por nodo: desplazamiento en x e y

print(f"Total de nodos: {len(nodes)}")

todos_los_nodos =  list(set(nodes) | set(grupos.get("BC_R1", [])) | set(grupos.get("BC_1", [])) | set(grupos.get("otros_grupos", [])))

x = np.array([node.coordenadas[0] for node in todos_los_nodos])
y = np.array([node.coordenadas[1] for node in todos_los_nodos])
node_index_map = {node.name: i for i, node in enumerate(todos_los_nodos)}

fig, ax = plt.subplots(figsize=(8, 6))

# Dibujar elementos
for elem in elements:
    indices = [node_index_map[node.name] for node in elem.node_list]
    coords = np.array([[x[i], y[i]] for i in indices + [indices[0]]])
    ax.plot(coords[:, 0], coords[:, 1], "k-", linewidth=1)

# Dibujar nodos Force
# Dibujar grupo BC_1 en verde si existe
if "BC_1" in grupos:
    for i, node in enumerate(grupos["BC_1"]):
        if node.name in node_index_map:
            idx = node_index_map[node.name]
            ax.plot(x[idx], y[idx], "ro", markersize=6, label="BC_1" if i == 0 else "")


# Dibujar nodos Nut
# Dibujar apoyos en azul desde los grupos físicos "BC_R1" y "BC_R2"
for grupo in ["BC_R1"]:
    if grupo in grupos:
        for i, node in enumerate(grupos[grupo]):
            if node.name in node_index_map:
                idx = node_index_map[node.name]
                ax.plot(x[idx], y[idx], 'bo', markersize=6, label='Restricción' if i == 0 else '')
                
                print(node, node.name, node.restrain, node.idx)




# Finalizar gráfico
ax.set_title("Estructura con elementos y nodos marcados")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")
ax.legend(
    handles=[
        plt.Line2D([0], [0], color="k", lw=1, label="Element"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="r", markersize=6, label="Force"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="b", markersize=6, label="Fixed Support"),
    ],
    loc="center left",
    bbox_to_anchor=(1.0, 0.5)
)
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

#Hasta acá funciona todo perfecto

def solve(nodes, elements):
    """
    Resuelve el sistema FEM considerando cargas distribuidas en líneas.

    Args:
        nodes (list): Lista de nodos con .idx, .nodalLoad, .restrain
        elements (list): Lista de elementos con .kg, .idx
        mesh: Objeto meshio con malla y etiquetas físicas
        load_dictionary (dict): Dict {id_físico: valor_carga}

    Returns:
        u (np.ndarray): Desplazamientos globales
        F_reactions (np.ndarray): Reacciones en DOFs restringidos
        Fj (np.ndarray): Fuerzas nodales aplicadas
        F_total (np.ndarray): Fuerzas totales incluyendo reacciones
    """
    nNodes = len(nodes)
    system_nDof = 2 * nNodes  # 2 DOFs por nodo

    # === CARGAS DISTRIBUIDAS EN LÍNEAS ===
    F_line_global = np.zeros(system_nDof)
    Fj = np.zeros(system_nDof)
    for node in nodes:
        for dof_idx, force in zip(node.idx, node.nodalLoad):
            
            Fj[dof_idx] = force
    
    


    # === MATRIZ DE RIGIDEZ GLOBAL ===
    Kg = np.zeros((system_nDof, system_nDof))
    for element in elements:
        Kg = matrix_replace(Kg, element.kg, element.idx, element.idx)


    # === RESTRICCIONES ===
    nodeIndex = np.full(system_nDof, '', dtype=str)
    for node in nodes:
        for dof_idx, status in zip(node.idx, node.restrain):
            nodeIndex[dof_idx] = status

    freeIndices = np.where(nodeIndex == 'f')[0]
    restrainedIndices = np.where(nodeIndex == 'r')[0]
    print("Índices restringidos: ", restrainedIndices)
    print()
    


    # === PARTICIÓN DEL SISTEMA ===
    Kff = matrix_extract(Kg, freeIndices, freeIndices)
    Kfr = matrix_extract(Kg, freeIndices, restrainedIndices)
    Krf = matrix_extract(Kg, restrainedIndices, freeIndices)
    Krr = matrix_extract(Kg, restrainedIndices, restrainedIndices)

    # === FUERZA TOTAL (aplicada + distribuida) ===
    F = Fj + F_line_global
    Ff = F[freeIndices]
    Fr = F[restrainedIndices]

    # === SOLUCIÓN ===
    uf = np.linalg.solve(Kff, Ff - Kfr @ Fr)
    ur = np.zeros(len(restrainedIndices))

    u = np.zeros(system_nDof)
    u[freeIndices] = uf
    u[restrainedIndices] = ur
    


    # === REACCIONES ===
    Fr = Krf @ uf + Krr @ ur
    F_reactions = np.zeros(system_nDof)
    F_reactions[restrainedIndices] = Fr
    


    return u


def apply_distributed_force(grupo_nodos, fuerza_total):
    """
    Aplica una fuerza distribuida a lo largo de una línea definida por nodos.
    La fuerza total se reparte proporcionalmente según la longitud de cada segmento.

    Args:
        grupo_nodos (list): Lista de nodos (Node) ordenados a lo largo de la línea.
        fuerza_total (float): Fuerza total a repartir (en N o unidades consistentes).
                              Puede estar en cualquier dirección (por ejemplo, solo vertical).
    """
    nodos = grupo_nodos
    n = len(nodos)

    if n < 2:
        print("Se requieren al menos dos nodos para aplicar fuerza distribuida.")
        return

    # Paso 1: calcular longitud total y longitudes de tramos
    longitudes = []
    total_length = 0.0

    for i in range(n - 1):
        dx = nodos[i+1].coordenadas[0] - nodos[i].coordenadas[0]
        dy = nodos[i+1].coordenadas[1] - nodos[i].coordenadas[1]
        L = np.sqrt(dx**2 + dy**2)
        longitudes.append(L)
        total_length += L

    if total_length == 0:
        print("La longitud total es cero, no se puede aplicar la fuerza.")
        return

    # Paso 2: calcular fuerza distribuida por unidad de longitud
    q_lineal = fuerza_total / total_length  # N/m

    # Paso 3: inicializar fuerzas nodales acumuladas
    nodal_forces = {node.name: np.array([0.0, 0.0]) for node in nodos}

    # Paso 4: recorrer cada segmento y distribuir fuerza
    for i in range(n - 1):
        ni = nodos[i]
        nj = nodos[i + 1]
        xi, yi = ni.coordenadas
        xj, yj = nj.coordenadas

        dx = xj - xi
        dy = yj - yi
        L = longitudes[i]

        # Vector unitario del tramo
        tx = dx / L
        ty = dy / L

        # Vector normal (perpendicular)
        nx = 1
        ny = 0  # sentido hacia abajo si horizontal

        # Fuerza total en el tramo
        F_total = q_lineal * L

        # Componentes globales de la fuerza
        fx = F_total * nx
        fy = F_total * ny

        # Distribuir mitad a cada nodo del segmento
        nodal_forces[ni.name] += np.array([fx / 2, fy / 2])
        nodal_forces[nj.name] += np.array([fx / 2, fy / 2])

    # Paso 5: asignar fuerzas a cada nodo
    for node in nodos:
        fx, fy = nodal_forces[node.name]
        node.set_nodalLoad([fx, fy])

def plot_deformed_shape(nodes, elements, u, contador, scale, folder="figuras"):
    """
    Dibuja la estructura deformada y original, marcando nodos cargados y restringidos.

    Args:
        nodes (list): lista de nodos con .coordenadas, .name, .idx
        elements (list): lista de elementos con .node_list
        u (np.ndarray): vector de desplazamientos global
        contador (int): índice del caso (0, 1, 2)
        scale (float): factor de amplificación de la deformación
    """
    filename = None
    # Inicializar vectores de coordenadas
    n = len(nodes)
    x = np.zeros(n)
    y = np.zeros(n)
    x_def = np.zeros(n)
    y_def = np.zeros(n)

    # Construcción coordenadas originales y deformadas usando node.idx
    for i, node in enumerate(nodes):
        x[i] = node.coordenadas[0]
        y[i] = node.coordenadas[1]
        x_def[i] = x[i] + scale * u[node.idx[0]]
        y_def[i] = y[i] + scale * u[node.idx[1]]

    # Mapeo name → índice en array
    node_index_map = {node.name: i for i, node in enumerate(nodes)}


    # Gráfico
    fig, ax = plt.subplots(figsize=(8, 6))

    # Estructura original (negro)
    for elem in elements:
        indices = [node_index_map[node.name] for node in elem.node_list]
        coords = np.array([[x[i], y[i]] for i in indices + [indices[0]]])
        ax.plot(coords[:, 0], coords[:, 1], 'k-', linewidth=1)

    # Estructura deformada (rojo punteado)
    for elem in elements:
        indices = [node_index_map[node.name] for node in elem.node_list]
        coords = np.array([[x_def[i], y_def[i]] for i in indices + [indices[0]]])
        ax.plot(coords[:, 0], coords[:, 1], 'r--', linewidth=1)



    # Título dinámico por caso
    if contador == 0:
        title="Deformed sctructure - Case 1 low resolution"
        ax.set_title(title)
    elif contador == 1:
        title="Deformed sctructure - Case 2 low resolution"
        ax.set_title(title)
    elif contador == 2:
        title="Deformed sctructure - Case 3 low resolution"
        ax.set_title(title)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")

    # Leyenda fuera del gráfico
    ax.legend(
        handles=[
            plt.Line2D([0], [0], color='k', lw=1, label='Elemento'),
            plt.Line2D([0], [0], linestyle='--', color='r', lw=1, label='Deformada')
        ],
        loc="center left",
        bbox_to_anchor=(1.02, 0.5)
    )

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()
    

nodos_fuerza = grupos.get("BC_1", [])
apply_distributed_force(nodos_fuerza, 1500 * kN)
u = solve(nodes, elements)
plot_deformed_shape(nodes, elements, u, 0, scale=1000, folder="Resultados")
