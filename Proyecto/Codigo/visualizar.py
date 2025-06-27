import meshio
import numpy as np
from collections import defaultdict
from openseespy import opensees as ops

# === Leer la malla
mesh = meshio.read("prueba_god.msh")
points = mesh.points

# Invertir field_data
field_data = mesh.field_data
tag_id_to_name = {v[0]: k for k, v in field_data.items()}

# Agrupar tetraedros
volumes = defaultdict(list)
tet_phys_tags = []

for cell_block, tag_block in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    if cell_block.type == "tetra":
        for conn, tag in zip(cell_block.data, tag_block):
            name = tag_id_to_name.get(tag, f"Tag_{tag}")
            volumes[name].append(conn)
            tet_phys_tags.append(tag)

# Nodos únicos por volumen
volume_nodes = {}
for name, tets in volumes.items():
    nodes_in_volume = np.unique(np.array(tets).flatten())
    volume_nodes[name] = nodes_in_volume

# Ver nodos
for name, nodes in volume_nodes.items():
    print(f"{name}: {len(nodes)} nodos")

print("hasta aquí todo bien")

import meshio
import numpy as np
from collections import defaultdict
from openseespy import opensees as ops

# === Leer malla ===
mesh = meshio.read("prueba_god.msh")
points = mesh.points

# Invertir field_data
field_data = mesh.field_data
tag_id_to_name = {v[0]: k for k, v in field_data.items()}

# Agrupar tetraedros
volumes = defaultdict(list)
tet_phys_tags = []

for cell_block, tag_block in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    if cell_block.type == "tetra":
        for conn, tag in zip(cell_block.data, tag_block):
            name = tag_id_to_name.get(tag, f"Tag_{tag}")
            volumes[name].append(conn)
            tet_phys_tags.append(tag)

# Nodos únicos por volumen
volume_nodes = {}
for name, tets in volumes.items():
    nodes_in_volume = np.unique(np.array(tets).flatten())
    volume_nodes[name] = nodes_in_volume

# Mostrar resumen
for name, nodes in volume_nodes.items():
    print(f"{name}: {len(nodes)} nodos")

print("hasta aquí todo bien")

# Crear modelo
ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)

# Crear nodos
for i, coord in enumerate(points):
    ops.node(i + 1, *coord)

# Materiales
ops.nDMaterial("ElasticIsotropic", 1, 3.5e9, 0.35, 1250)

ops.nDMaterial("ElasticIsotropic", 2, 1e20, 0, 7850)

# Crear elementos
element_id = 1
# Viga
for conn in volumes["Viga"]:
    node_tags = [int(n + 1) for n in conn]
    ops.element("FourNodeTetrahedron", element_id, *node_tags, 1)
    element_id += 1
# BC_R1
for conn in volumes["BC_R1"]:
    node_tags = [int(n + 1) for n in conn]
    ops.element("FourNodeTetrahedron", element_id, *node_tags, 2)
    element_id += 1
# BC_1
for conn in volumes["BC_1"]:
    node_tags = [int(n + 1) for n in conn]
    ops.element("FourNodeTetrahedron", element_id, *node_tags, 2)
    element_id += 1

# Restricciones
nodes_BC_R1 = volume_nodes["BC_R1"]
for n in nodes_BC_R1:
    ops.fix(int(n + 1), 1, 1, 1)

nodes_BC_1 = volume_nodes["BC_1"]
for n in nodes_BC_1:
    ops.fix(int(n + 1), 1, 1, 0)

# Nodos de la viga
nodes_viga = set(volume_nodes["Viga"])
# Nodos etiquetados como BC_1
nodes_bc1 = set(volume_nodes["BC_1"])

# Filtrar nodos de la viga que cumplan la condición
target_y = 110.0
target_z = 30.0
tolerance = 1e-3  # tolerancia por precisión numérica

# Nodos con y=110 y z=30
selected_nodes = []
for i, p in enumerate(points):
    if i in volume_nodes["Viga"]:
        if abs(p[1] - target_y) < tolerance and abs(p[2] - target_z) < tolerance:
            selected_nodes.append(i)

print(f"Se aplicará la carga en {len(selected_nodes)} nodos con y=110 y z=30.")

# Verificar si encontramos nodos
if len(selected_nodes) == 0:
    raise ValueError("No se encontraron nodos que cumplan la condición. Ajusta la tolerancia o revisa la malla.")

# Crear patrón de carga
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)

# Aplicar carga distribuida en esos nodos
for n in selected_nodes:
    ops.load(int(n + 1), 0, 0, -100000 / len(selected_nodes))



# Configuración análisis
ops.system("ProfileSPD")
ops.numberer("RCM")
ops.constraints("Plain")
ops.algorithm("Linear")
ops.integrator("LoadControl", 1.0)
ops.analysis("Static")

# Ejecutar análisis
print("Ejecutando análisis...")
ok = ops.analyze(1)
if ok == 0:
    print("Análisis exitoso.")
else:
    print("Error en el análisis.")

# Guardar desplazamientos por nodo
disp = {}
for i in range(points.shape[0]):
    u = ops.nodeDisp(i + 1)
    disp[i + 1] = u
    #print(f"Nodo {i + 1}: Ux={u[0]:.6e}, Uy={u[1]:.6e}, Uz={u[2]:.6e}")

import pyvista as pv

# === Construir la malla tetraédrica original ===
# Obtener todos los tetraedros en un solo array
all_tets = []
for name in ["Viga", "BC_R1", "BC_1"]:
    all_tets.extend(volumes[name])

cells = np.hstack([np.full((len(all_tets), 1), 4), np.array(all_tets)])
celltypes = np.full(len(all_tets), pv.CellType.TETRA)

# Crear grid
grid = pv.UnstructuredGrid(cells, celltypes, points)

# === Crear array de desplazamientos
displacements = np.zeros_like(points)
for i in range(points.shape[0]):
    u = ops.nodeDisp(i + 1)
    displacements[i, :] = u

# Escalar deformada (puedes ajustar)
scale_factor = 1000000.0

# === Crear nueva malla con desplazamiento aplicado
deformed_points = points + scale_factor * displacements
grid_deformed = pv.UnstructuredGrid(cells, celltypes, deformed_points)

# === Visualización
plotter = pv.Plotter()
plotter.add_mesh(grid, color="lightgray", opacity=0.3, show_edges=True, label="Original")
plotter.add_mesh(grid_deformed, color="red", opacity=0.6, show_edges=True, label="Deformada")
plotter.add_axes()
plotter.add_legend()
plotter.show()

