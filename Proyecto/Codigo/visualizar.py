import meshio
import numpy as np
from collections import defaultdict
from openseespy import opensees as ops
import pyvista as pv

# === Leer la malla
# === Leer la malla
mesh = meshio.read("prueba_god.msh")
points = mesh.points

print(f"=== Nodos totales en la malla: {points.shape[0]} ===")

# Invertir field_data
field_data = mesh.field_data
tag_id_to_name = {v[0]: k for k, v in field_data.items()}

# Agrupar tetraedros
volumes = defaultdict(list)
tet_phys_tags = []

# Mostrar todos los bloques de celdas
print("\n=== Tipos de celdas encontradas ===")
for c in mesh.cells:
    print(f"Celda: {c.type}, cantidad: {len(c.data)}")

# Extraer tetraedros con Physical Volume
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

# Mostrar resumen de nodos por volumen
print("\n=== Nodos por volumen ===")
for name, nodes in volume_nodes.items():
    print(f"{name}: {len(nodes)} nodos")

# Nodos usados en todos los tetraedros
all_used_nodes = np.unique(np.vstack([v for v in volumes.values()]).flatten())
print(f"\n=== Nodos usados en tetraedros: {len(all_used_nodes)}")

# Avisar si hay nodos no usados
n_unused = points.shape[0] - len(all_used_nodes)
if n_unused > 0:
    print(f"⚠️ Atención: {n_unused} nodos no están en ningún tetraedro (posible malla incompleta o volúmenes sin etiqueta).")
else:
    print("✅ Todos los nodos están usados en tetraedros.")


print("=== Preparando modelo ===")

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
viga_ids = []
# Viga
for conn in volumes["Viga"]:
    node_tags = [int(n + 1) for n in conn]
    ops.element("FourNodeTetrahedron", element_id, *node_tags, 1)
    viga_ids.append(element_id)
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
    ops.fix(int(n + 1), 1, 0, 1)

# Filtrar nodos de la viga donde se aplica la carga
target_x = 110.0
target_y = 30.0
tolerance = 1e-3

selected_nodes = []
for i, p in enumerate(points):
    if i in volume_nodes["Viga"]:
        if abs(p[1] - target_y) < tolerance and abs(p[0] - target_x) < tolerance:
            selected_nodes.append(i)

print(f"Se aplicará la carga en {len(selected_nodes)} nodos con x=110 y y=30.")

if len(selected_nodes) == 0:
    raise ValueError("No se encontraron nodos que cumplan la condición.")

# Crear patrón de carga
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)

# Aplicar carga distribuida
for n in selected_nodes:
    ops.load(int(n + 1), 0,-1e10/ len(selected_nodes),0 )

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
    print("✅ Análisis exitoso.")
else:
    raise RuntimeError("❌ Error en el análisis.")

# Guardar desplazamientos
displacements = np.zeros_like(points)
for i in range(points.shape[0]):
    u = ops.nodeDisp(i + 1)
    displacements[i, :] = u

max_disp = np.linalg.norm(displacements, axis=1).max()
print(f"Desplazamiento máximo: {max_disp:.6e} m")
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
scale_factor = 1.0

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

von_mises = []

print("=== Calculando tensiones de von Mises con 'stresses' ===")

for eid in viga_ids:
    stress = ops.eleResponse(eid, "stresses")
    if not stress or len(stress) < 6:
        print(f"⚠️ Elemento {eid} sin tensiones, asignando 0.")
        von_mises.append(0.0)
        continue
    sx, sy, sz, t_yz, t_xz, t_xy = stress
    vm = np.sqrt(
        0.5*((sx - sy)**2 + (sy - sz)**2 + (sz - sx)**2) + 3*(t_xy**2 + t_yz**2 + t_xz**2)
    )
    von_mises.append(vm)

print(f"Tensiones von Mises calculadas en {len(von_mises)} elementos.")

cells_viga = np.hstack([np.full((len(viga_ids), 1), 4), np.array(volumes["Viga"])])
celltypes_viga = np.full(len(viga_ids), pv.CellType.TETRA)
grid_viga = pv.UnstructuredGrid(cells_viga, celltypes_viga, points)
grid_viga.cell_data["von_mises"] = np.array(von_mises)

plotter_vm = pv.Plotter()
plotter_vm.add_mesh(
    grid_viga,
    scalars="von_mises",
    cmap="viridis",
    show_edges=True,
    opacity=0.7,
    scalar_bar_args={"title": "Von Mises [Pa]"}
)
plotter_vm.add_axes()
plotter_vm.add_title("Distribución de tensiones von Mises")
plotter_vm.show()
