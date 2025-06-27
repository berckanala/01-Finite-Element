import meshio
import pyvista as pv
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt  # Para obtener colores del cmap

# === Leer el archivo .msh ===
mesh = meshio.read("Pro_malla.msh")
points = mesh.points

# === Invertir field_data
field_data = mesh.field_data
tag_id_to_name = {v[0]: k for k, v in field_data.items()}

# === Agrupar tetraedros por volumen físico
volumes = defaultdict(list)
phys_tags = []

for cell_block, tag_block in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    if cell_block.type == "tetra":
        for conn, tag in zip(cell_block.data, tag_block):
            name = tag_id_to_name.get(tag, f"Tag_{tag}")
            volumes[name].append(conn)
            phys_tags.append(tag)

# === Preparar datos para PyVista
all_tets = []
tet_colors = []
tag_to_color = {}
color_map = plt.get_cmap("Set3")
unique_tags = sorted(set(phys_tags))

for i, tag in enumerate(unique_tags):
    tag_to_color[tag] = i

for conn, tag in zip([tet for conn_list in volumes.values() for tet in conn_list], phys_tags):
    all_tets.append(conn)
    tet_colors.append(tag_to_color[tag])

cells = np.hstack([np.full((len(all_tets), 1), 4), np.array(all_tets)])
celltypes = np.full(len(all_tets), pv.CellType.TETRA)

grid = pv.UnstructuredGrid(cells, celltypes, points)
grid.cell_data["volumen_fisico"] = np.array(tet_colors)

# === Crear colores válidos para la leyenda ===
legend_entries = []
for i, tag in enumerate(unique_tags):
    name = tag_id_to_name[tag]
    rgb = color_map(i / len(unique_tags))[:3]  # Obtener color como RGB
    legend_entries.append([name, rgb])

# === Visualizar
plotter = pv.Plotter()
plotter.add_mesh(grid, scalars="volumen_fisico", show_edges=True, cmap="Set2", opacity=0.6)
plotter.add_legend(legend_entries)
plotter.add_axes()
plotter.show()
