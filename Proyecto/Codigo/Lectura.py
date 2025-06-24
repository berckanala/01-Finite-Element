import openseespy.opensees as ops
import meshio
import numpy as np

# Leer malla .msh
mesh = meshio.read("prueba.msh")

# Extraer nodos y tetraedros
nodes = mesh.points
tets = [cell.data for cell in mesh.cells if cell.type == "tetra"][0]

# Iniciar modelo 3D
ops.wipe()
ops.model("BasicBuilder", "-ndm", 3, "-ndf", 3)

# Crear nodos
for i, coord in enumerate(nodes):
    ops.node(i + 1, *coord[:3])  # OpenSees indexa desde 1

# Crear material
ops.nDMaterial("ElasticIsotropic", 1, 2.1e11, 0.3, 7850)

# Revisar primeros elementos para verificar conectividad
# Crear elementos tetraédricos
for i, conn in enumerate(tets):
    if len(conn) != 4 or len(set(conn)) != 4:
        print(f"❌ Elemento inválido {i}: {conn}")
        continue
    node_ids = [int(nid + 1) for nid in conn]  # 🔧 conversión a int puro
    ops.element("FourNodeTetrahedron", i + 1, *node_ids, 1)


import pyvista as pv
import numpy as np
import meshio

# Leer la malla nuevamente (por si no tienes aún los datos)
mesh = meshio.read("prueba.msh")
points = mesh.points
tets = [cell.data for cell in mesh.cells if cell.type == "tetra"][0]

# Crear celdas para PyVista (formato: [4, n0, n1, n2, n3])
# Cada fila debe comenzar con el número de nodos del elemento (4)
cells = np.hstack([np.full((len(tets), 1), 4), tets])

# Crear tipo de celda
celltypes = np.full(len(tets), pv.CellType.TETRA)

# Crear la malla unstructured
grid = pv.UnstructuredGrid(cells, celltypes, points)

# Visualizar
plotter = pv.Plotter()
plotter.add_mesh(grid, show_edges=True, color="lightblue", opacity=0.5)
plotter.add_axes()
plotter.show()
