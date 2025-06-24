import pyvista as pv
import meshio
import numpy as np

# Leer malla
mesh = meshio.read("prueba.msh")
points = mesh.points
tets = [cell.data for cell in mesh.cells if cell.type == "tetra"][0]

# Formato para PyVista
cells = np.hstack([np.full((len(tets), 1), 4), tets])
celltypes = np.full(len(tets), pv.CellType.TETRA)
grid = pv.UnstructuredGrid(cells, celltypes, points)

# Visualización 3D en ventana externa
plotter = pv.Plotter()
plotter.add_mesh(grid, show_edges=True, color="lightblue", opacity=0.5)
plotter.add_axes()
plotter.show()
