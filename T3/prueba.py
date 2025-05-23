geo_file = "T3P.geo"
lc = 2000

gmsh.initialize()
gmsh.open(geo_file)

# Forzar elementos cuadriláteros
gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.Algorithm", 8)  # Frontal recombination (opcional, mejora calidad)

# Establecer tamaño de los elementos
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)

# Sincronizar modelo
gmsh.model.geo.synchronize()

# Generar malla 2D
gmsh.model.mesh.generate(2)

# Guardar la malla
output_file = 'malla_quad4.msh'
gmsh.write(output_file)

# Visualizar
gmsh.fltk.run()
gmsh.finalize()

ASTM_36=Material.Material(name='ASTM 36',
             E=200.0*GPa,
             nu=0.30,
             rho=7850*kgf/m**3)

# Definimos una seccion de la membrana
Sup=Membrane.Membrane(name='Steel',
              thickness=2*cm,
              material=ASTM_36)


# Definimos los grupos fisicos de las partes del modelo
# Map the physical group id to a section
section_dict={"Steel":Sup}

mesh = meshio.read(output_file)

for i, (cell_block, phys_ids) in enumerate(zip(mesh.cells, mesh.cell_data["gmsh:physical"])):
    print(f"Block {i} - Tipo: {cell_block.type}, Cantidad: {len(cell_block.data)}, Physical tags: {set(phys_ids)}")


mesh = meshio.read(output_file)

# Asociación de tag físico con nombre
tag_to_name = {v[0]: k for k, v in mesh.field_data.items()}

# Diccionario {nombre_grupo: [Node, Node, ...]}
grupos = {}
Materials = {}
# Procesar elementos tipo triangle
for cell_block, phys_tags in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    if cell_block.type != "quad":
        continue
    for tri, tag in zip(cell_block.data, phys_tags):
        nombre = tag_to_name.get(tag, f"{tag}")

        if nombre not in grupos:
            grupos[nombre] = []
        for node_id in tri:
            x, y = mesh.points[node_id][:2]
            grupos[nombre].append(Node(node_id+1, (x, y)))


# Procesar elementos tipo line (para grupos como "Fuerza")
for cell_block, phys_tags in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    if cell_block.type != "line":
        continue
    for line, tag in zip(cell_block.data, phys_tags):
        nombre = tag_to_name.get(tag, f"{tag}")
        if nombre not in grupos:
            grupos[nombre] = []
        for node_id in line:
            x, y = mesh.points[node_id][:2]
            restrain = ["f", "f"]
            if nombre in ["Nut"]:
                print(f"Grupo {nombre} - Nodo {node_id+1} - Coordenadas ({x}, {y})")
                restrain = ["r", "r"]
            grupos[nombre].append(Node(node_id+1, (x, y), restrain=restrain))


   


# Eliminar nodos duplicados por grupo
for nombre in grupos:
    nodos_unicos = {}
    for n in grupos[nombre]:
        nodos_unicos[n.name] = n
    grupos[nombre] = list(nodos_unicos.values())
    print(f"Grupo {nombre} - Cantidad de nodos: {len(grupos[nombre])}")

nodes_dict = {}
for group in grupos:
    for node in grupos[group]:
        nodes_dict[node.name] = node

def ordenar_nodos_quad(nodos):
    """
    Ordena los nodos de un Quad4 en sentido antihorario
    usando el centroide como referencia angular.
    """
    coords = np.array([n.coordenadas for n in nodos])
    cx, cy = np.mean(coords, axis=0)
    angulos = np.arctan2(coords[:, 1] - cy, coords[:, 0] - cx)
    orden = np.argsort(angulos)
    return [nodos[i] for i in orden]

