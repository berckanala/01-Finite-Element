import meshio
import numpy as np
from collections import defaultdict
from openseespy import opensees as ops
import pyvista as pv

# === Leer malla base ===
mesh = meshio.read("prueba_god.msh")
points = mesh.points

field_data = mesh.field_data
tag_id_to_name = {v[0]: k for k, v in field_data.items()}

volumes = defaultdict(list)
for cell_block, tag_block in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    if cell_block.type == "tetra":
        for conn, tag in zip(cell_block.data, tag_block):
            name = tag_id_to_name.get(tag, f"Tag_{tag}")
            volumes[name].append(conn)

volume_nodes = {}
for name, tets in volumes.items():
    nodes_in_volume = np.unique(np.array(tets).flatten())
    volume_nodes[name] = nodes_in_volume

# === Identificar nodos de carga
selected_nodes=[]
for i,p in enumerate(points):
    if i in volume_nodes["Viga"]:
        if abs(p[1]-30)<1e-3 and abs(p[0]-110)<1e-3:
            selected_nodes.append(i)

if not selected_nodes:
    raise ValueError("❌ No se encontraron nodos de carga.")

# === Función que reconstruye modelo y calcula tensiones
def calcular_von_mises(points, volumes, volume_nodes, selected_nodes, exclude_elements=None):
    if exclude_elements is None:
        exclude_elements = []

    ops.wipe()
    ops.model("basic", "-ndm",3,"-ndf",3)

    for i,coord in enumerate(points):
        ops.node(i+1,*coord)

    ops.nDMaterial("ElasticIsotropic",1,3.5e9,0.35,1250)
    ops.nDMaterial("ElasticIsotropic",2,1e20,0,7850)

    element_id=1
    idx_to_eid={}
    viga_idx_to_eid={}
    for idx,conn in enumerate(volumes["Viga"]):
        if idx in exclude_elements:
            continue
        node_tags = [int(n+1) for n in conn]
        ops.element("FourNodeTetrahedron",element_id,*node_tags,1)
        idx_to_eid[element_id]=idx
        viga_idx_to_eid[idx]=element_id
        element_id+=1

    for conn in volumes["BC_R1"]:
        node_tags=[int(n+1) for n in conn]
        ops.element("FourNodeTetrahedron",element_id,*node_tags,2)
        element_id+=1
    for conn in volumes["BC_1"]:
        node_tags=[int(n+1) for n in conn]
        ops.element("FourNodeTetrahedron",element_id,*node_tags,2)
        element_id+=1

    for n in volume_nodes["BC_R1"]:
        ops.fix(int(n+1),1,1,1)
    for n in volume_nodes["BC_1"]:
        ops.fix(int(n+1),1,0,1)

    ops.timeSeries("Linear",1)
    ops.pattern("Plain",1,1)
    for n in selected_nodes:
        ops.load(int(n+1),0,-1e10,0)

    ops.system("ProfileSPD")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.algorithm("Linear")
    ops.integrator("LoadControl",1.0)
    ops.analysis("Static")
    print("✅ Ejecutando análisis...")
    ok=ops.analyze(1)
    if ok!=0:
        return None,None,None

    resultado=[]
    for eid,idx in idx_to_eid.items():
        stress=ops.eleResponse(eid,"stresses")
        if not stress or len(stress)<6:
            vm=0.0
        else:
            sx,sy,sz,tyz,txz,txy=stress
            vm = np.sqrt(0.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2)+3*(txy**2+tyz**2+txz**2))
        resultado.append((idx,vm))

    return resultado, idx_to_eid, viga_idx_to_eid

# === Seleccionar qué eliminar
def seleccionar_a_eliminar(tensiones, percent_to_remove=0.10, protected_indices=None):
    if protected_indices is None:
        protected_indices=[]
    tensiones_ordenadas = sorted(tensiones,key=lambda x:x[1])
    n_eliminar = int(len(tensiones_ordenadas)*percent_to_remove)
    eliminados=[]
    for idx,vm in tensiones_ordenadas:
        if idx in protected_indices:
            continue
        eliminados.append(idx)
        if len(eliminados)>=n_eliminar:
            break
    return eliminados

# === Iteraciones
n_iteraciones=50
percent_to_remove=0.0009
exclude_indices=[]

for iteracion in range(n_iteraciones):
    print(f"\n=== Iteración {iteracion+1}/{n_iteraciones} ===")

    # Calcular tensiones
    tensiones, idx_to_eid, _ = calcular_von_mises(points, volumes, volume_nodes, selected_nodes, exclude_indices)

    if tensiones is None:
        print("❌ No se pudo resolver con los elementos actuales. Finalizando.")
        break

    vm_vals=[vm for _,vm in tensiones if vm>0]
    vmax = max(vm_vals)
    vmin = min(vm_vals)
    print(f"Von Mises min: {vmin:.2e}, max: {vmax:.2e}")

    nodos_criticos = set(volume_nodes["BC_R1"]).union(volume_nodes["BC_1"]).union(selected_nodes)
    protegidos=[]
    for idx,conn in enumerate(volumes["Viga"]):
        if any(n in nodos_criticos for n in conn):
            protegidos.append(idx)

    candidatos=seleccionar_a_eliminar(tensiones, percent_to_remove, protegidos)
    print(f"Intentando eliminar {len(candidatos)} candidatos...")

    confirmados=[]
    for idx in candidatos:
        prueba_exclude = exclude_indices + [idx]
        resultado, _, _ = calcular_von_mises(points, volumes, volume_nodes, selected_nodes, prueba_exclude)
        if resultado is not None:
            print(f"✅ Eliminado elemento {idx}")
            confirmados.append(idx)
        else:
            print(f"❌ No se puede eliminar elemento {idx} (matriz singular)")

    if not confirmados:
        print("⚠️ Ningún elemento se pudo eliminar. Finalizando.")
        break

    exclude_indices.extend(confirmados)

    activos = [idx for idx,_ in tensiones if idx not in exclude_indices]
    cells = np.hstack([
        np.full((len(activos),1),4),
        np.array([volumes["Viga"][i] for i in activos])
    ])
    vm_array = np.array([vm for idx,vm in tensiones if idx in activos])

    grid = pv.UnstructuredGrid(cells, np.full(len(activos),pv.CellType.TETRA), points)
    grid.cell_data["von_mises"]=vm_array

    p = pv.Plotter()
    p.add_mesh(grid, scalars="von_mises", cmap="viridis", show_edges=True, opacity=0.8)
    p.add_title(f"Iteración {iteracion+1}")
    if iteracion ==10 or iteracion == 20 or iteracion == 30 or iteracion == 40:
        p.show()
