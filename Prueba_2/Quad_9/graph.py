import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import os
import re

def clean_filename(text):
    """Convierte un string en nombre de archivo válido eliminando LaTeX y símbolos especiales."""
    return re.sub(r'[^\w\-_. ]', '', text.replace('$', '').replace('\\', '').replace('{', '').replace('}', ''))

import matplotlib.pyplot as plt
import numpy as np
import os

def plot_all_elements(elements, title, show_ids=True):
    """
    Dibuja todos los elementos (compatible con Quad2D y CST).
    Guarda la imagen como PNG en la carpeta GRAFICOS/.
    """
    all_x = []
    all_y = []

    for elem in elements:
        coords = np.array([node.coord for node in elem.node_list])
        coords = np.vstack([coords, coords[0]])  # cerrar polígono
        all_x.extend(coords[:, 0])
        all_y.extend(coords[:, 1])

    # Márgenes y tamaño gráfico
    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05
    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin
    fixed_width = 8
    aspect_ratio = y_range / x_range
    height = fixed_width * aspect_ratio

    fig, ax = plt.subplots(figsize=(fixed_width, height))

    for elem in elements:
        coords = np.array([node.coord for node in elem.node_list])
        coords = np.vstack([coords, coords[0]])  # cerrar polígono
        ax.plot(coords[:, 0], coords[:, 1], 'k-', linewidth=1)

        if show_ids:
            for node in elem.node_list:
                x, y = node.coord
                if np.all(node.restrain == [1, 1]):
                    pass  # Puedes activar esto si quieres mostrar nombres
                    # ax.text(x, y, f'N{node.index}', color='blue', fontsize=6, ha='center', va='center')

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("All Quad2D elements")
    ax.grid(True)

    os.makedirs("GRAFICOS", exist_ok=True)
    plt.savefig(f"GRAFICOS/{title}_elementos.png", dpi=300, bbox_inches='tight')
    plt.close()


import matplotlib.pyplot as plt
import numpy as np
import os

def plot_applied_forces(nodes, elements, title, f_vector, scale=1e-2):
    all_x = []
    all_y = []

    # Recolectar coordenadas de todos los elementos
    for elem in elements:
        coords = np.array([node.coord for node in elem.node_list])
        all_x.extend(coords[:, 0])
        all_y.extend(coords[:, 1])

    # También considerar nodos, por si hay fuerzas fuera del dominio de elementos
    for node in nodes:
        all_x.append(node.coord[0])
        all_y.append(node.coord[1])

    # Calcular límites con márgenes
    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05
    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin

    # Escalar altura según ancho fijo
    fixed_width = 8  # pulgadas
    aspect_ratio = y_range / x_range
    height = fixed_width * aspect_ratio

    fig, ax = plt.subplots(figsize=(fixed_width, height))
    ax.set_title("Fuerzas aplicadas en nodos")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True)

    # Dibujar contorno de los elementos
    for elem in elements:
        coords = np.array([node.coord for node in elem.node_list])
        coords = np.vstack([coords, coords[0]])  # cerrar polígono
        ax.plot(coords[:, 0], coords[:, 1], color='lightgray', linewidth=1, zorder=1)

    # Dibujar flechas de fuerzas
    for node in nodes:
        x, y = node.coord
        dof_x, dof_y = node.dofs
        fx = f_vector[dof_x][0] if dof_x < len(f_vector) else 0.0
        fy = f_vector[dof_y][0] if dof_y < len(f_vector) else 0.0

        if not np.isclose(fx, 0.0) or not np.isclose(fy, 0.0):
            ax.quiver(x, y, fx, fy,
                      angles='xy', scale_units='xy',
                      scale=1 / scale, color='red', width=0.005,
                      zorder=5)

    os.makedirs("GRAFICOS", exist_ok=True)
    plt.savefig(f"GRAFICOS/{title}_fuerzas.png", dpi=300, bbox_inches='tight')
    plt.close()

import matplotlib.tri as mtri

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import os

def plot_deformed_structure(elements, title, scale=1.0, show_ids=False):
    points = []
    displacements = []
    triangles = []
    point_index = {}
    idx_counter = 0

    for elem in elements:
        coords = np.array([node.coord for node in elem.node_list])  # (9,2)
        u = np.zeros_like(coords)

        # Obtener desplazamientos globales del nodo
        for i, nodo in enumerate(elem.node_list):
            ux = nodo.structure.u_global[nodo.dofs[0], 0]
            uy = nodo.structure.u_global[nodo.dofs[1], 0]
            u[i] = [ux, uy]

        def_coords = coords + scale * u

        tri = []
        for i in range(len(def_coords)):
            key = tuple(def_coords[i])
            if key not in point_index:
                point_index[key] = idx_counter
                points.append(key)
                mag = np.linalg.norm(u[i])
                displacements.append(mag)
                idx_counter += 1
            tri.append(point_index[key])

        # Triangulación para Quad9:
        # Divide el cuadrilátero en 4 triángulos usando:
        # nodos de esquina (0,1,2,3) y nodo central (8)
        if len(tri) == 9:
            triangles.append([tri[0], tri[1], tri[8]])
            triangles.append([tri[1], tri[2], tri[8]])
            triangles.append([tri[2], tri[3], tri[8]])
            triangles.append([tri[3], tri[0], tri[8]])
        elif len(tri) == 4:
            # Por si hay elementos Quad4
            triangles.append([tri[0], tri[1], tri[2]])
            triangles.append([tri[0], tri[2], tri[3]])
        elif len(tri) == 3:
            triangles.append(tri)
        else:
            print(f"⚠️ Elemento con {len(tri)} nodos no soportado.")

    points = np.array(points)
    displacements = np.array(displacements)
    triangles = np.array(triangles)

    if len(triangles) == 0:
        print("❌ No se generaron triángulos para graficar deformada.")
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles)

    tpc = ax.tripcolor(triang, displacements, shading='gouraud', cmap='viridis')
    cb = plt.colorbar(tpc, ax=ax)
    cb.set_label('Magnitud de desplazamiento [m]')

    ax.set_aspect('equal')
    ax.set_title(f"Estructura deformada (×{scale}) - Mapa de calor")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    os.makedirs("GRAFICOS", exist_ok=True)
    plt.savefig(f"GRAFICOS/{title}_deformada_colormap.png", dpi=300, bbox_inches='tight')
    plt.close()

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import os

def plot_deformed_with_reactions(title, elements, reactions, scale=1.0, reaction_scale=1e-3, show_ids=False):
    all_x = []
    all_y = []

    for elem in elements:
        coords = np.array([n.coord for n in elem.node_list])
        all_x.extend(coords[:, 0])
        all_y.extend(coords[:, 1])

    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05
    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin

    fixed_width = 16
    aspect_ratio = y_range / x_range
    height = fixed_width * aspect_ratio / 2

    fig, axs = plt.subplots(1, 2, figsize=(fixed_width, height))
    ax1, ax2 = axs

    # === GRÁFICO 1: Reacciones nodales ===
    for elem in elements:
        coords = np.array([n.coord for n in elem.node_list])
        coords_closed = np.vstack([coords, coords[0]])
        ax1.plot(coords_closed[:, 0], coords_closed[:, 1], 'lightgray', linewidth=0.5)

        for nodo in elem.node_list:
            x, y = nodo.coord
            rx = reactions[nodo.dofs[0]][0] if nodo.restrain[0] == 1 else 0.0
            ry = reactions[nodo.dofs[1]][0] if nodo.restrain[1] == 1 else 0.0

            if rx != 0.0 or ry != 0.0:
                ax1.quiver(x, y, rx, ry, angles='xy', scale_units='xy',
                           scale=1/reaction_scale, color='red', width=0.005)
            if show_ids:
                ax1.text(x, y, f'N{nodo.index}', fontsize=6)

    ax1.set_title("Nodal reactions", fontsize=14)
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.grid(True)
    ax1.set_xlim(x_min - x_margin, x_max + x_margin)
    ax1.set_ylim(y_min - y_margin, y_max + y_margin)
    ax1.set_aspect('equal', adjustable='box')

    # === GRÁFICO 2: Heatmap de reacciones ===
    nodal_values = {}
    for elem in elements:
        for nodo in elem.node_list:
            rx = reactions[nodo.dofs[0]][0] if nodo.restrain[0] == 1 else 0.0
            ry = reactions[nodo.dofs[1]][0] if nodo.restrain[1] == 1 else 0.0
            mag = np.sqrt(rx**2 + ry**2)
            nodal_values[nodo.index] = mag

    all_nodes = {n.index: n for e in elements for n in e.node_list}
    node_list = list(all_nodes.values())
    node_id_map = {node.index: i for i, node in enumerate(node_list)}
    points = np.array([node.coord for node in node_list])
    triangles = []
    values = []

    for elem in elements:
        ids = [n.index for n in elem.node_list]
        val = np.mean([nodal_values.get(nid, 0.0) for nid in ids])
        if val > 0:
            tri = [node_id_map[nid] for nid in ids]
            if len(tri) == 9:
                # Dividir Quad9 en 4 triángulos con el centro (nodo 8)
                triangles.append([tri[0], tri[1], tri[8]])
                triangles.append([tri[1], tri[2], tri[8]])
                triangles.append([tri[2], tri[3], tri[8]])
                triangles.append([tri[3], tri[0], tri[8]])
                values.extend([val, val, val, val])
            elif len(tri) == 4:
                # Dividir Quad4 en 2 triángulos
                triangles.append([tri[0], tri[1], tri[2]])
                triangles.append([tri[0], tri[2], tri[3]])
                values.extend([val, val])
            elif len(tri) == 3:
                triangles.append(tri)
                values.append(val)
            else:
                print(f"⚠️ Elemento con {len(tri)} nodos no soportado en heatmap.")

    triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles)

    for elem in elements:
        coords = np.array([n.coord for n in elem.node_list])
        coords_closed = np.vstack([coords, coords[0]])
        ax2.plot(coords_closed[:, 0], coords_closed[:, 1], color='lightgray', linewidth=0.5)

    if values:
        tpc = ax2.tripcolor(triang, facecolors=values, edgecolors='k', cmap='plasma', zorder=2)
        cbar = fig.colorbar(tpc, ax=ax2)
        cbar.set_label("Reaction force [N]", fontsize=14)

    ax2.set_title("Heatmap of reactions per element", fontsize=14)
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.grid(True)
    ax2.set_xlim(x_min - x_margin, x_max + x_margin)
    ax2.set_ylim(y_min - y_margin, y_max + y_margin)
    ax2.set_aspect('equal', adjustable='box')

    os.makedirs("GRAFICOS", exist_ok=True)
    fig.savefig(f"GRAFICOS/{title}_deformada_reacciones.png", dpi=300, bbox_inches='tight')
    plt.close()

    # === Sumatoria de reacciones ===
    total_rx = 0.0
    total_ry = 0.0
    for node in all_nodes.values():
        rx = reactions[node.dofs[0], 0] if node.restrain[0] == 1 else 0.0
        ry = reactions[node.dofs[1], 0] if node.restrain[1] == 1 else 0.0
        total_rx += rx
        total_ry += ry

    total_magnitude = np.sqrt(total_rx**2 + total_ry**2)
    print(f"📌 Sumatoria total de reacciones: Rx = {total_rx:.3f} N, Ry = {total_ry:.3f} N, ||R|| = {total_magnitude:.3f} N")


def plot_von_mises_field(nodes, elements, vm_nodal_dict, title, cmap='plasma'):
    node_id_to_index = {}
    xs, ys, vms = [], [], []

    # Asociar cada nodo con un índice para el triangulado
    for i, node in enumerate(nodes):
        node_id_to_index[node.index] = i
        xs.append(node.coord[0])
        ys.append(node.coord[1])

        vms.append(vm_nodal_dict.get(node.index, 0.0))  # Usar node.index aquí

    # Construir triangulación a partir de nodos de elementos
    triangles = []
    for elem in elements:
        n0 = node_id_to_index[elem.node_list[0].index]
        n1 = node_id_to_index[elem.node_list[1].index]
        n2 = node_id_to_index[elem.node_list[2].index]
        n3 = node_id_to_index[elem.node_list[3].index]

        # Dos triángulos: (n0, n1, n2) y (n0, n2, n3)
        triangles.append([n0, n1, n2])
        triangles.append([n0, n2, n3])


    triang = mtri.Triangulation(xs, ys, triangles)

    # Márgenes y proporciones
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin
    fixed_width = 8
    height = fixed_width * (y_range / x_range)

    fig, ax = plt.subplots(figsize=(fixed_width, height))
    tcf = ax.tricontourf(triang, vms, levels=20, cmap=cmap)

    cbar = fig.colorbar(tcf, ax=ax)
    cbar.set_label("Von Mises Stress (MPa)")

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title("Von Mises stress field over elements")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    # Guardar gráfico
    os.makedirs("GRAFICOS", exist_ok=True)
    fig.savefig(f"GRAFICOS/{title}_von_mises.png", dpi=300, bbox_inches='tight')
    plt.close()

from matplotlib.colors import BoundaryNorm

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm

from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import os
from matplotlib.colors import LogNorm

def plot_von_mises_per_element(nodes, elements, vm_nodal_dict, title, cmap='plasma'):
    node_id_to_index = {}
    xs, ys = [], []

    # Indexar nodos y almacenar coordenadas
    for i, node in enumerate(nodes):
        node_id_to_index[node.index] = i
        xs.append(node.coord[0])
        ys.append(node.coord[1])

    triangles = []
    element_colors = []

    for elem in elements:
        n_ids = [node.index for node in elem.node_list]
        try:
            # Para Quad9 esperamos 9 nodos
            if len(n_ids) == 9:
                idx = [node_id_to_index[nid] for nid in n_ids]
                # Dividir en 4 triángulos usando nodo central idx[8]
                triangles.append([idx[0], idx[1], idx[8]])
                triangles.append([idx[1], idx[2], idx[8]])
                triangles.append([idx[2], idx[3], idx[8]])
                triangles.append([idx[3], idx[0], idx[8]])
            elif len(n_ids) == 4:
                idx = [node_id_to_index[nid] for nid in n_ids]
                triangles.append([idx[0], idx[1], idx[2]])
                triangles.append([idx[0], idx[2], idx[3]])
            else:
                print(f"⚠️ Elemento con {len(n_ids)} nodos no soportado para Von Mises.")
                continue

            # Promedio Von Mises en nodos del elemento
            vms_elem = np.mean([vm_nodal_dict.get(nid, 0.0) for nid in n_ids])

            # Asegurar valor positivo para escala logarítmica
            if vms_elem <= 0:
                vms_elem = 1e-6

            # Repetir valor para cada triángulo generado
            if len(n_ids) == 9:
                element_colors.extend([vms_elem]*4)
            elif len(n_ids) == 4:
                element_colors.extend([vms_elem]*2)

        except KeyError:
            continue  # Ignora elementos con nodos faltantes

    triang = mtri.Triangulation(xs, ys, triangles)

    # Escalado gráfico
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin
    aspect_ratio = y_range / x_range
    fixed_width = 8
    height = fixed_width * aspect_ratio

    fig, ax = plt.subplots(figsize=(fixed_width, height))

    vmin = max(min(element_colors), 1e-6)
    vmax = max(element_colors)
    norm = LogNorm(vmin=vmin, vmax=vmax)

    tpc = ax.tripcolor(
        triang,
        facecolors=element_colors,
        cmap=cmap,
        norm=norm,
        edgecolors='k',
        shading='flat'
    )

    cbar = fig.colorbar(tpc, ax=ax)
    cbar.set_label("Von Mises Stress (MPa)")

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal')
    ax.set_title("Von Mises stress per element")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    os.makedirs("GRAFICOS", exist_ok=True)
    fig.savefig(f"GRAFICOS/{title}_von_mises_per_element.png", dpi=300, bbox_inches='tight')
    plt.close()


def compute_nodal_principal_fields(nodes, elements, u_global):
    """
    Calcula los esfuerzos y deformaciones principales σ1, σ2, ε1, ε2 por nodo.
    Devuelve 4 diccionarios: sigma1, sigma2, eps1, eps2
    """
    sigma1_map = {node.index: [] for node in nodes}
    sigma2_map = {node.index: [] for node in nodes}
    eps1_map = {node.index: [] for node in nodes}
    eps2_map = {node.index: [] for node in nodes}

    for elem in elements:
        σ = elem.get_stress(u_global)  # np.array([σx, σy, τxy])
        ε = elem.get_strain(u_global)  # np.array([εx, εy, γxy])

        σx, σy, τxy = σ
        εx, εy, γxy = ε

        # Esfuerzos principales
        σ_avg = 0.5 * (σx + σy)
        Rσ = np.sqrt(((σx - σy) / 2)**2 + τxy**2)
        σ1, σ2 = σ_avg + Rσ, σ_avg - Rσ

        # Deformaciones principales
        ε_avg = 0.5 * (εx + εy)
        Rε = np.sqrt(((εx - εy) / 2)**2 + (γxy / 2)**2)
        ε1, ε2 = ε_avg + Rε, ε_avg - Rε

        for node in elem.node_list:
            sigma1_map[node.index].append(σ1)
            sigma2_map[node.index].append(σ2)
            eps1_map[node.index].append(ε1)
            eps2_map[node.index].append(ε2)

    # Promediar por nodo
    sigma1 = {nid: np.mean(vals) for nid, vals in sigma1_map.items()}
    sigma2 = {nid: np.mean(vals) for nid, vals in sigma2_map.items()}
    eps1 = {nid: np.mean(vals) for nid, vals in eps1_map.items()}
    eps2 = {nid: np.mean(vals) for nid, vals in eps2_map.items()}

    return sigma1, sigma2, eps1, eps2

def plot_all_scalar_fields_separately(nodes, elements, nodal_fields, title_prefix):
    """
    Genera 6 gráficos individuales para: σxx, σyy, τxy, εxx, εyy, γxy
    """
    σxx, σyy, τxy, εxx, εyy, γxy = {}, {}, {}, {}, {}, {}

    for node_id, (sigma, epsilon) in nodal_fields.items():
        σxx[node_id] = sigma[0]
        σyy[node_id] = sigma[1]
        τxy[node_id] = sigma[2]
        εxx[node_id] = epsilon[0]
        εyy[node_id] = epsilon[1]
        γxy[node_id] = epsilon[2]

    plot_scalar_field(nodes, elements, σxx, r"$\sigma_{xx}$ (MPa)", f"{title_prefix} - sigma_xx")
    plot_scalar_field(nodes, elements, σyy, r"$\sigma_{yy}$ (MPa)", f"{title_prefix} - sigma_yy")
    plot_scalar_field(nodes, elements, τxy, r"$\tau_{xy}$ (MPa)", f"{title_prefix} - tau_xy")
    plot_scalar_field(nodes, elements, εxx, r"$\varepsilon_{xx}$", f"{title_prefix} - epsilon_xx")
    plot_scalar_field(nodes, elements, εyy, r"$\varepsilon_{yy}$", f"{title_prefix} - epsilon_yy")
    plot_scalar_field(nodes, elements, γxy, r"$\gamma_{xy}$", f"{title_prefix} - gamma_xy")

    plot_scalar_field_per_element(nodes, elements, σxx, r"$\sigma_{xx}$ (MPa)", f"{title_prefix} - sigma_xx_per_element")
    plot_scalar_field_per_element(nodes, elements, σyy, r"$\sigma_{yy}$ (MPa)", f"{title_prefix} - sigma_yy_per_element")
    plot_scalar_field_per_element(nodes, elements, τxy, r"$\tau_{xy}$ (MPa)", f"{title_prefix} - tau_xy_per_element")
    plot_scalar_field_per_element(nodes, elements, εxx, r"$\varepsilon_{xx}$", f"{title_prefix} - epsilon_xx_per_element")
    plot_scalar_field_per_element(nodes, elements, εyy, r"$\varepsilon_{yy}$", f"{title_prefix} - epsilon_yy_per_element")
    plot_scalar_field_per_element(nodes, elements, γxy, r"$\gamma_{xy}$", f"{title_prefix} - gamma_xy_per_element")

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import os

def plot_scalar_field(nodes, elements, nodal_values, field_title, filename_prefix, cmap='plasma'):
    """
    Grafica un campo escalar interpolado sobre una malla de elementos Quad9.
    """
    node_id_to_index = {node.index: i for i, node in enumerate(nodes)}
    xs = [node.coord[0] for node in nodes]
    ys = [node.coord[1] for node in nodes]

    triangles = []
    for elem in elements:
        ids = [node.index for node in elem.node_list]
        try:
            if len(ids) == 9:
                i0, i1, i2, i3, i4, i5, i6, i7, i8 = [node_id_to_index[nid] for nid in ids]
                # Dividir quad9 en 4 triángulos usando nodo central (i8)
                triangles.append([i0, i1, i8])
                triangles.append([i1, i2, i8])
                triangles.append([i2, i3, i8])
                triangles.append([i3, i0, i8])
            elif len(ids) == 4:
                i0, i1, i2, i3 = [node_id_to_index[nid] for nid in ids]
                triangles.append([i0, i1, i2])
                triangles.append([i0, i2, i3])
            else:
                print(f"⚠️ Elemento con {len(ids)} nodos no soportado en plot_scalar_field.")
                continue
        except KeyError:
            continue

    values = np.zeros(len(nodes))
    for node in nodes:
        values[node_id_to_index[node.index]] = nodal_values.get(node.index, 0.0)

    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_margin = 0.05 * (x_max - x_min)
    y_margin = 0.05 * (y_max - y_min)

    fig_width = 8
    aspect_ratio = (y_max - y_min + 2 * y_margin) / (x_max - x_min + 2 * x_margin)
    fig_height = fig_width * aspect_ratio

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    triang = mtri.Triangulation(xs, ys, triangles)

    tcf = ax.tricontourf(triang, values, levels=20, cmap=cmap)
    cbar = fig.colorbar(tcf, ax=ax)
    cbar.set_label(field_title)

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal')
    ax.set_title(field_title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    os.makedirs("GRAFICOS", exist_ok=True)
    fig.savefig(f"GRAFICOS/{filename_prefix}.png", dpi=300, bbox_inches='tight')
    plt.close()

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import Normalize
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import os
from matplotlib.colors import Normalize

def plot_scalar_field_per_element(nodes, elements, nodal_values, field_title, filename_prefix, cmap='plasma'):
    node_id_to_index = {node.index: i for i, node in enumerate(nodes)}
    xs = [node.coord[0] for node in nodes]
    ys = [node.coord[1] for node in nodes]

    triangles = []
    element_values = []

    for elem in elements:
        try:
            ids = [node.index for node in elem.node_list]
            if len(ids) == 9:
                i0, i1, i2, i3, i4, i5, i6, i7, i8 = [node_id_to_index[nid] for nid in ids]
                # Dividir quad9 en 4 triángulos con nodo central i8
                triangles.append([i0, i1, i8])
                triangles.append([i1, i2, i8])
                triangles.append([i2, i3, i8])
                triangles.append([i3, i0, i8])
            elif len(ids) == 4:
                i0, i1, i2, i3 = [node_id_to_index[nid] for nid in ids]
                triangles.append([i0, i1, i2])
                triangles.append([i0, i2, i3])
            else:
                print(f"⚠️ Elemento con {len(ids)} nodos no soportado.")
                continue

            vals = [nodal_values.get(nid, 0.0) for nid in ids]
            value = np.mean(vals)

            # Asignar valor promedio a cada triángulo generado
            if len(ids) == 9:
                element_values.extend([value]*4)
            else:
                element_values.extend([value]*2)

        except KeyError:
            continue

    triangles = np.array(triangles)
    element_values = np.array(element_values)

    triang = mtri.Triangulation(xs, ys, triangles)

    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05
    height = 8 * ((y_max - y_min + 2 * y_margin) / (x_max - x_min + 2 * x_margin))

    fig, ax = plt.subplots(figsize=(8, height))

    vmin = np.min(element_values)
    vmax = np.max(element_values)
    norm = Normalize(vmin=vmin, vmax=vmax)

    tpc = ax.tripcolor(triang, facecolors=element_values, cmap=cmap,
                       norm=norm, edgecolors='k', shading='flat')

    cbar = fig.colorbar(tpc, ax=ax)
    cbar.set_label(field_title)

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal')
    ax.set_title(field_title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    os.makedirs("GRAFICOS", exist_ok=True)
    fig.savefig(f"GRAFICOS/{filename_prefix}.png", dpi=300, bbox_inches='tight')
    plt.close()

from matplotlib.colors import LogNorm
import matplotlib.tri as mtri

from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import os

def plot_safety_field_per_element_log(nodes, elements, nodal_values, field_title, filename_prefix, cmap='viridis_r'):
    """
    Grafica un campo de factores de seguridad (FS) por elemento usando escala logarítmica.
    Admite elementos cuadriláteros dividiéndolos en dos triángulos.
    """
    node_id_to_index = {node.id: i for i, node in enumerate(nodes)}
    xs = [node.coord[0] for node in nodes]
    ys = [node.coord[1] for node in nodes]

    triangles = []
    element_values = []

    for elem in elements:
        try:
            n = elem.node_list
            # Dividir el cuadrado en dos triángulos: [n0, n1, n2] y [n0, n2, n3]
            tri1 = [node_id_to_index[n[0].id], node_id_to_index[n[1].id], node_id_to_index[n[2].id]]
            tri2 = [node_id_to_index[n[0].id], node_id_to_index[n[2].id], node_id_to_index[n[3].id]]

            triangles.append(tri1)
            triangles.append(tri2)

            # Usar el promedio de valores nodales para ambos triángulos
            value = np.mean([nodal_values[n[i].id] for i in range(4)])
            element_values.extend([value, value])
        except KeyError:
            continue

    triangles = np.array(triangles)
    element_values = np.array(element_values)

    triang = mtri.Triangulation(xs, ys, triangles)

    # Definir límites del gráfico
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05
    aspect_ratio = (y_max - y_min + 2 * y_margin) / (x_max - x_min + 2 * x_margin)
    height = 8 * aspect_ratio

    fig, ax = plt.subplots(figsize=(8, height))

    vmin = max(np.min(element_values), 1e-3)
    vmax = np.max(element_values)
    norm = LogNorm(vmin=vmin, vmax=vmax)

    tpc = ax.tripcolor(
        triang,
        facecolors=element_values,
        cmap=cmap,
        norm=norm,
        edgecolors='k',
        shading='flat'
    )

    cbar = fig.colorbar(tpc, ax=ax)
    cbar.set_label(field_title)

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal')
    ax.set_title(field_title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    os.makedirs("GRAFICOS", exist_ok=True)
    fig.savefig(f"GRAFICOS/{filename_prefix}.png", dpi=300, bbox_inches='tight')
    plt.close()



def plot_safety_factors(nodes, elements, u_global, title_prefix, sigma_y_tension, sigma_y_compression):
    """
    Genera dos mapas de factor de seguridad:
    - Tracción: FS = sigma_y_tension / abs(sigma_1)
    - Compresión: FS = sigma_y_compression / abs(sigma_2)
    También imprime y muestra el mínimo FS en cada caso.
    """
    # Obtener tensiones principales por nodo
    sigma1, sigma2, _, _ = compute_nodal_principal_fields(nodes, elements, u_global)

    # Calcular FS por nodo
    fs_tension = {nid: sigma_y_tension / max(abs(val), 1e-6) for nid, val in sigma1.items()}
    fs_compression = {nid: sigma_y_compression / max(abs(val), 1e-6) for nid, val in sigma2.items()}

    # Calcular mínimos
    min_fs_tension = min(fs_tension.values())
    min_fs_compression = min(fs_compression.values())

    print(f"🔧 Mínimo FS Tracción: {min_fs_tension:.3f}")
    print(f"🔧 Mínimo FS Compresión: {min_fs_compression:.3f}")

    # Graficar por elemento
    plot_safety_field_per_element_log(
        nodes, elements, fs_tension,
        field_title=f"FS Tracción (min = {min_fs_tension:.2f})",
        filename_prefix=f"{title_prefix}_FS_Tension_Log"
    )

    plot_safety_field_per_element_log(
        nodes, elements, fs_compression,
        field_title=f"FS Compresión (min = {min_fs_compression:.2f})",
        filename_prefix=f"{title_prefix}_FS_Compression_Log"
    )


def plot_principal_fields(nodes, elements, u_global, title_prefix, sigma_y_tension, sigma_y_compression):
        """
        Genera 4 gráficos: σ1, σ2, ε1, ε2.
        """
        sigma1, sigma2, eps1, eps2 = compute_nodal_principal_fields(nodes, elements, u_global)

        plot_scalar_field(nodes, elements, sigma1, r"$\sigma_1$ (MPa)", f"{title_prefix} - sigma_1")
        plot_scalar_field(nodes, elements, sigma2, r"$\sigma_2$ (MPa)", f"{title_prefix} - sigma_2")
        plot_scalar_field(nodes, elements, eps1, r"$\varepsilon_1$", f"{title_prefix} - epsilon_1")
        plot_scalar_field(nodes, elements, eps2, r"$\varepsilon_2$", f"{title_prefix} - epsilon_2")

        plot_scalar_field_per_element(nodes, elements, sigma1, r"$\sigma_1$ (MPa)", f"{title_prefix} - sigma_1_per_element")
        plot_scalar_field_per_element(nodes, elements, sigma2, r"$\sigma_2$ (MPa)", f"{title_prefix} - sigma_2_per_element")
        plot_scalar_field_per_element(nodes, elements, eps1, r"$\varepsilon_1$", f"{title_prefix} - epsilon_1_per_element")
        plot_scalar_field_per_element(nodes, elements, eps2, r"$\varepsilon_2$", f"{title_prefix} - epsilon_2_per_element")
        # ⚠️ Reparar `id` temporalmente antes de graficar FS
        for node in nodes:
            node.id = node.index
        
        plot_safety_factors(nodes, elements, u_global,title_prefix, sigma_y_tension, sigma_y_compression)

import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import numpy as np

def plot_elements_by_thickness(elements, title, cmap='viridis'):
    """
    Dibuja y guarda una figura coloreando los elementos según el espesor de su sección,
    ajustando automáticamente el alto para mantener proporciones.
    """
    node_ids = {}
    xs, ys = [], []
    triangles = []
    thicknesses = []

    counter = 0
    for elem in elements:
        triangle = []
        for node in elem.node_list:
            if node.id not in node_ids:
                node_ids[node.id] = counter
                xs.append(node.coord[0])
                ys.append(node.coord[1])
                counter += 1
            triangle.append(node_ids[node.id])
        triangles.append(triangle)
        thicknesses.append(elem.section.thickness)

    triang = mtri.Triangulation(xs, ys, triangles)

    # Cálculo de límites y dimensiones
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    x_range = (x_max - x_min) + 2 * x_margin
    y_range = (y_max - y_min) + 2 * y_margin

    fixed_width = 8
    aspect_ratio = y_range / x_range
    height = fixed_width * aspect_ratio

    fig, ax = plt.subplots(figsize=(fixed_width, height))
    tpc = ax.tripcolor(triang, facecolors=thicknesses, edgecolors='k', cmap=cmap)
    cbar = plt.colorbar(tpc, ax=ax, label="Thickness (mm)")
    ax.set_title("Thickness Distribution per Element")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True)

    plt.tight_layout()
    plt.savefig(f"GRAFICOS/{title}_espesores_post_topo.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Gráfico guardado como GRAFICOS/{title}_espesores.png")

def compute_nodal_von_mises(elements, u_global):
    """
    Promedia los esfuerzos de Von Mises en los nodos a partir de los elementos vecinos.

    Args:
        elements (list): Lista de elementos (Quad2D, CST, etc.)
        u_global (ndarray): Vector de desplazamientos global.

    Returns:
        dict: Diccionario {node.index: von Mises promedio}
    """
    nodal_vm = {}  # node.index : [vm1, vm2, ...]

    for elem in elements:
        vm = elem.von_mises_stress(u_global)
        for node in elem.node_list:
            if node.index not in nodal_vm:
                nodal_vm[node.index] = []
            nodal_vm[node.index].append(vm)

    # Promediar los esfuerzos por nodo
    nodal_vm_avg = {node_index: np.mean(vms) for node_index, vms in nodal_vm.items()}
    return nodal_vm_avg

def compute_nodal_stress_strain(nodes, elements, u_global):
    """
    Calcula σxx, σyy, σxy, εxx, εyy, εxy por nodo (promedio de elementos conectados).
    Retorna: diccionario {node.index: (σ_vec, ε_vec)}
    """
    stress_map = {node.index: [] for node in nodes}
    strain_map = {node.index: [] for node in nodes}

    for elem in elements:
        stress = elem.get_stress(u_global)  # Devuelve np.array([σx, σy, τxy])
        strain = elem.get_strain(u_global)  # Devuelve np.array([εx, εy, γxy])
        for node in elem.node_list:
            stress_map[node.index].append(stress)
            strain_map[node.index].append(strain)

    result = {}
    for node in nodes:
        σ_avg = np.mean(stress_map[node.index], axis=0) if stress_map[node.index] else np.zeros(3)
        ε_avg = np.mean(strain_map[node.index], axis=0) if strain_map[node.index] else np.zeros(3)
        result[node.index] = (σ_avg, ε_avg)

    return result  # node.index: (σ_vec, ε_vec)

def plot_results (estructure, elements, title, def_scale=1, force_scale=1e-2, reaction_scale=1e-2, sigma_y_tension=0, sigma_y_compression=0):
    f = estructure.f_original if hasattr(estructure, 'f_original') else estructure.f_global
    #plot_applied_forces(estructure.nodes, elements,title, f, scale=force_scale)

    # Importante: guardar los desplazamientos en cada nodo
    for node in estructure.nodes:
        node.structure = estructure  # para acceder a u_global desde cada nodo

    # Luego graficar
    #plot_deformed_structure(estructure.elements, title, scale=def_scale, show_ids=False)

    #reacciones = estructure.compute_reactions()

    #plot_deformed_with_reactions(title, estructure.elements, reacciones, scale=def_scale, reaction_scale=reaction_scale, show_ids=False)

    vm_nodal = compute_nodal_von_mises(estructure.elements, estructure.u_global)
    max_value = max(vm_nodal.values())
    min_value = min(vm_nodal.values())

    # Guardarlos en un archivo de texto
    with open('ENTREGA_2/QUAD9/resultados.txt', 'a') as file:
        file.write(f"Titulo: {title}\n")
        file.write(f"Máximo: {max_value}\n")
        file.write(f"Mínimo: {min_value}\n")
        file.write("\n")

    print("Máximo y Mínimo guardados en 'resultados.txt'.")
    plot_von_mises_field(estructure.nodes, estructure.elements, vm_nodal, title)
    plot_von_mises_field(estructure.nodes, estructure.elements, vm_nodal, title)
    #plot_von_mises_per_element(estructure.nodes, estructure.elements, vm_nodal, title)
    
    #nodal_fields = compute_nodal_stress_strain(estructure.nodes, estructure.elements, estructure.u_global)   
    
    #plot_all_scalar_fields_separately(estructure.nodes, estructure.elements, nodal_fields, title)
    plot_principal_fields(estructure.nodes, estructure.elements, estructure.u_global, title_prefix=title, sigma_y_tension=sigma_y_tension, sigma_y_compression=sigma_y_compression)