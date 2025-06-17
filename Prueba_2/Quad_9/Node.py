import numpy as np
import matplotlib.pyplot as plt
import os



class Node:
    def __init__(self, id, coord=[0,0], dofs=None, restrain=[0, 0], name=None):
        self.index = id
        self.coord = coord
        self.name = name if name is not None else f"Node_{id}"  # Asigna un nombre por defecto
        if dofs is None:
            self.dofs = np.array([(id * 2)-1, id * 2 ])
            

        else:
            self.dofs = np.array(dofs)
        self.restrain = np.array(restrain)



    def plot_nodes_por_grupo(grupo_nodos_dict, title, show_ids=False, save=True):
        """
        Plotea nodos por grupo con colores diferentes.

        Parámetros:
        - grupo_nodos_dict: dict con nombre del grupo como clave y lista de nodos como valor
        - title: nombre del gráfico y del archivo a guardar
        - show_ids: muestra los índices si True
        - save: guarda el gráfico como PNG si True, lo muestra en pantalla si False
        """
        all_x = []
        all_y = []

        for nodos in grupo_nodos_dict.values():
            all_x.extend([n.coord[0] for n in nodos])
            all_y.extend([n.coord[1] for n in nodos])

        x_min, x_max = min(all_x), max(all_x)
        y_min, y_max = min(all_y), max(all_y)
        x_margin = (x_max - x_min) * 0.05
        y_margin = (y_max - y_min) * 0.05

        x_range = (x_max - x_min) + 2 * x_margin
        y_range = (y_max - y_min) + 2 * y_margin

        fixed_width = 7  # pulgadas
        aspect_ratio = y_range / x_range
        height = fixed_width * aspect_ratio

        fig, ax = plt.subplots(figsize=(fixed_width, height))
        colors = plt.cm.get_cmap("tab10", len(grupo_nodos_dict))

        for i, (grupo, nodos) in enumerate(grupo_nodos_dict.items()):
            xs = [n.coord[0] for n in nodos]
            ys = [n.coord[1] for n in nodos]
            ax.scatter(xs, ys, s=5, color=colors(i), label=grupo)

            if show_ids:
                for node in nodos:
                    ax.text(node.coord[0] + 0.5, node.coord[1] + 0.5, str(node.index), fontsize=6)

        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("Nodes per Group")
        ax.grid(True)
        # ax.legend()

        if save:
            os.makedirs("INFORME/GRAFICOS", exist_ok=True)
            fig.savefig(f"INFORME/GRAFICOS/{title}_nodes_por_grupo.png", dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()



   