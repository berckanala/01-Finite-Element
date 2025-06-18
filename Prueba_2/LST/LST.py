import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class LST:
    def __init__(self, element_tag: int, node_list: list, section: object, load_direction=None, material=None, type: str = 'planeStress', print_summary=False):
        if len(node_list) < 3:
            raise ValueError("LST elements must have at least 3 nodes.")
        
        self.element_tag = element_tag
        self.node_list = node_list
        self.nodes = node_list
        self.section = section
        self.load_direction = load_direction
        self.type = type

        # Realizar cálculos para el área, rigidez y demás propiedades
        self.compute_area()
        self.idx = self.calculate_indices()
        self.kg = self.get_stiffness_matrix()
        
        if print_summary:
            self.printSummary()

    def __str__(self):
        return f"LST Element {self.element_tag}: Nodes {[node.name for node in self.nodes]}"

    def calculate_indices(self):
        """Devuelve los índices globales de los grados de libertad para el elemento LST."""
        idx = np.hstack([node.idx for node in self.nodes]).astype(int)
        return idx

    def get_xy_matrix(self):
        """Devuelve la matriz de coordenadas nodales del elemento LST."""
        xy = np.array([node.coordenadas for node in self.nodes])
        return xy

    def compute_area(self):
        """Calcula y almacena el área del elemento LST utilizando la fórmula del determinante."""
        x1, y1 = self.nodes[0].coordenadas
        x2, y2 = self.nodes[1].coordenadas
        x3, y3 = self.nodes[2].coordenadas

        # Calcular el área usando el determinante
        area_check = 0.5 * np.linalg.det(np.array([
            [1, x1, y1],
            [1, x2, y2],
            [1, x3, y3]
        ]))
        
        self.area=area_check
        if area_check == 0:
            raise ValueError(f"Element {self.element_tag} has collinear nodes, cannot compute area")

        if self.area <= 0:
            raise ValueError(f"Element {self.element_tag} has non-positive area: {self.area}")

    def get_B_matrix(self):
        """Calcula la matriz de deformación-desplazamiento B para el elemento LST."""
        x1, y1 = self.nodes[0].coordenadas
        x2, y2 = self.nodes[1].coordenadas
        x3, y3 = self.nodes[2].coordenadas

        b1 = y2 - y3
        b2 = y3 - y1
        b3 = y1 - y2
        c1 = x3 - x2
        c2 = x1 - x3
        c3 = x2 - x1

        B = (1 / (2 * self.area)) * np.array([
            [b1, 0, b2, 0, b3, 0],
            [0, c1, 0, c2, 0, c3],
            [c1, b1, c2, b2, c3, b3]
        ])
        return B

    def get_stiffness_matrix(self):
        """Calcula la matriz de rigidez local para el elemento LST."""
        D = self.section.get_Emat(self.type)
        B = self.get_B_matrix()
        t = self.section.thickness

        Ke = B.T @ D @ B * self.area * t
        return Ke

    def plotGeometry(self, ax=None, text=False, nodes=True, nodeLabels=False, facecolor='lightgray', edgecolor='k', alpha=0.5):
        """Visualiza la geometría del elemento LST como un triángulo sombreado."""
        if ax is None:
            fig, ax = plt.subplots()

        coords = self.get_xy_matrix()
        polygon = patches.Polygon(coords, closed=True, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)
        ax.add_patch(polygon)

        if nodes or nodeLabels:
            for node in self.nodes:
                node.plotGeometry(ax, text=nodeLabels)

        if text:
            x_c, y_c = self.get_centroid()
            ax.text(x_c, y_c, f'{self.element_tag}', fontsize=12, ha='center', va='center')

        return ax

    def printSummary(self):
        """Imprime un resumen detallado del elemento LST."""
        print(f'-------------------------------------------------------------')
        print(f"LST Element {self.element_tag}")
        print(f"Nodes: {[node.name for node in self.nodes]}")
        coords = self.get_xy_matrix()
        for i, node in enumerate(self.nodes):
            print(f"  Node {node.name}: ({coords[i,0]:.3f}, {coords[i,1]:.3f})")
        print(f"Thickness: {self.section.thickness}")
        print(f"Area: {self.area:.4f}")
        print(f"Element DoF indices: {self.calculate_indices()}")
        print(f'-------------------------------------------------------------\n')
