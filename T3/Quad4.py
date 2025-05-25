import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class Quad4:
    def __init__(self, element_tag, node_list, section, load_direction=None, material=None, type='planeStress', print_summary=False):
        if len(node_list) != 4:
            raise ValueError("Quad4 elements must have exactly 4 nodes.")
        
        self.node_list = node_list  # Asegúrate de asignar node_list a la propiedad de la clase
        self.area = self.compute_area()  # Calcular el área al inicializar el objeto
        self.element_tag = element_tag
        self.nodes = node_list
        self.section = section
        self.load_direction = load_direction
        self.type = type
        self.gauss_points = [(-1/np.sqrt(3), -1/np.sqrt(3)), (1/np.sqrt(3), -1/np.sqrt(3)),
                             (1/np.sqrt(3), 1/np.sqrt(3)), (-1/np.sqrt(3), 1/np.sqrt(3))]

        self.idx = self.calculate_indices()
        self.kg = self.get_stiffness_matrix()

        if print_summary:
            self.printSummary()

    def __str__(self):
        return f"Quad4 Element {self.element_tag}: Nodes {[node.name for node in self.nodes]}"
    
    def compute_area(self):
        """
        Calcula el área del cuadrilátero utilizando el determinante basado en las coordenadas de los nodos.
        Formula: Area = 0.5 * | x1(y2 - y3) + x2(y3 - y1) + x3(y1 - y2) + x4(y2 - y1) |
        """
        # Obtener las coordenadas de los nodos
        x1, y1 = self.node_list[0].coordenadas
        x2, y2 = self.node_list[1].coordenadas
        x3, y3 = self.node_list[2].coordenadas
        x4, y4 = self.node_list[3].coordenadas
        
        # Calcular el área usando la fórmula de determinante
        area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2) + x4 * (y2 - y1))
        
        return area
    
    def __repr__(self):
        return self.__str__()

    def calculate_indices(self):
        return np.hstack([node.idx for node in self.nodes]).astype(int)

    def get_xy_matrix(self):
        return np.array([node.coordenadas for node in self.nodes])

    def get_shape_function_derivatives(self, xi, eta):
        dN_dxi = np.array([
            [-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)],
            [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]
        ]) * 0.25
        return dN_dxi

    def get_B_matrix(self, xi, eta):
        dN_dxi = self.get_shape_function_derivatives(xi, eta)
        coords = self.get_xy_matrix()
        J = dN_dxi @ coords  # Jacobiano 2x2

        if np.linalg.det(J) <= 0:
            raise ValueError(f"Elemento {self.element_tag} tiene Jacobiano no positivo.")

        J_inv = np.linalg.inv(J)
        dN_dx = J_inv @ dN_dxi

        B = np.zeros((3, 8))
        for i in range(4):
            B[0, 2*i]     = dN_dx[0, i]
            B[1, 2*i + 1] = dN_dx[1, i]
            B[2, 2*i]     = dN_dx[1, i]
            B[2, 2*i + 1] = dN_dx[0, i]
        return B, np.linalg.det(J)

    def get_stiffness_matrix(self):
        D = self.section.get_Emat(self.type)
        t = self.section.thickness
        Ke = np.zeros((8, 8))

        for xi, eta in self.gauss_points:
            B, detJ = self.get_B_matrix(xi, eta)
            Ke += B.T @ D @ B * detJ * t
        return Ke

    def get_centroid(self):
        return np.mean(self.get_xy_matrix(), axis=0)

    def get_element_displacements(self, u):
        return u[self.idx]

    def get_element_strains(self, u):
        ue = self.get_element_displacements(u)
        B, _ = self.get_B_matrix(0.0, 0.0)  # en el centro
        epsilon = B @ ue
        return epsilon, ue

    def get_element_stress(self, u):
        epsilon, ue = self.get_element_strains(u)
        D = self.section.material.get_Emat(self.type)
        sigma = D @ epsilon
        return sigma, epsilon, ue

    def calculate_principal_stress(self, sigma):
        sx, sy, sxy = sigma
        S = np.array([[sx, sxy], [sxy, sy]])
        vals, _ = np.linalg.eig(S)
        return np.sort(vals)[::-1]

    def calculate_principal_strain(self, epsilon):
        ex, ey, exy = epsilon
        E = np.array([[ex, exy], [exy, ey]])
        vals, _ = np.linalg.eig(E)
        return np.sort(vals)[::-1]

    def get_element_internal_forces(self, u):
        ue = self.get_element_displacements(u)
        return self.kg @ ue

    def get_results(self, u):
        sigma, epsilon, ue = self.get_element_stress(u)
        fe = self.get_element_internal_forces(u)
        sigma_p = self.calculate_principal_stress(sigma)
        epsilon_p = self.calculate_principal_strain(epsilon)

        return {
            'stress': sigma,
            'strain': epsilon,
            'displacement': ue,
            'internal_forces': fe,
            'principal_stress': sigma_p,
            'principal_strain': epsilon_p
        }
    def body_weight_forces(self, x, y, force_vector):
        """
        Calcula las fuerzas internas por el peso propio para un elemento cuadrilátero.
        
        Args:
            x (float): Coordenada x del centroide del elemento
            y (float): Coordenada y del centroide del elemento
            force_vector (list): Fuerza total (en N) en la dirección [Fx, Fy]
        
        Returns:
            List: Fuerzas internas distribuidas para cada nodo
        """
        # Calculamos las fuerzas internas distribuidas (en N) para los 4 nodos
        # Se supone que la carga se distribuye uniformemente entre los nodos
        return np.array(force_vector)  # Retornamos el vector de fuerzas para el elemento
    def plotGeometry(self, ax=None, text=False, nodes=True, nodeLabels=False, facecolor='lightblue', edgecolor='k', alpha=0.5):
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
        print(f'-------------------------------------------------------------')
        print(f"Quad4 Element {self.element_tag}")
        print(f"Type: {self.type}")
        print(f"Nodes: {[node.name for node in self.nodes]}")
        
        coords = self.get_xy_matrix()
        for i, node in enumerate(self.nodes):
            print(f"  Node {node.name}: ({coords[i,0]:.3f}, {coords[i,1]:.3f})")
        
        print(f"Thickness: {self.section.thickness}")
        print(f"Element DoF indices: {self.calculate_indices()}")
        
        if self.load_direction is not None:
            print(f"Body force direction: {self.load_direction}")
        else:
            print(f"Body force direction: None")
        
        print(f"\nStiffness matrix (local):\n{self.get_stiffness_matrix()}")
        print(f'-------------------------------------------------------------\n')
