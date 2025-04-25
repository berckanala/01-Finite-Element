import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

globalParameters = {
    'nDoF': 2
}

class CST:
    def __init__(self, 
                 element_tag: int, 
                 node_list: list, 
                 section: object, 
                 load_direction=None,
                 point_load=None,
                 type: str = 'planeStress',
                 print_summary=False):
        if len(node_list) != 3:
            raise ValueError("CST elements must have exactly 3 nodes.")
        
        self.element_tag = element_tag
        self.node_list = node_list
        self.nodes = node_list
        self.section = section
        self.load_direction = load_direction
        self.point_load = point_load  # New variable for point load
        self.type = type

        self.compute_area()
        self.idx = self.calculate_indices()
        self.kg = self.get_stiffness_matrix()
        
        if print_summary is True:
            self.printSummary()

    def __str__(self):
        return f"CST Element {self.element_tag}: Nodes {[node.name for node in self.nodes]}"

    def calculate_indices(self):
        """Returns the global DoF indices for the CST element."""
        idx = np.hstack([node.idx for node in self.nodes])
        return idx

    def get_xy_matrix(self):
        """
        Returns the matrix of nodal coordinates of the CST element.

        Returns:
            X (np.ndarray): 3x2 array with node coordinates [[x1, y1], [x2, y2], [x3, y3]]
        """
        xy = np.array([node.coordenadas for node in self.nodes])
        return xy

    def get_centroid(self):
        """
        Computes the centroid of the triangular element using matrix operations.

        Returns:
            centroid (np.ndarray): (2,) array with centroid coordinates [x, y]
        """
        xy = self.get_xy_matrix()
        w = np.ones((1, 3)) / 3
        centroid = w @ xy
        return centroid.flatten()

    def compute_area(self):
        """
        Computes and stores the area of the CST element using determinant-based formula.

        Sets:
            self.area (float): Area of the triangle
        """
        x1, y1 = self.nodes[0].coordenadas
        x2, y2 = self.nodes[1].coordenadas
        x3, y3 = self.nodes[2].coordenadas

        self.area = 0.5 * np.linalg.det(np.array([
            [1, x1, y1],
            [1, x2, y2],
            [1, x3, y3]
        ]))

        if self.area <= 0:
            raise ValueError(f"Element {self.element_tag} has non-positive area: {self.area}")

    def get_interpolation_matrix(self, x: float, y: float):
        """
        Returns the interpolation matrix N at a given point (x, y)
        using matrix operations to compute barycentric coordinates.

        Returns:
            Nmat (2x6 np.array): Interpolation matrix
        """
        coords = self.get_xy_matrix()
        x1, y1 = coords[0]
        x2, y2 = coords[1]
        x3, y3 = coords[2]

        T = np.array([
            [1, x1, y1],
            [1, x2, y2],
            [1, x3, y3]
        ])

        p = np.array([1, x, y])

        try:
            lambdas = np.linalg.solve(T, p)
        except np.linalg.LinAlgError:
            raise ValueError("Singular element (zero area or colinear nodes)")

        Nmat = np.zeros((2, 6))
        Nmat[0, 0::2] = lambdas
        Nmat[1, 1::2] = lambdas

        return Nmat

    def get_B_matrix(self):
        """
        Computes the strain-displacement matrix B for the CST element.

        Returns:
            B (3x6 np.ndarray): Strain-displacement matrix
        """
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
        """
        Computes the local stiffness matrix for the CST element.

        Returns:
            Ke (6x6 np.array): Local stiffness matrix
        """
        D = self.section.get_Emat(self.type)
        B = self.get_B_matrix()
        t = self.section.thickness

        Ke = B.T @ D @ B * self.area * t
        return Ke
    
    def get_body_forces(self):
        """
        Calculates the equivalent body forces based on the type of load: distributed or point load.

        Returns:
            f_local (np.ndarray): The equivalent nodal forces due to the applied load (distributed or point load).
        """
        if self.load_direction is not None:
            # If there is a body force (distributed load), use the body_forces method
            return self.body_forces(self.load_direction)
        
        if self.point_load is not None:
            # If there is a point load, use apply_point_body_force
            x, y, fx, fy = self.point_load
            return self.apply_point_body_force(x, y, [fx, fy])

        return np.zeros(6)  # No load applied, return zero forces

    def body_forces(self, body_force_vector):
        """
        Calculates the equivalent nodal forces for a constant distributed body force.

        Args:
            body_force_vector (list): The body force vector [bx, by].

        Returns:
            f_body (np.ndarray): The equivalent nodal forces due to the distributed body force.
        """
        t = self.section.thickness
        A = self.area
        bx, by = body_force_vector

        # Vector of equivalent force due to the constant distributed load
        f_body = (t * A / 3) * np.array([bx, by, bx, by, bx, by])
        return f_body
    
    def apply_point_body_force(self, x, y, force_vector):
        """
        Calculates the equivalent nodal forces for a point load applied at a specific point in the element.

        Args:
            x (float): x-coordinate of the point load.
            y (float): y-coordinate of the point load.
            force_vector (list): The force vector [fx, fy].

        Returns:
            f_puntual (np.ndarray): The equivalent nodal forces due to the point load.
        """
        N = self.get_interpolation_matrix(x, y)
        fx, fy = force_vector
        f_puntual = (N.T @ np.array([fx, fy])).flatten()

        # Save for later use in visualization
        self.body_point = (x, y)
        self.body_vector = (fx, fy)

        print("\nPoint force applied in the element:")
        print(f"Location: ({x}, {y})")
        print(f"Force: {force_vector}")
        print(f"Equivalent nodal forces: {np.round(f_puntual, 2)}")

        return f_puntual
    
            
    def plotGeometry(self, ax=None, text=False, nodes=True, nodeLabels=False, facecolor='lightgray', edgecolor='k', alpha=0.5):

        if ax is None:
            fig, ax = plt.subplots()

        coords = self.get_xy_matrix()  # (3,2)

        # Crear y agregar el parche del triángulo
        polygon = patches.Polygon(coords, closed=True, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)
        ax.add_patch(polygon)

        # Dibujar los nodos
        if nodes or nodeLabels:
            for node in self.nodes:
                node.plotGeometry(ax, text=nodeLabels)

        # Mostrar el número del elemento en el centroide
        if text:
            x_c, y_c = self.get_centroid()
            ax.text(x_c, y_c, f'{self.element_tag}', fontsize=12, ha='center', va='center')

        return ax
    
    def printSummary(self):
        """
        Prints a detailed summary of the CST element.
        """
        print(f'-------------------------------------------------------------')
        print(f"CST Element {self.element_tag}")
        print(f"Type: {self.type}")
        print(f"Nodes: {[node.name for node in self.nodes]}")
        
        coords = self.get_xy_matrix()
        for i, node in enumerate(self.nodes):
            print(f"  Node {node.name}: ({coords[i,0]:.3f}, {coords[i,1]:.3f})")
        
        print(f"Thickness: {self.section.thickness}")
        print(f"Area: {self.area:.4f}")
        print(f"Element DoF indices: {self.calculate_indices()}")
        
        if self.load_direction is not None:
            print(f"Body force direction: {self.load_direction}")
        else:
            print(f"Body force direction: None")
        
        print(f"\nStiffness matrix (local):\n{self.get_stiffness_matrix()}")
        print(f'-------------------------------------------------------------\n')

