import numpy as np
from numpy.polynomial.legendre import leggauss

class Quad4:
    def __init__(self, 
                 element_tag: int, 
                 node_list: list, 
                 section: object, 
                 type: str = 'planeStress',
                 n_gauss: int = 2,
                 load_direction=None,
                 point_load=None,
                 print_summary: bool = False):

        if len(node_list) != 4:
            raise ValueError("Quad4 elements must have exactly 4 nodes.")
        
        self.element_tag = element_tag
        self.node_list = node_list
        self.nodes = node_list
        self.section = section
        self.type = type
        self.n_gauss = n_gauss
        self.load_direction = load_direction
        self.point_load = point_load

        self.coords = self.get_xy_matrix()
        self.area = self.compute_element_area()
        self.idx = self.calculate_indices()
        self.kg = self.get_stiffness_matrix()
        self.f_eq = self.get_body_forces()

        if print_summary:
            self.printSummary()

    def __str__(self):
        return f"Quad4 Element {self.element_tag}: Nodes {[node.name for node in self.nodes]}"

    def get_xy_matrix(self):
        return np.array([node.coordenadas for node in self.nodes])

    def calculate_indices(self):
        return np.hstack([node.idx for node in self.nodes])

    def get_centroid(self):
        return np.mean(self.coords, axis=0)

    def get_stiffness_matrix(self):
        E_mat = self.section.get_Emat(self.type)
        t = self.section.thickness
        points, weights = self.gauss_points(self.n_gauss)
        K = np.zeros((8, 8))

        for (xi, eta), w in zip(points, weights):
            B = self.strain_displacement_matrix(xi, eta)
            J = self.jacobian_matrix(xi, eta)
            detJ = self.jacobian_determinant(J)
            K += B.T @ E_mat @ B * detJ * t * w

        return K

    def plotGeometry(self, ax=None, text=False, nodes=True, nodeLabels=False, facecolor='lightgray', edgecolor='k', alpha=0.5):
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches

        if ax is None:
            fig, ax = plt.subplots()

        polygon = patches.Polygon(self.coords, closed=True, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)
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
        
        for i, node in enumerate(self.nodes):
            print(f"  Node {node.name}: ({self.coords[i, 0]:.3f}, {self.coords[i, 1]:.3f})")
        
        print(f"Thickness: {self.section.thickness}")
        print(f"Area: {self.area:.4f}")
        print(f"Element DoF indices: {self.idx}")
        print(f"\nStiffness matrix (local):\n{self.kg}")
        print(f'-------------------------------------------------------------\n')

    # ---------------------------
    # Métodos auxiliares internos
    # ---------------------------
    @staticmethod
    def shape_functions(xi, eta):
        N1 = 0.25 * (1 - xi) * (1 - eta)
        N2 = 0.25 * (1 + xi) * (1 - eta)
        N3 = 0.25 * (1 + xi) * (1 + eta)
        N4 = 0.25 * (1 - xi) * (1 + eta)
        return np.array([[N1], [N2], [N3], [N4]])

    @staticmethod
    def shape_function_derivatives(xi, eta):
        dN_dxi = np.array([
            [-0.25 * (1 - eta)],
            [ 0.25 * (1 - eta)],
            [ 0.25 * (1 + eta)],
            [-0.25 * (1 + eta)]
        ]).flatten()

        dN_deta = np.array([
            [-0.25 * (1 - xi)],
            [-0.25 * (1 + xi)],
            [ 0.25 * (1 + xi)],
            [ 0.25 * (1 - xi)]
        ]).flatten()

        return dN_dxi, dN_deta

    def jacobian_matrix(self, xi, eta):
        dN_dxi, dN_deta = self.shape_function_derivatives(xi, eta)
        J = np.zeros((2, 2))
        for i in range(4):
            J[0, 0] += dN_dxi[i] * self.coords[i, 0]
            J[0, 1] += dN_dxi[i] * self.coords[i, 1]
            J[1, 0] += dN_deta[i] * self.coords[i, 0]
            J[1, 1] += dN_deta[i] * self.coords[i, 1]
        return J

    def jacobian_determinant(self, J):
        return np.linalg.det(J)

    def shape_function_derivatives_global(self, xi, eta):
        dN_dxi, dN_deta = self.shape_function_derivatives(xi, eta)
        J = self.jacobian_matrix(xi, eta)
        J_inv = np.linalg.inv(J)

        dN_dx = []
        dN_dy = []
        for i in range(4):
            dN_nat = np.array([dN_dxi[i], dN_deta[i]])
            dN_phys = J_inv @ dN_nat
            dN_dx.append(dN_phys[0])
            dN_dy.append(dN_phys[1])
        return np.array(dN_dx), np.array(dN_dy)

    def strain_displacement_matrix(self, xi, eta):
        dN_dx, dN_dy = self.shape_function_derivatives_global(xi, eta)
        B = np.zeros((3, 8))
        for i in range(4):
            B[0, 2*i]     = dN_dx[i]
            B[1, 2*i + 1] = dN_dy[i]
            B[2, 2*i]     = dN_dy[i]
            B[2, 2*i + 1] = dN_dx[i]
        return B

    def generalized_interpolation(self, xi, eta, data_matrix):
        N = self.shape_functions(xi, eta)
        return (data_matrix @ N).flatten()

    def compute_element_area(self):
        x = self.coords[:, 0]
        y = self.coords[:, 1]
        area = 0.5 * abs(
            x[0]*y[1] + x[1]*y[2] + x[2]*y[3] + x[3]*y[0] -
            (y[0]*x[1] + y[1]*x[2] + y[2]*x[3] + y[3]*x[0])
        )
        return area

    @staticmethod
    def gauss_points(n):
        pts_1d, w_1d = leggauss(n)
        Ξ, Η = np.meshgrid(pts_1d, pts_1d)
        Wξ, Wη = np.meshgrid(w_1d, w_1d)
        points = np.column_stack((Ξ.flatten(), Η.flatten()))
        weights = (Wξ * Wη).flatten()
        return points, weights

    def get_body_forces(self):
        if self.load_direction is not None:
            return self.body_forces(self.load_direction)
        if self.point_load is not None:
            x, y, fx, fy = self.point_load
            return self.apply_point_load(x, y, [fx, fy])
        return np.zeros(8)

    def body_forces(self, body_force_vector):
        bx, by = body_force_vector
        t = self.section.thickness
        A = self.area
        f_eq = (t * A / 4) * np.array([bx, by] * 4)
        return f_eq

    def apply_point_load(self, x, y, force_vector):
        xi = eta = 0.0  # Aproximación: punto cerca del centro
        N = self.shape_functions(xi, eta).flatten()
        fx, fy = force_vector
        f_p = np.zeros(8)
        for i in range(4):
            f_p[2*i]     = fx * N[i]
            f_p[2*i + 1] = fy * N[i]

        print("\nCarga puntual aplicada:")
        print(f"Ubicación estimada: (xi, eta) = (0, 0)")
        print(f"Fuerza: {force_vector}")
        print(f"Fuerzas nodales equivalentes: {np.round(f_p, 3)}")
        return f_p
