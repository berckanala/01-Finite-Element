import numpy as np

class Element:
    def __init__(self, id, node_ids, n_nodes):
        self.id = id
        self.node_ids = node_ids  # Lista de 6 nodos (orden: 1-2-3-4-5-6)
        self.K = None  # Se calculará con get_stiffness_matrix

    def shape_function_derivatives(self, xi, eta):
        L1 = 1 - xi - eta
        L2 = xi
        L3 = eta

        dN_dxi = np.array([
            -1,          # ∂N1/∂xi
            1-2*L3,           # ∂N2/∂xi
            4*L3,                  # ∂N3/∂xi
            -2*L3       # ∂N4/∂xi
            
        ])

        dN_deta = np.array([
            -1,          # ∂N1/∂eta
            -2*L2,                  # ∂N2/∂eta
            4*L2,           # ∂N3/∂eta
            1-2*L2           # ∂N4/∂eta
           
        ])

        return dN_dxi, dN_deta

    def get_B_matrix(self, nodes, xi, eta):
        coords = np.array([[nodes[i - 1].x, nodes[i - 1].y] for i in self.node_ids])  # (4,2)
        dN_dxi, dN_deta = self.shape_function_derivatives(xi, eta)

        J = np.zeros((2, 2))
        for i in range(4):
            J[0, 0] += dN_dxi[i] * coords[i, 0]
            J[0, 1] += dN_dxi[i] * coords[i, 1]
            J[1, 0] += dN_deta[i] * coords[i, 0]
            J[1, 1] += dN_deta[i] * coords[i, 1]

        detJ = np.linalg.det(J)
        
        if abs(detJ) < 1e-12:
            return np.zeros((2, 4)), 0.0

        Jinv = np.linalg.inv(J)

        dN_dx = Jinv[0, 0] * dN_dxi + Jinv[0, 1] * dN_deta
        dN_dy = Jinv[1, 0] * dN_dxi + Jinv[1, 1] * dN_deta

        B = np.vstack((dN_dx, dN_dy))

        return B, detJ

    def get_stiffness_matrix(self, nodes):
        gauss_points = [
            (1/6, 1/6, 1/6),
            (2/3, 1/6, 1/6),
            (1/6, 2/3, 1/6)
        ]

        K = np.zeros((4, 4))
        I = np.identity(2)

        for xi, eta, w in gauss_points:
            B, detJ = self.get_B_matrix(nodes, xi, eta)
            if detJ == 0:
                continue
            K += w * detJ * (B.T @ I @ B)

        self.K = K
        return self.K

    def get_load_vector(self, nodes, alpha):
        coords = np.array([[nodes[i - 1].x, nodes[i - 1].y] for i in self.node_ids])
        gauss_points = [
            (1/6, 1/6, 1/6),
            (2/3, 1/6, 1/6),
            (1/6, 2/3, 1/6)
        ]

        f_local = np.zeros(4)

        for xi, eta, w in gauss_points:
            N = self.shape_functions(xi, eta)
            x = np.dot(N, coords[:, 0])
            y = np.dot(N, coords[:, 1])
            r2 = x**2 + y**2

            f_val = 0.0
            if not (r2 == 0 and alpha < 1):
                f_val = -alpha**2 * r2**(alpha / 2 - 1)

            _, detJ = self.get_B_matrix(nodes, xi, eta)
            f_local += f_val * w * detJ * N

        return f_local

    def shape_functions(self, xi, eta):
        # Funciones de forma para un elemento cuadrilátero de 4 nodos
        N1 = (1 - xi) * (1 - eta) / 4
        N2 = (1 + xi) * (1 - eta) / 4
        N3 = (1 + xi) * (1 + eta) / 4
        N4 = (1 - xi) * (1 + eta) / 4

        N = np.array([N1, N2, N3, N4])
        return N