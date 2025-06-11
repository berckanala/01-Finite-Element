import numpy as np

class LST:
    def __init__(self, id, node_ids):
        self.id = id
        self.node_ids = node_ids  # Lista de 6 nodos (orden: 1-2-3-4-5-6)
        self.K = None  # Se calculará con get_stiffness_matrix

    def shape_function_derivatives(self, xi, eta):
        """
        Derivadas de las funciones de forma cuadráticas respecto a xi y eta.
        """
        L1 = 1 - xi - eta
        L2 = xi
        L3 = eta

        dN_dxi = np.array([
            4*L1 - 1,
            0,
            -4*L3 + 1,
            4*(L2 - L1),
            4*L3,
            -4*L2
        ])

        dN_deta = np.array([
            4*L1 - 1,
            -4*L2 + 1,
            0,
            -4*L1,
            4*L2,
            4*(L3 - L1)
        ])

        return dN_dxi, dN_deta

    def get_stiffness_matrix(self, nodes):
        """
        Calcula la matriz de rigidez del elemento LST usando integración numérica (orden 3),
        incluyendo validaciones de nodos repetidos y elementos degenerados.
        """
        # Extraer coordenadas de los nodos del elemento (base 1 → base 0)
        coords = np.array([[nodes[i - 1].x, nodes[i - 1].y] for i in self.node_ids])  # (6x2)

        # ⚠️ Verifica que los nodos no estén repetidos
        if len(set(self.node_ids)) < 6:
            print(f"⚠️ Elemento {self.id} tiene nodos repetidos: {self.node_ids}")
            self.K = np.zeros((6, 6))
            return self.K

        # ⚠️ Verifica que el triángulo principal no esté colapsado (área ~ 0)
        x1, y1 = coords[0]
        x2, y2 = coords[1]
        x3, y3 = coords[2]
        area = 0.5 * abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))
        if area < 1e-10:
            print(f"⚠️ Elemento {self.id} con área degenerada (≈ 0)")
            self.K = np.zeros((6, 6))
            return self.K

        # Puntos de Gauss y pesos para integración triangular (orden 3)
        gauss_pts = np.array([
            [1/6, 1/6],
            [2/3, 1/6],
            [1/6, 2/3]
        ])
        weights = np.array([1/6, 1/6, 1/6])

        K = np.zeros((6, 6))
        I = np.identity(2)

        for (xi, eta), w in zip(gauss_pts, weights):
            dN_dxi, dN_deta = self.shape_function_derivatives(xi, eta)

            # Matriz Jacobiana
            J = np.zeros((2, 2))
            for a in range(6):
                J[0, 0] += dN_dxi[a] * coords[a, 0]
                J[0, 1] += dN_dxi[a] * coords[a, 1]
                J[1, 0] += dN_deta[a] * coords[a, 0]
                J[1, 1] += dN_deta[a] * coords[a, 1]

            detJ = np.linalg.det(J)
            if detJ <= 0:
                print(f"⚠️ Elemento {self.id} tiene Jacobiano singular o negativo (det = {detJ:.2e})")
                self.K = np.zeros((6, 6))
                return self.K

            J_inv = np.linalg.inv(J)

            # Derivadas globales ∇N
            dN_dx = np.zeros((6, 2))
            for a in range(6):
                grad = J_inv @ np.array([dN_dxi[a], dN_deta[a]])
                dN_dx[a, :] = grad

            # Matriz B y ensamblaje
            B = dN_dx.T  # 2x6
            K += (B.T @ I @ B) * detJ * w

        if not np.any(K):
            print(f"⚠️ Elemento {self.id} matriz de rigidez completamente nula.")

        self.K = K
        return self.K

    def get_load_vector(self, nodes, alpha):
        """
        Calcula el vector de carga elemental ∫_T f(x,y) N_i(x,y) dA
        donde f(x,y) = -Δu = -α(α+2)(x² + y²)^{α/2 - 1}
        """
        node_coords = np.array([[nodes[nid - 1].x, nodes[nid - 1].y] for nid in self.node_ids])

        # Puntos y pesos para integración en triángulo (orden 3, exacto para polinomios grado ≤ 2)
        gauss_points = [
            (1/6, 1/6, 1/6),
            (2/3, 1/6, 1/6),
            (1/6, 2/3, 1/6)
        ]

        fe = np.zeros(6)

        for ξ, η, w in gauss_points:
            L1 = 1 - ξ - η
            L2 = ξ
            L3 = η

            # Coordenadas físicas (x, y)
            x = L1 * node_coords[0, 0] + L2 * node_coords[1, 0] + L3 * node_coords[2, 0]
            y = L1 * node_coords[0, 1] + L2 * node_coords[1, 1] + L3 * node_coords[2, 1]

            r2 = x**2 + y**2
            f = -alpha * (alpha + 2) * r2**(alpha / 2 - 1)

            # Funciones de forma de orden 2 (Lagrange 6 nodos)
            N = np.array([
                L1 * (2 * L1 - 1),
                L2 * (2 * L2 - 1),
                L3 * (2 * L3 - 1),
                4 * L1 * L2,
                4 * L2 * L3,
                4 * L3 * L1
            ])

            # Área del triángulo
            x1, y1 = node_coords[0]
            x2, y2 = node_coords[1]
            x3, y3 = node_coords[2]
            area = 0.5 * abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))

            fe += f * N * w * area * 2  # multiplicar por 2 ya que los pesos están normalizados

        return fe

