import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

class Solve:
    def __init__(self, nodes, elements, alpha):
        self.nodes = nodes
        self.elements = elements  # LST elements
        self.alpha = alpha
        self.nnodes = len(nodes)

        self.make_elements_stiffness_matrices()
        self.K_global = self.assemble_global_matrix()
        self.f = self.assemble_global_load_vector()  # <--- Aquí usamos el vector f real
        self.u = np.zeros(self.nnodes)

    def make_elements_stiffness_matrices(self):
        for element in self.elements:
            element.get_stiffness_matrix(self.nodes)

    def assemble_global_matrix(self):
        row_idx = []
        col_idx = []
        data = []

        for elem in self.elements:
            K_local = elem.K
            node_ids = elem.node_ids  # deben ser enteros base-1

            for i_local in range(6):
                for j_local in range(6):
                    global_i = node_ids[i_local] - 1
                    global_j = node_ids[j_local] - 1
                    row_idx.append(global_i)
                    col_idx.append(global_j)
                    data.append(K_local[i_local, j_local])

        return coo_matrix((data, (row_idx, col_idx)), shape=(self.nnodes, self.nnodes)).tocsr()

    def assemble_global_load_vector(self):
        f_global = np.zeros(self.nnodes)

        for elem in self.elements:
            f_local = elem.get_load_vector(self.nodes, self.alpha)  # <-- requiere estar definido en LST
            for i_local, node_id in enumerate(elem.node_ids):
                global_i = node_id - 1
                f_global[global_i] += f_local[i_local]

        return f_global

    def solve_matrix(self):
        fixed_dofs = []
        fixed_values = []

        for node in self.nodes:
            if hasattr(node, "boundary_label") and any("Dirichlet" in label for label in node.boundary_label):
                fixed_dofs.append(node.id - 1)
                fixed_values.append(node.u)

        fixed_dofs = np.array(fixed_dofs, dtype=int)
        fixed_values = np.array(fixed_values, dtype=float)

        all_dofs = np.arange(self.nnodes)
        free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

        K = self.K_global
        f = self.f

        f_reduced = f[free_dofs] - K[free_dofs][:, fixed_dofs] @ fixed_values
        u_free = spsolve(K[free_dofs][:, free_dofs], f_reduced)

        u_full = np.zeros(self.nnodes)
        u_full[fixed_dofs] = fixed_values
        u_full[free_dofs] = u_free

        self.u = u_full

        for node in self.nodes:
            node.u_fem = u_full[node.id - 1]

    def real_solution(self):
        for node in self.nodes:
            node.solve_u(self.alpha)

    def semi_norm_H1_0(self, alpha, orden=5):
        """
        Calcula |u|^2_{H^1_0(Ω)} = ∫_Ω |∇u(x,y)|² dxdy, donde u = (x² + y²)^{α/2},
        y Ω = [0,1] × [0,1] usando cuadratura de Gauss-Legendre.
        """
        puntos, pesos = np.polynomial.legendre.leggauss(orden)
        puntos = 0.5 * (puntos + 1)
        pesos = 0.5 * pesos

        total = 0.0

        for i in range(orden):
            for j in range(orden):
                x = puntos[i]
                y = puntos[j]
                w = pesos[i] * pesos[j]

                r2 = x**2 + y**2
                if r2 == 0 and alpha < 1:
                    grad2 = 0.0  # evitar singularidad
                else:
                    grad2 = alpha**2 * r2**(alpha - 1)

                total += grad2 * w

        return total  # ya es la semi-norma al cuadrado

    def femm_solution(self):
        x = np.zeros(self.nnodes)
        for node in self.nodes:
            if hasattr(node, "boundary_label") and any("Dirichlet" in label for label in node.boundary_label):
                x[node.id - 1] = node.u
            else:
                x[node.id - 1] = node.u_fem

        return x @ self.K_global @ x