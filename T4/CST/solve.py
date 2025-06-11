import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

class Solver:
    def __init__(self, mesh_nodes, triangle_elements, alpha_param):
        # Guardar nodos, elementos y parámetro alpha
        self.mesh_nodes = mesh_nodes
        self.triangle_elements = triangle_elements
        self.alpha_param = alpha_param
        self.num_nodes = len(mesh_nodes)

        # Ensamblar matrices de rigidez locales y globales
        self.compute_local_stiffness_matrices()
        self.global_stiffness_matrix = self.assemble_global_matrix()

        # Inicializar solución y vector de fuerzas
        self.solution_vector = np.zeros(self.num_nodes)
        self.force_vector = np.zeros(self.num_nodes)

    def compute_local_stiffness_matrices(self):
        # Calcular la matriz de rigidez local para cada elemento triangular
        for element in self.triangle_elements:
            element.get_stiffness_matrix(self.mesh_nodes)

    def assemble_global_matrix(self):
        # Ensamblar la matriz de rigidez global usando formato disperso
        row_idx = []
        col_idx = []
        data = []

        for elem in self.triangle_elements:
            node_ids = [int(i) for i in elem.node_ids]  # IDs base 1
            K_local = elem.K

            for i_local in range(3):
                for j_local in range(3):
                    global_i = node_ids[i_local] - 1
                    global_j = node_ids[j_local] - 1

                    row_idx.append(global_i)
                    col_idx.append(global_j)
                    data.append(K_local[i_local, j_local])

        return coo_matrix((data, (row_idx, col_idx)), shape=(self.num_nodes, self.num_nodes)).tocsr()

    def solve_system(self):
        # Resolver el sistema de ecuaciones aplicando condiciones de Dirichlet
        fixed_dofs = []
        fixed_values = []

        for node in self.mesh_nodes:
            if hasattr(node, "boundary_label") and any("Diritchlet" in label for label in node.boundary_label):
                fixed_dofs.append(node.id - 1)
                fixed_values.append(node.u)

        fixed_dofs = np.array(fixed_dofs, dtype=int)
        fixed_values = np.array(fixed_values, dtype=float)

        all_dofs = np.arange(self.num_nodes)
        free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

        K = self.global_stiffness_matrix
        f = self.force_vector.copy()

        f_reduced = f[free_dofs] - K[free_dofs][:, fixed_dofs] @ fixed_values
        u_free = spsolve(K[free_dofs][:, free_dofs], f_reduced)

        u_full = np.zeros(self.num_nodes)
        u_full[fixed_dofs] = fixed_values
        u_full[free_dofs] = u_free

        self.solution_vector = u_full
        for node in self.mesh_nodes:
            node.u_fem = u_full[node.id - 1]

    def assign_exact_solution(self):
        # Asignar solución exacta a cada nodo (basada en la función conocida)
        for node in self.mesh_nodes:
            node.solve_u(self.alpha_param)

    def compute_H1_seminorm_exact(self):
        # Calcular ∫Ω |∇u|² para la solución exacta, usando cuadratura de orden 3
        alpha = self.alpha_param
        order = 3

        puntos, pesos = np.polynomial.legendre.leggauss(order)
        puntos = 0.5 * (puntos + 1)
        pesos = 0.5 * pesos

        total = 0.0

        for i in range(order):
            for j in range(order):
                x = puntos[i]
                y = puntos[j]
                w = pesos[i] * pesos[j]

                r2 = x**2 + y**2
                if r2 == 0 and alpha < 1:
                    grad2 = 0.0  # evitar singularidad en el origen
                else:
                    grad2 = alpha**2 * r2**(alpha - 1)

                total += grad2 * w

        return total

    def compute_fem_energy(self):
        # Calcular a(u_h, u_h) = u_hᵀ K u_h
        x = np.zeros(self.num_nodes)

        for node in self.mesh_nodes:
            if hasattr(node, "boundary_label") and any("Diritchlet" in label for label in node.boundary_label):
                x[node.id - 1] = node.u
            else:
                x[node.id - 1] = node.u_fem

        K = self.global_stiffness_matrix
        fem_energy = x.T @ K @ x

        return fem_energy
