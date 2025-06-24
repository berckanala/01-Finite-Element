import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

class Solve:
    def __init__(self, nodes, elements, alpha):
        self.nodes = nodes
        self.elements = elements  # Lista de elementos
        self.alpha = alpha
        self.nnodes = len(nodes)

        # Crear la matriz global de rigidez y el vector de carga
        self.make_elements_stiffness_matrices()
        self.K_global = self.assemble_global_matrix()
        self.f = self.assemble_global_load_vector()  # Vector de carga global
        self.u = np.zeros(self.nnodes)

    def make_elements_stiffness_matrices(self):
        # Calcular las matrices de rigidez para cada elemento
        for element in self.elements:
            element.get_stiffness_matrix(self.nodes)

    def assemble_global_matrix(self):
        # Ensamblar la matriz global de rigidez
        row_idx = []
        col_idx = []
        data = []

        for elem in self.elements:
            K_local = elem.K
            node_ids = elem.node_ids  # Deben ser enteros base-1

            # Añadir los valores de la matriz de rigidez local a la matriz global
            for i_local in range(len(node_ids)):
                for j_local in range(len(node_ids)):
                    global_i = node_ids[i_local] - 1  # Convertir a índice base-0
                    global_j = node_ids[j_local] - 1  # Convertir a índice base-0
                    row_idx.append(global_i)
                    col_idx.append(global_j)
                    data.append(K_local[i_local, j_local])

        # Crear la matriz dispersa en formato COO y convertir a CSR
        return coo_matrix((data, (row_idx, col_idx)), shape=(self.nnodes, self.nnodes)).tocsr()

    def assemble_global_load_vector(self):
        # Ensamblar el vector global de carga
        f_global = np.zeros(self.nnodes)

        for elem in self.elements:
            f_local = elem.get_load_vector(self.nodes, self.alpha)  # Vector de carga local
            for i_local, node_id in enumerate(elem.node_ids):
                global_i = node_id - 1  # Convertir a índice base-0
                f_global[global_i] += f_local[i_local]  # Añadir la carga local al global

        return f_global

    def solve_matrix(self):
        # Resolver el sistema de ecuaciones considerando condiciones de frontera
        fixed_dofs = []
        fixed_values = []

        # Identificar grados de libertad fijos por las condiciones de frontera
        for node in self.nodes:
            if hasattr(node, "boundary_label") and any("Dirichlet" in label for label in node.boundary_label):
                #print(f"🔒 Nodo {node.id} tiene condición de Dirichlet: {node.u}")
                fixed_dofs.append(node.id - 1)  # Convertir a índice base-0
                fixed_values.append(node.u)

        fixed_dofs = np.array(fixed_dofs, dtype=int)
        fixed_values = np.array(fixed_values, dtype=float)
        all_dofs = np.arange(self.nnodes)
        free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

        # Resolver el sistema de ecuaciones
        K = self.K_global
        f = self.f

        # Reducir el sistema para los grados de libertad libres
        f_reduced = f[free_dofs] - K[free_dofs][:, fixed_dofs] @ fixed_values
        u_free = spsolve(K[free_dofs][:, free_dofs], f_reduced)

        # Reconstituir la solución completa
        u_full = np.zeros(self.nnodes)
        u_full[fixed_dofs] = fixed_values
        u_full[free_dofs] = u_free

        self.u = u_full

        # Asignar la solución a los nodos
        for node in self.nodes:
            node.u_fem = u_full[node.id - 1]

    def real_solution(self):
        # Calcular la solución analítica
        for node in self.nodes:
            node.solve_u(self.alpha)

    def semi_norm_H1_0(self, alpha, orden=5):
        """
        Calcula la semi-norma H1_0 en el espacio de Sobolev.
        """
        puntos, pesos = np.polynomial.legendre.leggauss(orden)
        puntos = 0.5 * (puntos + 1)  # Ajustar para el dominio [0,1]
        pesos = 0.5 * pesos

        total = 0.0

        for i in range(orden):
            for j in range(orden):
                x = puntos[i]
                y = puntos[j]
                w = pesos[i] * pesos[j]

                r2 = x**2 + y**2
                if r2 == 0 and alpha < 1:
                    grad2 = 0.0  # Evitar singularidad
                else:
                    grad2 = alpha**2 * r2**(alpha - 1)

                total += grad2 * w

        return total  # Ya es la semi-norma al cuadrado

    def femm_solution(self):
        # Calcular la solución FEM
        x = np.zeros(self.nnodes)
        for node in self.nodes:
            if hasattr(node, "boundary_label") and any("Dirichlet" in label for label in node.boundary_label):
                x[node.id - 1] = node.u  # Usar solución de frontera
            else:
                x[node.id - 1] = node.u_fem  # Usar solución FEM

        return x @ self.K_global @ x
