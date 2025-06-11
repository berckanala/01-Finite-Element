import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import matrix_rank

class Pro_Solver:
    def __init__(self, nodes, elements, alpha, verbose=False):
        self.nodes = nodes
        self.elements = elements  # LST o CST
        self.alpha = alpha
        self.nnodes = len(nodes)
        self.verbose = verbose

        self._apply_dirichlet_conditions()
        self.compute_local_stiffness_matrices()
        self.K_global = self.assemble_global_matrix()
        self.f = self.assemble_force_vector()

        self.u = np.zeros(self.nnodes)
        from collections import defaultdict
        node_uses = defaultdict(int)
        for elem in self.elements:
            for nid in elem.node_ids:
                node_uses[nid] += 1
        unused = [nid for nid in range(1, self.nnodes + 1) if node_uses[nid] == 0]
        if unused:
            print("⚠️ Nodos no conectados a ningún elemento:", unused)
        for i, elem in enumerate(self.elements):
            if len(set(elem.node_ids)) < len(elem.node_ids):
                print(f"⚠️ Elemento {i+1} tiene nodos repetidos:", elem.node_ids)



    def _apply_dirichlet_conditions(self):
        for node in self.nodes:
            if hasattr(node, "boundary_label") and any("Dirichlet" in label for label in node.boundary_label):
                node.solve_u(self.alpha)

    def compute_local_stiffness_matrices(self):
        for i, element in enumerate(self.elements):
            element.get_stiffness_matrix(self.nodes)
            if element.detJ < 0:
                print(f"⚠️ Elemento {i+1} tiene Jacobiano singular o negativo (det = {element.detJ:.2e})")
                element.K = np.zeros((6, 6))  # Evita que lo use
            elif not np.any(element.



    def assemble_global_matrix(self):
        row_idx = []
        col_idx = []
        data = []

        for elem in self.elements:
            node_ids = [int(i) for i in elem.node_ids]
            K_local = elem.K
            for i_local in range(len(node_ids)):
                for j_local in range(len(node_ids)):
                    global_i = node_ids[i_local] - 1
                    global_j = node_ids[j_local] - 1
                    row_idx.append(global_i)
                    col_idx.append(global_j)
                    data.append(K_local[i_local, j_local])

        return coo_matrix((data, (row_idx, col_idx)), shape=(self.nnodes, self.nnodes)).tocsr()

    def assemble_force_vector(self):
        f_global = np.zeros(self.nnodes)

        for element in self.elements:
            f_local = element.get_load_vector(self.nodes, self.alpha)
            for i_local, node_id in enumerate(element.node_ids):
                f_global[node_id - 1] += f_local[i_local]

        return f_global

    def solve_system(self):
        fixed_dofs = []
        fixed_values = []

        for node in self.nodes:
            if hasattr(node, "boundary_label") and any("Dirichlet" in label for label in node.boundary_label):
                fixed_dofs.append(node.id - 1)
                fixed_values.append(node.u)

        fixed_dofs = np.array(fixed_dofs, dtype=np.intp)
        fixed_values = np.array(fixed_values)

        if len(fixed_dofs) == 0:
            raise ValueError("❌ No hay nodos con condición de Dirichlet. El sistema es indeterminado.")

        all_dofs = np.arange(self.nnodes)
        free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

        if self.verbose:
            print(f"🧱 Total Dirichlet nodes: {len(fixed_dofs)}")
            print(f"🌀 Nodos libres: {len(free_dofs)} → {free_dofs.tolist()}")

        K = self.K_global
        f = self.f.copy()

        f_reduced = f[free_dofs] - K[free_dofs][:, fixed_dofs] @ fixed_values
        K_free = K[free_dofs][:, free_dofs].todense()

        # Diagnóstico
        print("🔍 Tamaño sistema libre:", K_free.shape)
        print("📐 Rango matriz libre:", matrix_rank(K_free))
        print("🔁 Matriz simétrica:", np.allclose(K_free, K_free.T, atol=1e-10))

        try:
            u_free = spsolve(K[free_dofs][:, free_dofs], f_reduced)
        except Exception as e:
            raise RuntimeError(f"⚠️ Error al resolver el sistema: {e}")

        u_full = np.zeros(self.nnodes)
        u_full[fixed_dofs] = fixed_values
        u_full[free_dofs] = u_free

        self.u = u_full

        for node in self.nodes:
            node.u_fem = u_full[node.id - 1]

        if self.verbose:
            print("✅ Sistema resuelto con", len(free_dofs), "grados de libertad libres.")

    def compute_H1_seminorm_exact(self, order=3):
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
                grad2 = 0.0 if r2 == 0 and self.alpha < 1 else self.alpha**2 * r2**(self.alpha - 1)
                total += grad2 * w

        return total

    def compute_fem_energy(self):
        x = np.zeros(self.nnodes)
        for node in self.nodes:
            x[node.id - 1] = node.u if hasattr(node, "boundary_label") and any("Dirichlet" in label for label in node.boundary_label) else node.u_fem

        return x.T @ self.K_global @ x

    def export_nodal_solution(self):
        return np.array([[node.x, node.y, node.u_fem] for node in self.nodes])
