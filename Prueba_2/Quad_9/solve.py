import numpy as np

class Solve:
    def __init__(self, nodes, elements):
        self.nodes = nodes
        self.elements = elements
        self.ndof = max(dof for node in nodes for dof in node.dofs)  # máximo DOF usado
        self.K_global = np.zeros((self.ndof + 1, self.ndof + 1))
        self.f_global = np.zeros((self.ndof + 1, 1))
        self.u_global = np.zeros((self.ndof + 1, 1))

    def assemble(self):
        #print("🔧 Ensamblando matriz global...")

        for elem in self.elements:
            ke = elem.Kg
            idx = elem.calculate_indices()

            if ke.shape != (len(idx), len(idx)):
                raise ValueError(f"❌ Dimensión inconsistente: Ke {ke.shape}, idx {len(idx)}")

            for i in range(len(idx)):
                for j in range(len(idx)):
                    if idx[i] >= self.K_global.shape[0] or idx[j] >= self.K_global.shape[1]:
                        #print(f"❌ DOF fuera de rango: idx[{i}] = {idx[i]}, idx[{j}] = {idx[j]}")
                        continue
                    self.K_global[idx[i], idx[j]] += ke[i, j]


    def apply_force(self, dof_index, value):
        self.f_global[dof_index] += value

    def apply_forces_vector(self, force_vector):
        self.f_global += force_vector.reshape(-1, 1)

    def apply_boundary_conditions(self):
        self.fixed_dofs = []
        self.free_dofs = []

        for node in self.nodes:
            for dof_val, dof_idx in zip(node.restrain, node.dofs):
                if dof_val == 1:
                    self.fixed_dofs.append(dof_idx)
                else:
                    self.free_dofs.append(dof_idx)

        self.fixed_dofs = np.array(self.fixed_dofs)
        self.free_dofs = np.array(self.free_dofs)

        for dof in self.fixed_dofs:
            self.K_global[dof, :] = 0
            self.K_global[:, dof] = 0
            self.K_global[dof, dof] = 1
            self.f_global[dof] = 0

    def check_zero_rows(self):
        zero_rows = np.where(~self.K_global.any(axis=1))[0]
        return zero_rows

    def solve(self):
        self.assemble()

        self.K_original = self.K_global.copy()
        self.f_original = self.f_global.copy()

        self.apply_boundary_conditions()

        used_dofs = sorted(set(dof for node in self.nodes for dof in node.dofs))
        #print(f"DOFs used: {used_dofs}")
        K_reduced = self.K_global[np.ix_(used_dofs, used_dofs)]
        #print(K_reduced)
        f_reduced = self.f_global[used_dofs]

        rowsums = np.sum(np.abs(K_reduced), axis=1)
        zero_rows = np.where(rowsums == 0)[0]
        if len(zero_rows) > 0:
            #print(f"❌ Filas completamente nulas en K_reduced: {zero_rows}")
            raise ValueError("Sistema subdeterminado: nodos con DOF sin rigidez.")


        u_reduced = np.linalg.solve(K_reduced, f_reduced)

        self.u_global = np.zeros_like(self.f_global)
        self.u_global[used_dofs] = u_reduced
        return self.u_global

    def get_displacement_at_node(self, node_index):
        node = self.nodes[node_index - 1]  # cuidado: node.index parte en 1
        ux = self.u_global[node.dofs[0], 0]
        uy = self.u_global[node.dofs[1], 0]
        return ux, uy

    def compute_reactions(self):
        R_total = self.K_original @ self.u_global - self.f_original

        self.reactions = np.zeros_like(self.f_global)
        self.reactions[self.fixed_dofs] = R_total[self.fixed_dofs]

        return self.reactions

    def print_applied_forces(self):
        f = self.f_original if hasattr(self, 'f_original') else self.f_global
        for node in self.nodes:
            dof_x, dof_y = node.dofs
            fx = f[dof_x][0] if dof_x < len(f) else 0.0
            fy = f[dof_y][0] if dof_y < len(f) else 0.0
            #print(f"Nodo {node.index}: Fx = {fx:.2f}, Fy = {fy:.2f}")

    def print_summary(self):
        #print("Desplazamientos por nodo:")
        for node in self.nodes:
            ux, uy = self.get_displacement_at_node(node.index)
            #print(f"Nodo {node.index}: ux = {ux:.6e}, uy = {uy:.6e}")