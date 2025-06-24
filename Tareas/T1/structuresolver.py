import numpy as np

class StructureSolver:
    def __init__(self, nodos, K_global, nDoF):
        self.nodos = nodos
        self.K_global = K_global
        self.nDoF = nDoF
        self.total_dof = len(nodos) * nDoF

        # Inicialización de vectores y matrices
        self.free_dofs = []
        self.fixed_dofs = []

        self._clasificar_grados_de_libertad()
        self._construir_vectores_de_carga()
        self._particionar_matrices()
        self._resolver_sistema()

    def _clasificar_grados_de_libertad(self):
        for nodo in self.nodos:
            for i, restriccion in enumerate(nodo.restrain):
                dof = nodo.idx[i]
                if restriccion == 'r':
                    self.fixed_dofs.append(dof)
                else:
                    self.free_dofs.append(dof)

        self.free_dofs = np.array(sorted(self.free_dofs))
        self.fixed_dofs = np.array(sorted(self.fixed_dofs))

    def _construir_vectores_de_carga(self):
        self.P = np.zeros(self.total_dof)

        # Sumar cargas nodales
        for nodo in self.nodos:
            self.P[nodo.idx] += nodo.nodalLoad

        self.P_f = self.P[self.free_dofs]
        self.P_c = self.P[self.fixed_dofs]
        self.U_c = np.zeros(len(self.fixed_dofs))  # Se asume desplazamiento nulo en restricciones

    def _particionar_matrices(self):
        K = self.K_global
        fd = self.free_dofs
        rd = self.fixed_dofs

        self.K_ff = K[np.ix_(fd, fd)]
        self.K_fc = K[np.ix_(fd, rd)]
        self.K_cf = K[np.ix_(rd, fd)]
        self.K_cc = K[np.ix_(rd, rd)]

    def _resolver_sistema(self):
        # Resolver desplazamientos libres
        self.U_f = np.linalg.solve(self.K_ff, self.P_f - self.K_fc @ self.U_c)

        # Vector completo de desplazamientos
        self.U = np.zeros(self.total_dof)
        self.U[self.free_dofs] = self.U_f
        self.U[self.fixed_dofs] = self.U_c

        # Calcular reacciones en los apoyos
        self.R = self.K_cf @ self.U_f + self.K_cc @ self.U_c - self.P_c

    def get_displacements(self):
        return self.U

    def get_reactions(self):
        return self.R

    def resumen(self):
        print("\n🔹 Desplazamientos [m, rad]:")
        for i, u in enumerate(self.U):
            print(f"  DoF {i}: {u:.6e}")

        print("\n🔹 Reacciones en Apoyos [N, Nm]:")
        for i, r in zip(self.fixed_dofs, self.R):
            print(f"  DoF {i}: {r:.6e}")
