import numpy as np
import matplotlib.pyplot as plt
from scipy.special import roots_legendre
import matplotlib.patches as patches
import sympy as sp
x, y = sp.symbols('x y')

class Quad9:
    def __init__(self, elementTag, node_list, membrane, type='planeStress', samplingPoints=3, load_direction=None, eval_points=[0, 0], nDof=2):
        if len(node_list) != 9:
            raise ValueError("node_list must contain exactly 9 nodes.")
        
        self.nDof = nDof
        self.node_list = node_list
        self.elementTag = elementTag
        self.thickness = membrane.thickness
        self.material = membrane.material
        self.type = type
        self.samplingPoints = samplingPoints
        self.load_direction = load_direction or [0, 0]
        self.eval_points = eval_points
        self.C = membrane.material.Emat

        self.xy = np.array([node.coord for node in self.node_list])
        self.index = np.hstack([node.dofs for node in self.node_list])
        self.Kg, self.A, self.F_fe_global = self.calculate_K0()
        self.B = None

    def calculate_interpolation_functions(self, zeta, eta):
        # Funciones de forma
        N = np.zeros((2, 18))
        dNnatural = np.zeros((2, 9))  # [fila 0: dN/dζ, fila 1: dN/dη]

        # Alias
        ξ = zeta
        η = eta

        # Funciones de forma (orden estándar)
        N1 = 0.25 * ξ * η * (ξ - 1) * (η - 1)
        N2 = 0.25 * ξ * η * (ξ + 1) * (η - 1)
        N3 = 0.25 * ξ * η * (ξ + 1) * (η + 1)
        N4 = 0.25 * ξ * η * (ξ - 1) * (η + 1)
        N5 = 0.5 * (1 - ξ ** 2) * η * (η - 1)
        N6 = 0.5 * ξ * (ξ + 1) * (1 - η ** 2)
        N7 = 0.5 * (1 - ξ ** 2) * η * (η + 1)
        N8 = 0.5 * ξ * (ξ - 1) * (1 - η ** 2)
        N9 = (1 - ξ ** 2) * (1 - η ** 2)

        Ns = [N1, N2, N3, N4, N5, N6, N7, N8, N9]

        for i in range(9):
            N[0, 2 * i] = Ns[i]
            N[1, 2 * i + 1] = Ns[i]

        # Derivadas exactas respecto a ζ (zeta)
        dNnatural[0, 0] = 0.25 * η * (2 * ξ - 1) * (η - 1)
        dNnatural[0, 1] = 0.25 * η * (2 * ξ + 1) * (η - 1)
        dNnatural[0, 2] = 0.25 * η * (2 * ξ + 1) * (η + 1)
        dNnatural[0, 3] = 0.25 * η * (2 * ξ - 1) * (η + 1)
        dNnatural[0, 4] = -ξ * η * (η - 1)
        dNnatural[0, 5] = 0.5 * (2 * ξ + 1) * (1 - η ** 2)
        dNnatural[0, 6] = -ξ * η * (η + 1)
        dNnatural[0, 7] = 0.5 * (2 * ξ - 1) * (1 - η ** 2)
        dNnatural[0, 8] = -2 * ξ * (1 - η ** 2)

        # Derivadas exactas respecto a η (eta)
        dNnatural[1, 0] = 0.25 * ξ * (ξ - 1) * (2 * η - 1)
        dNnatural[1, 1] = 0.25 * ξ * (ξ + 1) * (2 * η - 1)
        dNnatural[1, 2] = 0.25 * ξ * (ξ + 1) * (2 * η + 1)
        dNnatural[1, 3] = 0.25 * ξ * (ξ - 1) * (2 * η + 1)
        dNnatural[1, 4] = 0.5 * (1 - ξ ** 2) * (2 * η - 1)
        dNnatural[1, 5] = -η * ξ * (ξ + 1)
        dNnatural[1, 6] = 0.5 * (1 - ξ ** 2) * (2 * η + 1)
        dNnatural[1, 7] = -η * ξ * (ξ - 1)
        dNnatural[1, 8] = -2 * η * (1 - ξ ** 2)

        return N, dNnatural

    def calculate_B_matrix(self, zeta=0.0, eta=0.0):
        """
        Calcula la matriz B, Jacobiano (J), determinante de Jacobiano (J_det) y las funciones de forma (N)
        para un punto de integración dado (zeta, eta).
        Si no se pasan valores de zeta y eta, utiliza los valores predeterminados (0.0, 0.0).
        """
        N, dNnatural = self.calculate_interpolation_functions(zeta, eta)
        
        # Calcular el Jacobiano
        J = dNnatural @ self.xy
        J_det = np.linalg.det(J)

        if J_det <= 0:
            raise ValueError("Jacobiano no positivo.")

        # Calculamos las derivadas en coordenadas cartesianas
        dNcartesian = np.linalg.solve(J, dNnatural)
        
        # Crear la matriz B
        B = np.zeros((3, 18))  # 3 filas, 18 columnas (9 nodos * 2 grados de libertad)
        
        B[0, 0::2] = dNcartesian[0, :]  # Derivadas respecto a x (dN/dx)
        B[1, 1::2] = dNcartesian[1, :]  # Derivadas respecto a y (dN/dy)
        B[2, 0::2] = dNcartesian[1, :]  # Derivadas cruzadas (dN/dy respecto a x)
        B[2, 1::2] = dNcartesian[0, :]  # Derivadas cruzadas (dN/dx respecto a y)
        
        return B, J, J_det, N

    def calculate_K0(self):
        roots, weights = roots_legendre(self.samplingPoints)
        t = self.thickness
        b = np.array(self.load_direction).reshape(-1, 1)
        Ke = np.zeros((18, 18))
        fe = np.zeros((18, 1))
        A = 0

        for r, wr in zip(roots, weights):
            for s, ws in zip(roots, weights):
                B, _, J_det, N = self.calculate_B_matrix(r, s)
                Ke += wr * ws * t * B.T @ self.C @ B * J_det
                fe += wr * ws * N.T @ b * J_det
                A += wr * ws * abs(J_det)

        fe *= t * self.material.gamma
        return Ke, A, fe.flatten()

    def get_element_displacements(self, u_global):
        return u_global[self.index].flatten()

    def von_mises_stress(self, u_global):
        ue = self.get_element_displacements(u_global)
        B, _, _, _ = self.calculate_B_matrix(0, 0)
        self.B = B
        sigma = self.C @ (B @ ue)
        σx, σy, τxy = sigma
        return float(np.sqrt(σx**2 - σx * σy + σy**2 + 3 * τxy**2))

    def get_stress(self, u_global):
        # Verificar si B ya está calculada, si no, calcularla
        if self.B is None:
            # Calcular la matriz B si no está calculada
            self.B, J, J_det, N = self.calculate_B_matrix()  # Usamos el valor predeterminado de zeta y eta (0, 0)
        
        # Verificar que B tenga las dimensiones correctas
        if self.B.shape != (3, 18):
            raise ValueError(f"Dimensiones incorrectas para la matriz B: {self.B.shape}. Se esperaba (3, 18).")
        
        # Obtener los desplazamientos del elemento
        ue = self.get_element_displacements(u_global)
        
        # Calcular la tensión
        sigma = self.C @ (self.B @ ue)  # Ahora B está calculada y tiene las dimensiones correctas
        return sigma


    def get_strain(self, u_global):
        ue = self.get_element_displacements(u_global)
        return self.B @ ue

    def calculate_indices(self):
        return self.index

    def element_visualization(self, offset=0):
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches

        fig, ax = plt.subplots()
        ax.set_aspect('equal', adjustable='box')

        for n, node in enumerate(self.xy):
            ax.plot(node[0], node[1], 'ko', ms=6)
            label = f'{self.node_list[n].name}'
            ax.text(node[0] + offset, node[1] + offset, label, fontsize=10)

        polygon = patches.Polygon(xy=self.xy[:4], edgecolor='black', facecolor='grey', alpha=0.3)
        ax.add_patch(polygon)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'Quad9 Element {self.elementTag}')
        plt.grid(True)
        plt.show()

    def get_centroid(self):
        """
        Devuelve el centroide del elemento cuadrilateral de 9 nodos
        como el promedio de las coordenadas de todos los nodos.
        """
        coords = np.array([node.coord for node in self.node_list])
        return np.mean(coords, axis=0)
    
    def apply_point_body_force(self, x, y, force_vector):
        """
        Aplica una fuerza puntual en (x, y) interpolándola con las
        funciones de forma del Quad9.

        Returns:
            f_puntual (ndarray): Vector de fuerza equivalente (18 x 1)
        """
        # Asumimos aplicación centrada para peso propio → centro natural
        zeta = 0.0
        eta = 0.0

        N, _ = self.calculate_interpolation_functions(zeta, eta)  # (2 x 18)
        fx, fy = force_vector
        fuerza = np.array([[fx], [fy]])  # (2 x 1)

        f_puntual = (N.T @ fuerza).flatten()  # (18 x 1) → flatten
        return f_puntual

