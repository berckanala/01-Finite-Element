import numpy as np
import sympy as sp

class Element:
    def __init__(self, id, node_ids, n_nodes):
        self.id = id
        self.node_ids = node_ids  # Lista de nodos
        self.num_nodes = n_nodes  # Número de nodos en el elemento
        self.K = None  # Se calculará con get_stiffness_matrix
        


    def gausspoints (self):

        if self.id in [500, 501, 502, 503]:
            gauss_points = [
            (1/6, 1/6, 1/6),
            (2/3, 1/6, 1/6),
            (1/6, 2/3, 1/6)
        ]
            
        else:
            #Ojo, hay que modificar los puntos de gauss dependiendo la forma
            sqrt3 = np.sqrt(3 / 5)
            points = [-sqrt3, 0, sqrt3]
            weights = [5/9, 8/9, 5/9]
            gauss_points = [(xi, eta, wx * wy) for xi, wx in zip(points, weights)
                                                for eta, wy in zip(points, weights)]
            
        return gauss_points

    def shape_functions(self, xi, eta):
        if self.id in [500, 501, 502, 503]:
            # Definir funciones de forma para nice_quad (p.ej. P1 en triángulo con 4 nodos)
            N1 = 1 - xi - eta
            N2 = xi * (1 - 2 * eta)
            N3 = 4 * xi * eta
            N4 = eta * (1 - 2 * xi)
            return np.array([N1, N2, N3, N4])
        else:
            # Funciones estándar (Quad9 degenerado a Quad4)
            N1 = eta*xi/4 - eta/4 - xi/4 + 1/4
            N2 = -eta*xi/4 - eta/4 + xi/4 + 1/4
            N3 =  eta*xi/4 + eta/4 + xi/4 + 1/4
            N4 = -eta*xi/4 + eta/4 - xi/4 + 1/4
            return np.array([N1, N2, N3, N4])

        

    def shape_function_derivatives(self, xi, eta):
        # Derivadas de las funciones de forma respecto a xi y eta

        # Las derivadas simbólicas con sympy
        xi, eta = sp.symbols('xi eta')
        N_sym = self.shape_functions(xi, eta)

        # Derivadas respecto a xi
        dN_dxi = np.array([sp.diff(N_i, xi) for N_i in N_sym])
        
        # Derivadas respecto a eta
        dN_deta = np.array([sp.diff(N_i, eta) for N_i in N_sym])

        # Convertir de sympy a numpy para la evaluación numérica
        dN_dxi_func = np.array([sp.lambdify((xi, eta), d, "numpy") for d in dN_dxi])
        dN_deta_func = np.array([sp.lambdify((xi, eta), d, "numpy") for d in dN_deta])

        return dN_dxi_func, dN_deta_func

    def get_B_matrix(self, nodes, xi, eta):
        coords = np.array([[node.x, node.y] for node in self.node_ids], dtype=float)

        
        dN_dxi_func, dN_deta_func = self.shape_function_derivatives(xi, eta)

        J = np.zeros((2, 2))
        for i in range(self.num_nodes):
            J[0, 0] += dN_dxi_func[i](xi, eta) * coords[i, 0]
            J[0, 1] += dN_dxi_func[i](xi, eta) * coords[i, 1]
            J[1, 0] += dN_deta_func[i](xi, eta) * coords[i, 0]
            J[1, 1] += dN_deta_func[i](xi, eta) * coords[i, 1]

        detJ = np.linalg.det(J)
        if abs(detJ) < 1e-12:
            return np.zeros((2, self.num_nodes)), 0.0

        Jinv = np.linalg.inv(J)

        dN_dx = Jinv[0, 0] * np.array([dN_dxi_func[i](xi, eta) for i in range(self.num_nodes)]) + \
                Jinv[0, 1] * np.array([dN_deta_func[i](xi, eta) for i in range(self.num_nodes)])

        dN_dy = Jinv[1, 0] * np.array([dN_dxi_func[i](xi, eta) for i in range(self.num_nodes)]) + \
                Jinv[1, 1] * np.array([dN_deta_func[i](xi, eta) for i in range(self.num_nodes)])

        B = np.vstack([dN_dx, dN_dy])  # 2 x num_nodes
        
        return B, detJ


    def get_stiffness_matrix(self, nodes):
        # Puntos de Gauss

        gauss_points = self.gausspoints()
        

        # Matriz de rigidez global
        K = np.zeros((self.num_nodes, self.num_nodes))


        # Integración numérica
        for xi, eta, w in gauss_points:
            B, detJ = self.get_B_matrix(nodes, xi, eta)
            if detJ == 0:
                continue
            
            K += w * detJ * (B.T @ B)


        self.K = K
        return self.K
    
    def get_load_vector(self, nodes, alpha):
        coords = np.array([[node.x, node.y] for node in self.node_ids], dtype=float)


        gauss_points = self.gausspoints()


        f_local = np.zeros(self.num_nodes)  # Vector de carga local

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