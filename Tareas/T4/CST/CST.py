import numpy as np

class CST:
    def __init__(self, id, node_ids):
        self.id = id
        self.node_ids = node_ids  # Lista de 3 enteros
        self.K = None  # Se calcula luego con get_stiffness_matrix

    def get_B_matrix(self, nodes):
        """
        Calcula la matriz B del elemento CST. Si el área es demasiado pequeña,
        retorna una matriz nula para evitar NaNs.
        """
        n1, n2, n3 = [nodes[i - 1] for i in self.node_ids]

        x1, y1 = n1.x, n1.y
        x2, y2 = n2.x, n2.y
        x3, y3 = n3.x, n3.y

        # Área del triángulo (puede ser negativa si nodos están en orden horario)
        area = 0.5 * abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))

        B = np.array([
            [y2 - y3, y3 - y1, y1 - y2],
            [x3 - x2, x1 - x3, x2 - x1]
        ]) / (2 * area)


        return B, area

    def get_stiffness_matrix(self, nodes):

        B, area = self.get_B_matrix(nodes)

        if area == 0.0:
            self.K = np.zeros((3, 3))
        else:
            I = np.identity(2)
            self.K = area * B.T @ I @ B  # Matriz 3x3

        return self.K