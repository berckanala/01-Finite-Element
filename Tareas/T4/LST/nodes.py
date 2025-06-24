# nodes.py
class Node:
    def __init__(self, id, x, y):
        self.id = id
        self.x = x
        self.y = y
        self.u = 0.0  # Inicializar u en 0.0
        self.u_fem = 0.0
        self.boundary_label = None

    def solve_u (self, alpha):
        # Método para calcular u basado en la posición del nodo
        self.u = (self.x**2 + self.y**2)**(alpha / 2)