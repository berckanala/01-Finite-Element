class Node:
    def __init__(self, id, x, y):
        self.id = id
        self.x = x
        self.y = y
        self.u = 0.0
        self.u_fem = 0.0
        self.boundary_label = []

    def solve_u(self, alpha):
        self.u = (self.x**2 + self.y**2)**(alpha / 2)
