import numpy as np
import sympy as sp
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

# Define the Node class with u_fem and boundary labels
class Node:
    def __init__(self, id, x, y):
        self.id = id
        self.x = x
        self.y = y
        self.u = 0.0  # Initial displacement (analytical)
        self.u_fem = 0.0  # Displacement (FEM solution)
        self.boundary_label = []

    def solve_u(self, alpha):
        # Calculate u based on the node's position (analytical solution)
        self.u = (self.x**2 + self.y**2)**(alpha / 2)

# Known manufactured solution g(x, y) = x * y * (1 - x) * (1 - y)
def g(x, y):
    return x * y * (1 - x) * (1 - y)

# Compute the source term f(x, y) using the Poisson equation
def f(x, y):
    d2g_dx2 = 2 * (y - 1) * (y) * (1 - 2 * x)
    d2g_dy2 = 2 * (x - 1) * (x) * (1 - 2 * y)
    return -(d2g_dx2 + d2g_dy2)

# Shape functions for a Quad7 element (7 nodes)
def shape_functions(xi, eta):
    N1 = eta*xi/4 - eta/4 - xi/4 + 1/4
    N2 = -eta*xi/4 - eta/4 + xi/4 + 1/4
    N3 = eta*xi/4 + eta/4 + xi/4 + 1/4
    N4 = -eta*xi/4 + eta/4 - xi/4 + 1/4
    return np.array([N1, N2, N3, N4])

# Derivatives of shape functions with respect to xi and eta
def shape_function_derivatives(xi, eta):
    xi, eta = sp.symbols('xi eta')
    N_sym = shape_functions(xi, eta)

    # Derivatives with respect to xi and eta
    dN_dxi = np.array([sp.diff(N_i, xi) for N_i in N_sym])
    dN_deta = np.array([sp.diff(N_i, eta) for N_i in N_sym])

    # Convert from sympy to numpy for numerical evaluation
    dN_dxi_func = np.array([sp.lambdify((xi, eta), d, "numpy") for d in dN_dxi])
    dN_deta_func = np.array([sp.lambdify((xi, eta), d, "numpy") for d in dN_deta])

    return dN_dxi_func, dN_deta_func

# Gauss points for numerical integration
def gausspoints():
    sqrt3 = np.sqrt(3 / 5)
    points = [-sqrt3, 0, sqrt3]
    weights = [5/9, 8/9, 5/9]
    gauss_points = [(xi, eta, wx * wy) for xi, wx in zip(points, weights)
                                            for eta, wy in zip(points, weights)]
    return gauss_points

# Get B matrix for an element at Gauss points
def get_B_matrix(node_ids, nodes, xi, eta):
    coords = np.array([[nodes[i - 1].x, nodes[i - 1].y] for i in node_ids])
    dN_dxi_func, dN_deta_func = shape_function_derivatives(xi, eta)

    # Jacobian matrix calculation
    J = np.zeros((2, 2))
    for i in range(4):  # Assuming 4 nodes for simplicity
        J[0, 0] += dN_dxi_func[i](xi, eta) * coords[i, 0]
        J[0, 1] += dN_dxi_func[i](xi, eta) * coords[i, 1]
        J[1, 0] += dN_deta_func[i](xi, eta) * coords[i, 0]
        J[1, 1] += dN_deta_func[i](xi, eta) * coords[i, 1]

    detJ = np.linalg.det(J)
    if abs(detJ) < 1e-12:
        return np.zeros((2, 4)), 0.0

    Jinv = np.linalg.inv(J)

    dN_dx = Jinv[0, 0] * np.array([dN_dxi_func[i](xi, eta) for i in range(4)]) + \
            Jinv[0, 1] * np.array([dN_deta_func[i](xi, eta) for i in range(4)])

    dN_dy = Jinv[1, 0] * np.array([dN_dxi_func[i](xi, eta) for i in range(4)]) + \
            Jinv[1, 1] * np.array([dN_deta_func[i](xi, eta) for i in range(4)])

    B = np.vstack([dN_dx, dN_dy])
    return B, detJ

# Assemble global stiffness matrix for the mesh
def assemble_global_matrix(nodes, elements):
    row_idx = []
    col_idx = []
    data = []

    for element in elements:
        K_local = np.zeros((4, 4))  # Assuming 4 nodes per element for simplicity
        node_ids = element['node_ids']

        # Assemble the local stiffness matrix
        gauss_points = gausspoints()
        for xi, eta, w in gauss_points:
            B, detJ = get_B_matrix(node_ids, nodes, xi, eta)
            if detJ == 0:
                continue
            K_local += w * detJ * (B.T @ B)

        for i_local in range(len(node_ids)):
            for j_local in range(len(node_ids)):
                global_i = node_ids[i_local] - 1
                global_j = node_ids[j_local] - 1
                row_idx.append(global_i)
                col_idx.append(global_j)
                data.append(K_local[i_local, j_local])

    # Return the global stiffness matrix in CSR format for efficient solving
    return coo_matrix((data, (row_idx, col_idx)), shape=(len(nodes), len(nodes))).tocsr()

# Assemble global load vector for the mesh
def assemble_global_load_vector(nodes, elements, alpha):
    f_global = np.zeros(len(nodes))

    for element in elements:
        f_local = np.zeros(4)  # Assuming 4 nodes per element
        node_ids = element['node_ids']
        gauss_points = gausspoints()

        for xi, eta, w in gauss_points:
            N = shape_functions(xi, eta)
            x = np.dot(N, [nodes[i - 1].x for i in node_ids])
            y = np.dot(N, [nodes[i - 1].y for i in node_ids])
            r2 = x**2 + y**2

            f_val = -alpha**2 * r2**(alpha / 2 - 1) if r2 != 0 else 0.0
            _, detJ = get_B_matrix(node_ids, nodes, xi, eta)
            f_local += f_val * w * detJ * N

        for i_local, node_id in enumerate(node_ids):
            global_i = node_id - 1
            f_global[global_i] += f_local[i_local]

    return f_global

# Solve the FEM system
def solve_fem(nodes, elements, alpha):
    K_global = assemble_global_matrix(nodes, elements)
    f_global = assemble_global_load_vector(nodes, elements, alpha)

    # Apply boundary conditions (Dirichlet)
    fixed_dofs = [0]  # Example: Fix node 1 (index 0)
    fixed_values = [nodes[0][2]]  # Example: Set the displacement for node 1

    # Apply the boundary conditions
    for i in fixed_dofs:
        K_global[i, :] = 0
        K_global[:, i] = 0
        K_global[i, i] = 1
        f_global[i] = fixed_values[fixed_dofs.index(i)]

    # Solve the system
    u = spsolve(K_global, f_global)

    # Update the nodes with the FEM solution
    for i, node in enumerate(nodes):
        node.u_fem = u[i]

    return u

# Visualization function
def plot_solution_3d(nodes):
    x_coords = [node[0] for node in nodes]
    y_coords = [node[1] for node in nodes]
    z_fem = [node.u_fem for node in nodes]  # FEM solution
    z_analytical = [node.u for node in nodes]  # Analytical solution

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x_coords, y_coords, z_fem, c='b', marker='o', label='FEM Solution')
    ax.scatter(x_coords, y_coords, z_analytical, c='r', marker='^', label='Analytical Solution')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Solution (u)')
    ax.legend()
    ax.set_title('Comparison of FEM and Analytical Solutions')
    plt.show()

# Main function to run the problem
def main(alpha):
    nodes = [
        Node(1, 0.0, 0.0),  # Node 1
        Node(2, 1.0, 0.0),  # Node 2
        Node(3, 0.5, 0.5),  # Node 3 (analytical solution)
        Node(4, 0.0, 1.0),  # Node 4
        Node(5, 1.0, 1.0)   # Node 5 (analytical solution)
    ]

    elements = [
        {'node_ids': [1, 2, 3, 4]},
        {'node_ids': [2, 5, 4, 3]}
    ]

    # Solve the FEM system
    u_fem = solve_fem(nodes, elements, alpha)

    # Plot the results
    plot_solution_3d(nodes)

if __name__ == "__main__":
    alpha = 3  # Change as needed
    main(alpha)
