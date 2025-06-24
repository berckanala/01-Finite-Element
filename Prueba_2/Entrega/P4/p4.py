import numpy as np

def shape_functions(xi, eta):
    """
    Shape functions for the NiceQuad (quadrilateral with four nodes) using triangular coordinates.
    The nodes are defined as 0, 1, and 3 for the corner nodes.
    """
    N0 = (1 - xi - eta)  # Corner node 0
    N1 = xi               # Corner node 1
    N3 = eta              # Corner node 3
    return np.array([N0, N1, N3])

def shape_function_derivatives(xi, eta):
    """
    Derivatives of the shape functions with respect to xi and eta for a NiceQuad
    """
    dN0_dxi = -1
    dN0_deta = -1
    dN1_dxi = 1
    dN1_deta = 0
    dN3_dxi = 0
    dN3_deta = 1

    return np.array([[dN0_dxi, dN1_dxi, dN3_dxi],
                     [dN0_deta, dN1_deta, dN3_deta]])

def jacobian(xi, eta, coords):
    """
    Calculate the Jacobian matrix for a quadrilateral element using triangular coordinates
    """
    # Extracting the coordinates of the nodes (x, y) for nodes 0, 1, 3
    x0, y0 = coords[0]
    x1, y1 = coords[1]
    x3, y3 = coords[2]

    # Shape function derivatives
    dN_dxi_eta = shape_function_derivatives(xi, eta)
    
    # Jacobian matrix calculation
    J = np.array([
        [np.dot(dN_dxi_eta[0], [x0, x1, x3]), np.dot(dN_dxi_eta[0], [y0, y1, y3])],
        [np.dot(dN_dxi_eta[1], [x0, x1, x3]), np.dot(dN_dxi_eta[1], [y0, y1, y3])]
    ])
    
    return J

def jacobian_determinant(J):
    """
    Compute the determinant of the Jacobian matrix
    """
    return np.linalg.det(J)

def main():
    # Coordinates of the quadrilateral element (NiceQuad) using triangular coordinates
    coords = np.array([
        [0.0, 0.0],  # Node 0
        [1.0, 0.0],  # Node 1
        [0.5, 0.5],  # Node 3
        [0.0, 1.0]   # Mid-node (this is the edge midpoint between nodes 0 and 1)
    ])
    
    # Generate test points using linspace for xi and eta
    xi_values = np.linspace(0, 1, 10)  # Generate 10 values from 0 to 1 for xi
    eta_values = np.linspace(0, 1, 10)  # Generate 10 values from 0 to 1 for eta
    
    naughty_points = []  # List to store points where det(J) <= 0

    # Loop through all combinations of xi and eta
    for xi in xi_values:
        for eta in eta_values:
            J = jacobian(xi, eta, coords)
            det_J = jacobian_determinant(J)
            if det_J <= 0:
                naughty_points.append((xi, eta, det_J))  # Store points where det(J) <= 0

    # Print results
    if naughty_points:
        print("Points where the Jacobian determinant is <= 0:")
        for point in naughty_points:
            print(f"xi = {point[0]:.2f}, eta = {point[1]:.2f}, det(J) = {point[2]:.4f}")
    else:
        print("All points have a positive Jacobian determinant (element is nice).")

if __name__ == "__main__":
    main()
