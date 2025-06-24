import numpy as np

def shape_functions(xi, eta):
    """
    Shape functions for the NiceQuad (quadrilateral with four nodes) using triangular coordinates.
    The nodes are defined as 0, 1, and 3 for the corner nodes.
    """
    # For the triangular part:
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
    # Extracting the coordinates of the nodes (x, y) for nodes 0, 1, 3, and mid-node 0-1
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

def confirmation(shape):
    j=0
    for i in shape:
        j+=i
    if j!= 1:
        print(f"La suma de las funciones de interpolación no suman 1: {j}")
    else:
        print(f"La suma de las funciones de interpolación suman {j} correctamente.")



def main():
    # Coordinates of the quadrilateral element (NiceQuad) using triangular coordinates
    # Updated coordinates (x, y) for nodes 0, 1, 3, and mid-node 0-1
    coords = np.array([
        [0.0, 0.0],  # Node 0
        [1.0, 0.0],  # Node 1
        [0.5, 0.5],  # Node 3
        [0.0, 1.0]   # Mid-node (this is the edge midpoint between nodes 0 and 1)
    ])
    
    # Natural coordinates for the center of the element
    xi = 0.7
    eta = 0.3
    
    # Confirmation of shape functions
    shape = shape_functions(xi, eta)
    confirmation(shape)

    # Calculate Jacobian matrix
    J = jacobian(xi, eta, coords)
    
    # Calculate the determinant of the Jacobian
    det_J = jacobian_determinant(J)
    
    print(f"Jacobian Matrix at (xi={xi}, eta={eta}):\n{J}")
    print(f"Determinant of the Jacobian: {det_J}")

if __name__ == "__main__":
    main()
