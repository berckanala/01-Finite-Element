import numpy as np

class Material:
    def __init__(self, E, nu, gamma, type='planeStress'):
    
        self.E = E  # Young's modulus
        self.nu = nu  # Poisson's ratio
        self.type = type
        self.gamma = gamma
        self.Emat = self._compute_Emat()

    def _compute_Emat(self):

        E, nu = self.E, self.nu
        if isinstance(E, np.ndarray):
            return E

        if self.type == 'planeStress':
            return (E / (1 - nu**2)) * np.array([
                [1, nu, 0],
                [nu, 1, 0],
                [0, 0, (1 - nu) / 2]
            ])
        elif self.type == 'planeStrain':
            coef = E / ((1 + nu)*(1 - 2*nu))
            return coef * np.array([
                [1 - nu, nu, 0],
                [nu, 1 - nu, 0],
                [0, 0, (1 - 2*nu) / 2]
            ])
        else:
            raise ValueError(f"Invalid type: {self.type}")
        
    def get_properties(self):
        return {
            'E': self.E,
            'nu': self.nu,
            'gamma': self.gamma,
            'type': self.type
        }
   