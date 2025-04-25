import numpy as np

class Section:
    def __init__(self, thickness, Ex, Ey, nuxy, Gxy):
        self.thickness = thickness
        self.Ex = Ex
        self.Ey = Ey
        self.nuxy = nuxy
        self.Gxy = Gxy
        self.D = self._compute_D_orthotropic()

    def _compute_D_orthotropic(self):
        Ex, Ey, nuxy, Gxy = self.Ex, self.Ey, self.nuxy, self.Gxy
        nuyx = nuxy * Ey / Ex
        denom = 1 - nuxy * nuyx
        return np.array([
            [Ex / denom, nuyx * Ey / denom, 0],
            [nuxy * Ex / denom, Ey / denom, 0],
            [0, 0, Gxy]
        ])