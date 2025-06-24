import numpy as np
import matplotlib.pyplot as plt
globalParameters = {'nDoF': 3} 
class Node:
    def __init__(self, name, coordenadas, nodalLoad=None, restrain=None, printSummary=True):
        global globalParameters  # Parámetros globales del sistema
        
        self.name = name  # Identificación del nodo
        self.coordenadas = np.array(coordenadas)  # Coordenadas (x, y)
        
        # Calcular índices en el sistema global
        self.idx = self._indices()
        
        # Asignar carga nodal
        if nodalLoad is not None:
            if len(nodalLoad) == globalParameters['nDoF']:
                self.nodalLoad = np.array(nodalLoad)
            else:
                raise ValueError(f'La carga nodal debe tener {globalParameters["nDoF"]} valores.')
        else:
            self.nodalLoad = np.zeros(globalParameters['nDoF'])
        
        # Asignar restricciones
        if restrain is not None:
            if isinstance(restrain, list) and len(restrain) == globalParameters['nDoF']:
                self.restrain = np.array(restrain)
            else:
                raise ValueError(f'Las restricciones deben ser una lista con {globalParameters["nDoF"]} valores.')
        else:
            self.restrain = np.full(globalParameters['nDoF'], 'f')  # Por defecto, sin restricciones
        
        if printSummary:
            self.printSummary()
    
    def __str__(self):
        return f"Node {self.name} at {self.coordenadas}"
        
    def _indices(self):
        """
        Calcula los índices de los grados de libertad en el sistema global.
        """
        base_idx = (self.name - 1) * globalParameters['nDoF']  # Ajuste por indexación 0 en Python
        return np.array([base_idx, base_idx + 1, base_idx + 2])
    
    def plotGeometry(self, ax=None, text=False):
        """
        Dibuja la ubicación del nodo en un gráfico.
        """
        if ax is None:
            fig, ax = plt.subplots()
        
        ax.plot(self.coordenadas[0], self.coordenadas[1], 'ro')  # Dibujar nodo
        if text:
            ax.text(self.coordenadas[0], self.coordenadas[1], f'{self.name}', fontsize=12)
        
        return ax
    
    def printSummary(self):
        """
        Imprime un resumen con la información del nodo.
        """
        print(f'--------------------------------------------')
        print(f"Node {self.name} at {self.coordenadas}")
        print(f"Indices en matriz global: {self.idx}")
        print(f"Carga nodal: {self.nodalLoad}")
        print(f"Restricciones: {self.restrain}")
        print(f'--------------------------------------------\n')