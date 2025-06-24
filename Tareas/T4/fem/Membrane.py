class Membrane:
    def __init__(self,
                 name:str,
                 thickness:float, 
                 material:object):
        
        self.name = name
        self.thickness = thickness
        self.material = material
        
    def __str__(self):
        return f"Membrane: {self.name} - Mat:{self.material}"
    
    def get_Emat(self, type: str):
        return self.material.get_Emat(type)
        