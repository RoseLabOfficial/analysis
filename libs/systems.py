from . import *

class WholeCellRecording:
    def __init__(self, data: pd.DataFrame, parameters: pd.DataFrame):
        self.data = data
        self.parameters = parameters
        pass

    def scale_data(self):
        for idx, inj in enumerate(self.parameters["Iinj"]):
            self.data[inj] = ((self.data[inj] - self.data[inj][0])*1e-3)+self.parameters["Ess"][idx]
        return self.data
    
    def compute_polarizations(self, clamp_idx: int):
        if clamp_idx == None:
            min_clamp = self.parameters["Iinj"].min()
            


