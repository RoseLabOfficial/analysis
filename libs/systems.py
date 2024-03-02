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
    
    def compute_active_conductance_constants(self):
        self.parameters["alpha"] = (1.0/self.parameters["Rin"])/(self.parameters["Et"] - self.parameters["Eact"])
        self.parameters["beta"] = self.parameters["alpha"]*(self.parameters["Et"] - self.parameters["Er"])
        self.parameters["alpha"] = self.parameters["alpha"]*self.parameters["xalpha"]
        self.parameters["beta"] = self.parameters["beta"]*self.parameters["xbeta"]
        return self.parameters
    
    def compute_polarizations(self):
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            self.data["depolarization_"+str(clamp)] = np.where(self.data[clamp] > self.parameters["Ess"][idx], self.data[clamp] - self.parameters["Ess"][idx], 0)
            self.data["hyperpolarization_"+str(clamp)] = np.where(self.data[clamp] < self.parameters["Ess"][idx], self.data[clamp] - self.parameters["Ess"][idx], 0)
        return self.data
    
    def compute_leakage_currents(self):
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            self.data["Ileak_"+str(clamp)] = (1/self.parameters["Rin"][idx])*(self.data[clamp] - self.parameters["Er"][idx])
            print(idx)
        return self.data
    
    def compute_active_currents(self):
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            self.data["Iactive_"+str(clamp)] = self.parameters["alpha"][idx]*(self.data[clamp] - self.parameters["Er"][idx])*(self.data[clamp] - self.parameters["Et"][idx]) + \
                                            self.parameters["beta"][idx]*(self.data[clamp] - self.parameters["Er"][idx])
        return self.data
    
    def compute_membrane_currents(self):
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            self.data["Im_"+str(clamp)] = self.parameters["Cm"][idx]*self.data[clamp].diff()
            self.data["Im_"+str(clamp)][0] = 0.0
        return self.data
    
    