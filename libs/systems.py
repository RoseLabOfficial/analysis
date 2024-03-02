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
    
    def filter_data(self, cutoff: float = 100.0, stopband_attenuation: float = 80, passband_ripple: float = 0.01):
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        filter = LowPassFilter(sampling_rate)
        for idx, inj in enumerate(self.parameters["Iinj"]):
            self.data["filtered_"+str(inj)] = filter.filter(self.data[inj], cutoff, stopband_attenuation, passband_ripple)
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
            self.data.at[0, "Im_"+str(clamp)] = 0.0
        return self.data
    
class LowPassFilter:
    def __init__(self, sampling_rate: float):
        self.sampling_rate = sampling_rate
        pass

    def compute_minimum_order(self, cutoff: float, stopband_attenuation: float, passband_ripple: float):
        passband_edge_frequency = cutoff/(self.sampling_rate/2)
        stopband_edge_frequency = 1.5*passband_edge_frequency
        order, normalized_cutoff_frequency = filters.buttord(passband_edge_frequency, stopband_edge_frequency, 
                                                            passband_ripple, stopband_attenuation)
        return order, normalized_cutoff_frequency
    
    def filter(self, input: float, cutoff: float, stopband_attenuation: float, passband_ripple: float):
        order, normalized_frequency = self.compute_minimum_order(cutoff, stopband_attenuation, passband_ripple)
        second_order_sections = filters.butter(order, normalized_frequency, output="sos")
        return filters.sosfilt(second_order_sections, input)