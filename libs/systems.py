from . import *

class WholeCellRecording:
    def __init__(self, data: pd.DataFrame, parameters: pd.DataFrame):
        self.data = data
        self.parameters = parameters
        pass

    def scale_data(self):
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            self.data[clamp] = ((self.data[clamp] - self.data[clamp][0])*1e-3)+self.parameters["Ess"][idx]
        return self.data
    
    def filter_data(self, cutoff: float = 200.0, stopband_attenuation: float = 80, passband_ripple: float = 0.01):
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        filter = LowPassFilter(sampling_rate)
        for idx, inj in enumerate(self.parameters["Iinj"]):
            self.data[inj] = filter.filter(self.data[inj], cutoff, stopband_attenuation, passband_ripple)
        return self.data
    
    def compute_active_conductance_constants(self):
        self.parameters["alpha"] = (1.0/self.parameters["Rin"])/(self.parameters["Eact"] - self.parameters["Ess"])
        self.parameters["beta"] = self.parameters["alpha"]*(self.parameters["Et"] - self.parameters["Ess"])
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
    
    def filter_membrane_currents(self, cutoff: float = 200.0, stopband_attenuation: float = 80, passband_ripple: float = 0.01):
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        filter = LowPassFilter(sampling_rate)
        for idx, inj in enumerate(self.parameters["Iinj"]):
            self.data["filtered_Im_"+str(inj)] = filter.filter(self.data["Im_"+str(inj)], cutoff, stopband_attenuation, passband_ripple)
        return self.data
    
    def compute_passive_conductances(self):
        ntimesteps = self.data.shape[0]
        A = np.zeros((ntimesteps, 2, 2))
        B = np.zeros((ntimesteps, 2, 1))
        self.data["excitation"] = self.data["times"]*0.0
        self.data["inhibition"] = self.data["times"]*0.0
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            A[:, 0, 0] = A[:, 0, 0] + np.square(self.data[clamp] - self.parameters["Ee"][idx])
            A[:, 0, 1] = A[:, 0, 1] + (self.data[clamp] - self.parameters["Ee"][idx]) * (self.data[clamp] - self.parameters["Ei"][idx])
            A[:, 1, 1] = A[:, 1, 1] + np.square(self.data[clamp] - self.parameters["Ei"][idx])
            B[:, 0, 0] = B[:, 0, 0] + (self.data["filtered_Im_"+str(clamp)] - self.data["Iactive_"+str(clamp)] - clamp + self.data["Ileak_"+str(clamp)])*(self.data[clamp] - self.parameters["Ee"][idx])
            B[:, 1, 0] = B[:, 1, 0] + (self.data["filtered_Im_"+str(clamp)] - self.data["Iactive_"+str(clamp)] - clamp + self.data["Ileak_"+str(clamp)])*(self.data[clamp] - self.parameters["Ei"][idx])
        A[:, 1, 0] = A[:, 0, 1]
        for tdx in range(ntimesteps):
            G = np.linalg.inv(A[tdx, ...]) @ B[tdx, ...]
            self.data.at[tdx, "excitation"] = G[0, :]
            self.data.at[tdx, "inhibition"] = G[1, :]
        return self.data

    def estimate(self):
        self.scale_data()
        self.filter_data()
        self.compute_polarizations()
        self.compute_active_conductance_constants()
        self.compute_active_currents()
        self.compute_leakage_currents()
        self.compute_membrane_currents()
        self.filter_membrane_currents()
        self.compute_passive_conductances()
        return self.data

    def plot(self):
        fig, ax = plt.subplots(5, 1, sharex=True)
        print(self.data)
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            ax[0].plot(self.data["times"], self.data[clamp])
            ax[1].plot(self.data["times"], self.data["Ileak_"+str(clamp)])
            ax[2].plot(self.data["times"], self.data["filtered_Im_"+str(clamp)])
            ax[3].plot(self.data["times"], self.data["Iactive_"+str(clamp)])
        ax[4].plot(self.data["times"], self.data["excitation"], c="r")
        ax[4].plot(self.data["times"], self.data["inhibition"], c="b")
        ax[4].set_xlabel("times(sec)")
        ax[0].set_ylabel("Vm(V)")
        ax[1].set_ylabel("Ileak(A)")
        ax[2].set_ylabel("Im(A)")
        ax[3].set_ylabel("Iactive(A)")
        ax[4].set_ylabel("G(s)")
        plt.show()
        pass

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
        return filters.sosfiltfilt(second_order_sections, input)