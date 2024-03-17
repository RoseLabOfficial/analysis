import numpy as np
import pandas as pd
import scipy.signal as filters
import matplotlib.pyplot as plt
import libs.signals as sigs
from pathlib import Path
from libs.sysops import SystemOperations
from libs.readers import XLReader

class NonLinearModel:
    def __init__(self, parameters: pd.DataFrame, tstart: float, tend: float, sampling_rate: float, dtype: np.dtype = np.float64):
        self.parameters = parameters
        self.timeseries = sigs.TimeSeries(tstart, tend, sampling_rate, dtype)
        self.times = self.timeseries.make_time()
        pass

    def make_voltage_grid(self):
        nsteps = self.timeseries.get_total_samples()
        membrane_potential = np.zeros((nsteps, self.parameters["Iinj"].shape[0]), dtype=self.timeseries.dtype)
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            membrane_potential[:, idx] = np.linspace(self.parameters["Er"][idx], self.parameters["Et"][idx], nsteps, dtype=self.timeseries.dtype)
        return membrane_potential
    
    def compute_active_conductance_constants(self):
        self.parameters["alpha"] = (1.0/self.parameters["Rin"])/(2.0*(self.parameters["Eact"] - self.parameters["Ess"]))
        self.parameters["beta"] = self.parameters["alpha"]*(self.parameters["Et"] - self.parameters["Ess"])
        self.parameters["alpha"] = self.parameters["alpha"]*self.parameters["xalpha"]
        self.parameters["beta"] = self.parameters["beta"]*self.parameters["xbeta"]
        return self.parameters
    
    def compute_current(self, membrane_potential: np.ndarray):
        self.compute_active_conductance_constants()
        alpha = self.parameters["alpha"].to_numpy()
        beta = self.parameters["beta"].to_numpy()
        resting_potential = self.parameters["Er"].to_numpy()
        threshold_potential = self.parameters["Et"].to_numpy()
        input_resistance = self.parameters["Rin"].to_numpy()
        return -alpha*(membrane_potential - resting_potential)*(membrane_potential - threshold_potential) - \
                    beta*(membrane_potential - resting_potential) + (1.0/input_resistance)*(membrane_potential - resting_potential)

class LowPassFilter:
    def __init__(self, sampling_rate: float):
        self.sampling_rate = sampling_rate
        pass

    def compute_minimum_order(self, cutoff: float, stopband_attenuation: float, passband_ripple: float):
        passband_edge_frequency = cutoff/(self.sampling_rate/2)
        stopband_edge_frequency = 5*passband_edge_frequency
        order, normalized_cutoff_frequency = filters.buttord(passband_edge_frequency, stopband_edge_frequency, 
                                                            passband_ripple, stopband_attenuation)
        return order, normalized_cutoff_frequency
    
    def filter(self, input: float, cutoff: float, stopband_attenuation: float, passband_ripple: float):
        order, normalized_frequency = self.compute_minimum_order(cutoff, stopband_attenuation, passband_ripple)
        second_order_sections = filters.butter(order, normalized_frequency, output="sos")
        return filters.sosfiltfilt(second_order_sections, input)

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
    
    def compute_activation_conductance_constants(self):
        self.parameters["alpha"] = (1.0/self.parameters["Rin"])/(2.0*(self.parameters["Eact"] - self.parameters["Ess"]))
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
    
    def compute_activation_currents(self):
        self.compute_activation_conductance_constants()
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            alpha_current = self.parameters["alpha"][idx]*(self.data[clamp] - self.parameters["Er"][idx])*(self.data[clamp] - self.parameters["Et"][idx])
            beta_current = self.parameters["beta"][idx]*(self.data[clamp] - self.parameters["Er"][idx])
            Iact = alpha_current + beta_current
            Iact[self.data[clamp] < self.parameters["Ess"][idx]] = 0.0
            Iact[self.data[clamp] > self.parameters["Eact"][idx]] = 0.0
            self.data["Iactive_"+str(clamp)] = Iact
        return self.data
    
    def compute_membrane_currents(self):
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            self.data["Im_"+str(clamp)] = self.parameters["Cm"][idx]*(self.data[clamp].diff()/self.data["times"].diff())
            self.data.at[0, "Im_"+str(clamp)] = 0.0
            self.data["Im_"+str(clamp)] = self.data["Im_"+str(clamp)] - self.data["Im_"+str(clamp)][0]
        return self.data
    
    def filter_membrane_currents(self, cutoff: float = 40.0, stopband_attenuation: float = 80, passband_ripple: float = 0.01):
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        filter = LowPassFilter(sampling_rate)
        for idx, inj in enumerate(self.parameters["Iinj"]):
            self.data["filtered_Im_"+str(inj)] = filter.filter(self.data["Im_"+str(inj)], cutoff, stopband_attenuation, passband_ripple)
        return self.data
    
    def filter_activation_currents(self, cutoff: float = 40.0, stopband_attenuation: float = 80, passband_ripple: float = 0.01):
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        filter = LowPassFilter(sampling_rate)
        for idx, inj in enumerate(self.parameters["Iinj"]):
            Iact = filter.filter(self.data["Iactive_"+str(inj)], cutoff, stopband_attenuation, passband_ripple)
            # Iact = Iact - Iact[0]
            # Iact[Iact<0] = 0.0
            self.data["filtered_Iactive_"+str(inj)] = Iact
        return self.data
    
    def compute_passive_conductances(self):
        ntimesteps = self.data.shape[0]
        A = np.zeros((ntimesteps, 2, 2))
        B = np.zeros((ntimesteps, 2, 1))
        Vm = self.data[[x for x in self.parameters["Iinj"]]].to_numpy()
        Im = self.data[["filtered_Im_"+str(x) for x in self.parameters["Iinj"]]].to_numpy()
        Iact = self.data[["filtered_Iactive_"+str(x) for x in self.parameters["Iinj"]]].to_numpy()
        Ileak = self.data[["Ileak_"+str(x) for x in self.parameters["Iinj"]]].to_numpy()
        Ee = self.parameters["Ee"].to_numpy()
        Ei = self.parameters["Ei"].to_numpy()
        Iinj = self.parameters["Iinj"].to_numpy()
        A[:, 0, 0] = np.sum(np.square(Vm - Ee), axis=1)
        A[:, 0, 1] = np.sum((Vm - Ee) * (Vm- Ei), axis=1)
        A[:, 1, 0] = A[:, 0, 1]
        A[:, 1, 1] = np.sum(np.square(Vm - Ei), axis=1)
        B[:, 0, 0] = -1.0*np.sum((Im - Iact - Iinj + Ileak)*(Vm - Ee), axis=1)
        B[:, 1, 0] = -1.0*np.sum((Im - Iact - Iinj + Ileak)*(Vm - Ei), axis=1)
        G = np.linalg.pinv(A) @ B
        self.data["excitation"] = G[:, 0, 0]
        self.data["inhibition"] = G[:, 1, 0]
        return self.data

    def estimate_conductances(self, membrane_voltage_filter_cutoff: float = 200, activation_currents_filter_cutoff: float = 100, membrane_currents_filter_cuoff: float=100):
        self.scale_data()
        self.filter_data(cutoff=membrane_voltage_filter_cutoff)
        self.compute_polarizations()
        self.compute_activation_currents()
        self.filter_activation_currents(cutoff=activation_currents_filter_cutoff)
        self.compute_leakage_currents()
        self.compute_membrane_currents()
        self.filter_membrane_currents(cutoff=membrane_currents_filter_cuoff)
        self.compute_passive_conductances()
        return self.data

    def plot(self, plot_save_path: str = "./plot.png"):
        fig, ax = plt.subplots(5, 1, sharex=True)
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            ax[0].plot(self.data["times"], self.data[clamp])
            ax[1].plot(self.data["times"], self.data["Ileak_"+str(clamp)])
            ax[2].plot(self.data["times"], self.data["filtered_Im_"+str(clamp)])
            ax[3].plot(self.data["times"], self.data["filtered_Iactive_"+str(clamp)])
        ax[4].plot(self.data["times"], self.data["excitation"], c="r")
        ax[4].plot(self.data["times"], self.data["inhibition"], c="b")
        ax[4].plot(self.data["times"], self.data["times"]*0.0, "--k")
        ax[4].set_xlabel("times(sec)")
        ax[0].set_ylabel("Vm(V)")
        ax[1].set_ylabel("Ileak(A)")
        ax[2].set_ylabel("Im(A)")
        ax[3].set_ylabel("Iactive(A)")
        ax[4].set_ylabel("G(s)")
        plt.savefig(plot_save_path)
        pass

class Analyzer:
    def __init__(self, filepaths: Path, output_path: Path):
        self.sysops = SystemOperations()
        self.filepaths = self.sysops.split_file_address(filepaths)
        self.output_path = output_path

    def analyze(self, fileaddress):
        filepath = fileaddress[0]
        filename = fileaddress[1]
        fileformat = fileaddress[2]
        print(f"Reading: {filename}")
        reader = XLReader(Path(filepath) / (filename + fileformat))
        store_path = self.sysops.make_timed_directory(self.output_path, filename)
        recordings = {}
        for paradigm in reader.get_paradigms():
            recordings[paradigm] = WholeCellRecording(
                reader.get_paradigm_data(paradigm), 
                reader.get_paradigm_parameters(paradigm)
            )
            recordings[paradigm].estimate_conductances()
            self.sysops.make_directory(store_path / paradigm)
            recordings[paradigm].data.to_csv(store_path / (paradigm+"/analysis.csv"), index=False)
            recordings[paradigm].parameters.to_csv(store_path / (paradigm+"/parameters.csv"), index=False)
            print(f"{paradigm} done")
        print(f"{filename} done")
        return recordings

    def plot_dev(self, recordings):
        fig, axs = plt.subplots(nrows = 5, ncols = len(recordings), sharex="all", sharey="row", figsize=(15, 10), constrained_layout=True)
        for idx, paradigm in enumerate(recordings):
            Vm = recordings[paradigm].data[[x for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            Im = recordings[paradigm].data[["filtered_Im_"+str(x) for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            Ileak = recordings[paradigm].data[["Ileak_"+str(x) for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            Iact = recordings[paradigm].data[["filtered_Iactive_"+str(x) for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            G = recordings[paradigm].data[["excitation", "inhibition"]].to_numpy()
            times = recordings[paradigm].data["times"].to_numpy()
            Er = times*0 + recordings[paradigm].parameters["Er"][0]
            axs[0, idx].plot(times, Vm)
            axs[0, idx].plot(times, Er, '--k')
            axs[0, idx].set_ylabel("Vm (V)")
            axs[0, idx].grid(True)
            axs[1, idx].plot(times, Im)
            axs[1, idx].set_ylabel("Im (A)")
            axs[1, idx].grid(True)
            axs[2, idx].plot(times, Ileak)
            axs[2, idx].set_ylabel("Ileak (A)")
            axs[2, idx].grid(True)
            axs[3, idx].plot(times, Iact)
            axs[3, idx].set_ylabel("Iact (A)")
            axs[3, idx].grid(True)
            axs[4, idx].plot(times, G[:, 0], c='r')
            axs[4, idx].plot(times, G[:, 1], c='b')
            axs[4, idx].plot(times, times*0, '--k')
            axs[4, idx].set_ylabel("G (S)")
            axs[4, idx].set_xlabel("times(sec)")
            axs[4, idx].grid(True)
        plt.show()
        pass

    def plot_prod(self, recordings):
        fig, axs = plt.subplots(nrows = 5, ncols = len(recordings), sharex="all", sharey="row", figsize=(15, 10), constrained_layout=True)
        for idx, paradigm in enumerate(recordings):
            Vm = recordings[paradigm].data[[x for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            Im = recordings[paradigm].data[["filtered_Im_"+str(x) for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            Ileak = recordings[paradigm].data[["Ileak_"+str(x) for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            Iact = recordings[paradigm].data[["filtered_Iactive_"+str(x) for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            G = recordings[paradigm].data[["excitation", "inhibition"]].to_numpy()
            times = recordings[paradigm].data["times"].to_numpy()
            Er = times*0 + recordings[paradigm].parameters["Er"][0]
            axs[0, idx].plot(times, Vm)
            axs[0, idx].plot(times, Er, '--k')
            axs[0, idx].set_ylabel("Vm (V)")
            axs[0, idx].grid(True)
            axs[1, idx].plot(times, Im)
            axs[1, idx].set_ylabel("Im (A)")
            axs[1, idx].grid(True)
            axs[2, idx].plot(times, Ileak)
            axs[2, idx].set_ylabel("Ileak (A)")
            axs[2, idx].grid(True)
            axs[3, idx].plot(times, Iact)
            axs[3, idx].set_ylabel("Iact (A)")
            axs[3, idx].grid(True)
            axs[4, idx].plot(times, G[:, 0], c='r')
            axs[4, idx].plot(times, G[:, 1], c='b')
            axs[4, idx].plot(times, times*0, '--k')
            axs[4, idx].set_ylabel("G (S)")
            axs[4, idx].set_xlabel("times(sec)")
            axs[4, idx].grid(True)
        plt.show()
        pass

    def run(self):
        for filepath in self.filepaths:
            recordings = self.analyze(filepath)
            self.plot_dev(recordings)
        pass