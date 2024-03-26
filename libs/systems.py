import numpy as np
import pandas as pd
import scipy.signal as filters
from scipy.optimize import linprog, minimize, Bounds
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
    def __init__(self, parameters):
        self.passband = parameters["passband"]
        self.stopband = parameters["stopband"]
        self.attenuation = parameters["attenuation"]
        self.ripple = parameters["ripple"]
        pass

    def compute_minimum_order(self, sampling_rate: float):
        normalizing_frequency = sampling_rate/2
        normalized_passband = self.passband/normalizing_frequency
        normalized_stopband = self.stopband/normalizing_frequency
        order, normalized_cutoff_frequency = filters.buttord(normalized_passband, 
                                                            normalized_stopband, 
                                                            self.ripple, 
                                                            self.attenuation)
        return order, normalized_cutoff_frequency
    
    def append_samples(self, signal):
        nsamples = signal.shape[0]
        nappend = int(nsamples/2)
        append_samples = np.zeros((nappend,), dtype=signal.dtype)+signal[0]
        return np.concatenate([append_samples, signal, append_samples], axis=0), nappend
    
    def deppend_samples(self, signal, nappend):
        return signal[nappend:-nappend, ...]
    
    def propagate(self, input: float, sampling_rate: float):
        order, normalized_frequency = self.compute_minimum_order(sampling_rate)
        second_order_sections = filters.butter(order, normalized_frequency, output="sos")
        input_signal, nappend = self.append_samples(input)
        output_signal = filters.sosfiltfilt(second_order_sections, input_signal)
        return self.deppend_samples(output_signal, nappend)

class WholeCellRecording:
    def __init__(self, data: pd.DataFrame, parameters: pd.DataFrame, filter_parameters):
        self.data = data
        self.parameters = parameters
        self.filters = {}
        self.filters["membrane_potentials"] = LowPassFilter(filter_parameters["membrane_potentials"])
        self.filters["membrane_currents"] = LowPassFilter(filter_parameters["membrane_currents"])
        self.filters["activation_currents"] = LowPassFilter(filter_parameters["activation_currents"])
        self.scale_data()
        pass

    def scale_data(self):
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            self.data[clamp] = ((self.data[clamp] - self.data[clamp][0])*1e-3)+self.parameters["Ess"][idx]
        return self.data
    
    def filter_membrane_potentials(self):
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        for inj in self.parameters["Iinj"]:
            self.data[inj] = self.filters["membrane_potentials"].propagate(self.data[inj], sampling_rate)
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
    
    def filter_membrane_currents(self):
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        for inj in self.parameters["Iinj"]:
            self.data["filtered_Im_"+str(inj)] = self.filters["membrane_currents"].propagate(self.data["Im_"+str(inj)], sampling_rate)
        return self.data
    
    def filter_activation_currents(self):
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        for inj in self.parameters["Iinj"]:
            Iact = self.filters["activation_currents"].propagate(self.data["Iactive_"+str(inj)], sampling_rate)
            self.data["filtered_Iact_"+str(inj)] = Iact
        return self.data

    def objective(self, solution):
        self.parameters["Eact"] = solution
        self.compute_activation_currents()
        self.filter_activation_currents()
        self.compute_passive_conductances()
        self.compute_stats()
        neg_excitation = self.data["excitation"][self.data["excitation"] < 0]
        neg_excitation = neg_excitation/np.amin(neg_excitation)
        neg_inhibition = self.data["inhibition"][self.data["inhibition"] < 0]
        neg_inhibition = neg_inhibition/np.amin(neg_inhibition)
        return np.square(np.sum(neg_excitation) + np.sum(neg_inhibition))
    
    def optimize(self):
        self.estimate_conductances()
        def obj_wrapper(solution):
            return self.objective(solution)
        
        Emax = [max_val for max_val in self.data[[x for x in self.parameters["Iinj"]]].max()]
        Ess = np.asarray([x for x in self.parameters["Ess"]])
        # Emax = np.asarray([x for x in self.parameters["Et"]])
        # initial_guess = (Ess + Emax)/2
        bounds = Bounds(Ess, Emax)
        result = minimize(obj_wrapper, Emax, method='trust-constr', bounds=bounds)
        self.parameters["Eact"] = result.x
        return result

    
    def compute_passive_conductances(self):
        ntimesteps = self.data.shape[0]
        A = np.zeros((ntimesteps, 2, 2))
        B = np.zeros((ntimesteps, 2, 1))
        Vm = self.data[[x for x in self.parameters["Iinj"]]].to_numpy()
        Im = self.data[["filtered_Im_"+str(x) for x in self.parameters["Iinj"]]].to_numpy()
        Iact = self.data[["filtered_Iact_"+str(x) for x in self.parameters["Iinj"]]].to_numpy()
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
        self.data["positive_excitation"] = self.data["excitation"]
        self.data["positive_inhibition"] = self.data["inhibition"]
        self.data.loc[self.data["positive_excitation"] < 0, "positive_excitation"] = 0.0
        self.data.loc[self.data["positive_inhibition"] < 0, "positive_inhibition"] = 0.0
        self.data["resultant_excitation"] = self.data["positive_excitation"] - self.data["positive_inhibition"]
        self.data["resultant_inhibition"] = self.data["positive_inhibition"] - self.data["positive_excitation"]
        self.data.loc[self.data["resultant_excitation"] < 0, "resultant_excitation"] = 0.0
        self.data.loc[self.data["resultant_inhibition"] < 0, "resultant_inhibition"] = 0.0
        return self.data
    
    def get_clamp_near_0(self):
        min_Iinj_index = self.parameters["Iinj"].abs().argmin()
        min_Iinj = self.parameters["Iinj"][min_Iinj_index]
        return min_Iinj_index, min_Iinj
    
    def compute_stats(self):
        min_Iinj_index, min_Iinj = self.get_clamp_near_0()
        stats = pd.DataFrame()
        paradigm_all_var_stats = self.data.mean(numeric_only=True).to_frame().T
        stats["depolarization"] = paradigm_all_var_stats[f"depolarization_{min_Iinj}"]
        stats["hyperpolarization"] = paradigm_all_var_stats[f"hyperpolarization_{min_Iinj}"]
        stats["Im"] = paradigm_all_var_stats[f"filtered_Im_{min_Iinj}"]
        stats["Ileak"] = paradigm_all_var_stats[f"Ileak_{min_Iinj}"]
        stats["Iact"] = paradigm_all_var_stats[f"filtered_Iact_{min_Iinj}"]
        stats["mean_excitation"] = paradigm_all_var_stats["positive_excitation"]
        stats["mean_inhibition"] = paradigm_all_var_stats["positive_inhibition"]
        stats["net_excitation"] = paradigm_all_var_stats["resultant_excitation"]
        stats["net_inhibition"] = paradigm_all_var_stats["resultant_inhibition"]
        stats["sps"] = self.parameters["sps"][min_Iinj_index]
        return stats

    def estimate_conductances(self):
        self.filter_membrane_potentials()
        self.compute_polarizations()
        self.compute_activation_currents()
        self.filter_activation_currents()
        self.compute_leakage_currents()
        self.compute_membrane_currents()
        self.filter_membrane_currents()
        self.compute_passive_conductances()
        self.stats = self.compute_stats()
        return self.data

    def plot(self, plot_save_path: str = "./plot.png"):
        fig, ax = plt.subplots(5, 1, sharex=True)
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            ax[0].plot(self.data["times"], self.data[clamp])
            ax[1].plot(self.data["times"], self.data["Ileak_"+str(clamp)])
            ax[2].plot(self.data["times"], self.data["filtered_Im_"+str(clamp)])
            ax[3].plot(self.data["times"], self.data["filtered_Iact_"+str(clamp)])
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
    def __init__(self, filepaths: Path, output_path: Path, filter_configurations: dict):
        self.sysops = SystemOperations()
        self.filepaths = self.sysops.split_file_address(filepaths)
        self.output_path = output_path
        self.filter_configurations = filter_configurations

    def analyze(self, fileaddress):
        filepath = fileaddress[0]
        filename = fileaddress[1]
        fileformat = fileaddress[2]
        print(f"Reading: {filename}")
        reader = XLReader(Path(filepath) / (filename + fileformat))
        recordings = {}
        for idx, paradigm in enumerate(reader.get_paradigms()):
            recordings[paradigm] = WholeCellRecording(
                reader.get_paradigm_data(paradigm), 
                reader.get_paradigm_parameters(paradigm),
                self.filter_configurations
            )
        max_depols = np.zeros(len(recordings))
        for idx, paradigm in enumerate(recordings):
            min_idx_Iinj, min_Iinj = recordings[paradigm].get_clamp_near_0()
            vms = recordings[paradigm].data[min_Iinj]
            depols = vms - recordings[paradigm].parameters["Ess"][min_idx_Iinj]
            max_depols[idx] = np.max(depols)
        max_index = np.argmax(max_depols)
        max_paradigm = list(recordings.keys())[max_index]
        print(f"Optimizing Eact for {max_paradigm}")
        recordings[max_paradigm].optimize()
        Eact = recordings[max_paradigm].parameters["Eact"]
        print(f"New Eacts: {Eact.to_numpy().tolist()}")
        recordings[max_paradigm].stats.insert(0, "paradigm", max_paradigm)
        overall_stats = recordings[max_paradigm].stats.copy()
        for idx, paradigm in enumerate(recordings):
            if not paradigm == max_paradigm:
                recordings[paradigm].parameters["Eact"] = Eact
                recordings[paradigm].estimate_conductances()
                recordings[paradigm].stats.insert(0, "paradigm", paradigm)
                overall_stats = pd.concat([overall_stats, recordings[paradigm].stats], axis=0)
            print(f"{paradigm} done")
        overall_stats = overall_stats.sort_values(by="paradigm", ascending=True)
        result_filename = self.output_path / f"{filename}_analyzed{self.sysops.get_time_as_string()}"
        # self.write_to_excel(f"{result_filename}.xlsx", recordings, overall_stats)
        print(f"{filename} done \n")
        return recordings, overall_stats, result_filename
    
    def write_to_excel(self, filepath: Path, recordings, stats):
        if not self.sysops.check_directory(filepath):
            pd.DataFrame().to_excel(filepath)
        with pd.ExcelWriter(filepath, mode='a', engine='openpyxl', if_sheet_exists='new') as writer:
            for paradigm in recordings:
                recordings[paradigm].data.to_excel(writer, sheet_name=paradigm, index=False)
                recordings[paradigm].parameters.to_excel(writer, sheet_name="parameters_"+paradigm, index = False)
            stats.to_excel(writer, sheet_name="stats", index=False)
        pass

    def write_analysis_to_excel(self, filepath: Path, paradigm: str, paradigm_data: pd.DataFrame):
        if self.sysops.check_directory(filepath):
            with pd.ExcelWriter(filepath, mode='a', engine='openpyxl', if_sheet_exists='replace') as writer:
                paradigm_data.to_excel(writer, sheet_name=paradigm, index=False)
        else:
            # File does not exist, create a new file
            paradigm_data.to_excel(filepath, sheet_name=paradigm, index=False)
        pass

    def plot_dev(self, recordings, filename: Path):
        fig, axs = plt.subplots(nrows = 7, ncols = len(recordings), sharex="all", sharey="row", figsize=(15, 10), constrained_layout=True)
        for idx, paradigm in enumerate(recordings):
            if "representative" in recordings[paradigm].data:
                rep = recordings[paradigm].data["representative"].to_numpy()
            Vm = recordings[paradigm].data[[x for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            Im = recordings[paradigm].data[["filtered_Im_"+str(x) for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            Ileak = recordings[paradigm].data[["Ileak_"+str(x) for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            Iact = recordings[paradigm].data[["filtered_Iact_"+str(x) for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            G = recordings[paradigm].data[["excitation", "inhibition"]].to_numpy()
            times = recordings[paradigm].data["times"].to_numpy()
            Er = times*0 + recordings[paradigm].parameters["Er"][0]
            Et = times*0 + recordings[paradigm].parameters["Et"][0]
            Eact = Vm*0 + recordings[paradigm].parameters["Eact"].to_numpy()
            if "stimulus" in recordings[paradigm].data:
                stim = recordings[paradigm].data["stimulus"].to_numpy()
            axs[0, idx].set_title(paradigm)
            if "representative" in recordings[paradigm].data:
                axs[0, idx].plot(times, rep)
                axs[0, idx].plot(times, Er, '--k')
            axs[0, idx].set_ylabel("Rep. Vm (V)")
            axs[0, idx].set_title(paradigm)
            curves = axs[1, idx].plot(times, Vm)
            colors = [x.get_color() for x in curves]
            axs[1, idx].plot(times, Er, '--k')
            axs[1, idx].plot(times, Et, linestyle='--', color=(0.5, 0.5, 0.5))
            for i in range(Eact.shape[-1]):
                axs[1, idx].plot(times, Eact[:, i], linestyle='--', color=colors[i])
            axs[1, idx].set_ylabel("Vm (V)")
            axs[1, idx].grid(True)
            axs[2, idx].plot(times, Im)
            axs[2, idx].set_ylabel("Im (A)")
            axs[2, idx].grid(True)
            axs[3, idx].plot(times, Ileak)
            axs[3, idx].set_ylabel("Ileak (A)")
            axs[3, idx].grid(True)
            axs[4, idx].plot(times, Iact)
            axs[4, idx].set_ylabel("Iact (A)")
            axs[4, idx].grid(True)
            # axs[4, idx].set_title(", ".join([f"{val:.3f}" for val in recordings[paradigm].parameters["Eact"]]))
            axs[5, idx].plot(times, G[:, 0], c='r')
            axs[5, idx].plot(times, G[:, 1], c='b')
            axs[5, idx].plot(times, times*0, '--k')
            axs[5, idx].set_ylabel("G (S)")
            axs[5, idx].grid(True)
            if "stimulus" in recordings[paradigm].data:
                axs[6, idx].plot(times, stim)
            axs[6, idx].set_ylabel("Stimulus")
            axs[6, idx].set_xlabel("times(sec)")
            axs[6, idx].grid(True)
        plt.savefig(str(filename)+f"_dev_traces.png")
        pass

    def set_stats_scale(self, ax, scale_max, margin=0.1):
        scale_max = scale_max + margin*scale_max
        ax.set_ylim([-1*scale_max, scale_max])
        pass

    def plot_stats_dev(self, recordings, filename: Path):
        fig, axs = plt.subplots(nrows = 5, ncols = 1, sharex="all", figsize=(15, 10), constrained_layout=True)
        mean_depolarizations = [recordings[paradigm].stats["depolarization"] for paradigm in recordings]
        mean_hyperpolarizations = [recordings[paradigm].stats["hyperpolarization"] for paradigm in recordings]
        mean_Im = [recordings[paradigm].stats["Im"] for paradigm in recordings]
        mean_Ileak = [recordings[paradigm].stats["Ileak"] for paradigm in recordings]
        mean_Iactivation = [recordings[paradigm].stats["Iact"] for paradigm in recordings]
        mean_excitation = [recordings[paradigm].stats["mean_excitation"] for paradigm in recordings]
        mean_inhibition = [recordings[paradigm].stats["mean_inhibition"] for paradigm in recordings]
        net_excitation = [recordings[paradigm].stats["net_excitation"] for paradigm in recordings]
        net_inhibition = [recordings[paradigm].stats["net_inhibition"] for paradigm in recordings]
        sps = [recordings[paradigm].stats["sps"] for paradigm in recordings]
        paradigms = [paradigm for paradigm in recordings]
        xlocations = np.arange(len(paradigms))
        # for paradigm in paradigms:
        #     mean_depolarizations.append(recordings[paradigm].stats["depolarization"])
        #     mean_hyperpolarizations.append(recordings[paradigm].stats["hyperpolarization"])
        #     mean_Im.append(recordings[paradigm].stats["Im"])
        #     mean_Iactivation.append(recordings[paradigm].stats["Iact"])
        #     mean_Ileak.append(recordings[paradigm].stats["Ileak"])
        #     mean_excitation.append(recordings[paradigm].stats["mean_excitation"])
        #     mean_inhibition.append(recordings[paradigm].stats["mean_inhibition"])
        #     net_excitation.append(recordings[paradigm].stats["net_excitation"])
        #     net_inhibition.append(recordings[paradigm].stats["net_inhibition"])
        #     sps.append(recordings[paradigm].parameters["sps"])
        # mean_depolarizations = np.asarray(mean_depolarizations)
        # mean_hyperpolarizations = np.asarray(mean_hyperpolarizations)
        # mean_Im = np.asarray(mean_Im)
        # mean_Iactivation = np.asarray(mean_Iactivation)
        # mean_Ileak = np.asarray(mean_Ileak)
        # mean_excitation = np.asarray(mean_excitation)
        # mean_inhibition = np.asarray(mean_inhibition)
        # sps = np.asarray(sps)
        axs[0].bar(xlocations, mean_depolarizations, align='center', color="red")
        axs[0].bar(xlocations, mean_hyperpolarizations, align='center', color="blue")
        axs[0].axhline(0, color='grey', linewidth=0.8)
        axs[0].set_ylabel("polarizations (V)")
        scale_max = np.amax([np.amax(mean_depolarizations), np.amax(mean_hyperpolarizations)])
        self.set_stats_scale(axs[0], scale_max, 0.1)
        axs[1].bar(xlocations, mean_Im, align='center', color="black")
        axs[1].axhline(0, color='grey', linewidth=0.8)
        axs[1].set_ylabel("Im (A)")
        scale_max = np.amax(mean_Im)
        self.set_stats_scale(axs[1], scale_max, 0.1)
        axs[2].bar(xlocations, mean_Iactivation, align='center', color="black")
        axs[2].axhline(0, color='grey', linewidth=0.8)
        axs[2].set_ylabel("Iact (A)")
        scale_max = np.amax(mean_Iactivation)
        self.set_stats_scale(axs[2], scale_max, 0.1)
        axs[3].bar(xlocations, mean_Ileak, align='center', color="black")
        axs[3].axhline(0, color='grey', linewidth=0.8)
        axs[3].set_ylabel("Ileak (A)")
        scale_max = np.amax(mean_Ileak)
        self.set_stats_scale(axs[3], scale_max, 0.1)
        axs[4].bar(xlocations, mean_excitation, align='center', color="red")
        axs[4].bar(xlocations, -1*mean_inhibition, align='center', color="blue")
        axs[4].axhline(0, color='grey', linewidth=0.8)
        scale_max = np.amax([np.amax(mean_excitation), np.amax(mean_inhibition)])
        self.set_stats_scale(axs[4], scale_max, 0.1)
        ax = axs[4].twinx()
        ax.plot(xlocations, sps, color='k', marker = 'o')
        scale_max = np.amax(sps)
        self.set_stats_scale(ax, scale_max, 0.1)
        axs[4].set_ylabel("G (S)")
        plt.savefig(str(filename)+f"_dev_stats.png")
        pass

    def run(self):
        for filepath in self.filepaths:
            recordings, stats, result_filename = self.analyze(filepath)
            # self.write_to_excel(f"{result_filename}.xlsx", recordings, stats)
            self.plot_dev(recordings, result_filename)
            self.plot_stats_dev(recordings, result_filename)
        pass