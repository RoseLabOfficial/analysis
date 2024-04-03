from . import *
import libs.signals as sigs
from pathlib import Path
from libs.readers import XLReader
from libs.utils import *

class NonLinearModel:
    def __init__(self, parameters: pd.DataFrame, tstart: float, tend: float, sampling_rate: float, dtype: np.dtype = np.float64):
        self.parameters = parameters
        self.timeseries = sigs.TimeSeries(tstart, tend, sampling_rate, dtype)
        self.times = self.timeseries.make_time()
        pass

    def make_voltage_grid(self):
        nsteps = self.timeseries.get_total_samples()
        membrane_potential = np.zeros((nsteps, self.parameters["Iinj"].shape[0]), dtype=self.timeseries.dtype)
        for idx, _ in enumerate(self.parameters["Iinj"]):
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
    def __init__(self, parameters, name:str = "generic"):
        self.passband = parameters["passband"]
        self.stopband = parameters["stopband"]
        self.attenuation = parameters["attenuation"]
        self.ripple = parameters["ripple"]
        if self.stopband < self.passband:
            lpf_logger.debug(f"For {name} low pass filter stopband cannot be less than passband")
            raise ValueError(f"For {name} low pass filter stopband cannot be less than passband")
        self.name = name
        pass        

    def compute_minimum_order(self, sampling_rate: float, log=False):
        normalizing_frequency = sampling_rate/2
        normalized_passband = self.passband/normalizing_frequency
        normalized_stopband = self.stopband/normalizing_frequency
        if log:
            lpf_logger.info(f"{self.name} low pass filter ")
        order, normalized_cutoff_frequency = filters.buttord(normalized_passband, 
                                                            normalized_stopband, 
                                                            self.ripple, 
                                                            self.attenuation)
        if log:
            lpf_logger.info(f"{self.name} low pass filter computed minimum order: {order} with transition gap: {self.stopband - self.passband}")
        return order, normalized_cutoff_frequency
    
    def append_samples(self, signal):
        nsamples = signal.shape[0]
        nappend = int(nsamples/2)
        append_samples = np.zeros((nappend,), dtype=signal.dtype)+signal[0]
        return np.concatenate([append_samples, signal, append_samples], axis=0), nappend
    
    def deppend_samples(self, signal, nappend):
        return signal[nappend:-nappend, ...]
    
    def propagate(self, input: float, sampling_rate: float, log=False):
        order, normalized_frequency = self.compute_minimum_order(sampling_rate, log)
        second_order_sections = filters.butter(order, normalized_frequency, output="sos")
        input_signal, nappend = self.append_samples(input)
        output_signal = filters.sosfiltfilt(second_order_sections, input_signal)
        return self.deppend_samples(output_signal, nappend)

class WholeCellRecording:
    def __init__(self, data: pd.DataFrame, parameters: pd.DataFrame, filter_parameters):
        self.data = data
        self.parameters = parameters
        self.filters = {}
        self.filters["membrane_potentials"] = LowPassFilter(filter_parameters["membrane_potentials"], name="membrane_potentials")
        self.filters["membrane_currents"] = LowPassFilter(filter_parameters["membrane_currents"], name="membrane_currents")
        self.filters["activation_currents"] = LowPassFilter(filter_parameters["activation_currents"], name="activation_currents")
        self.complier = Compliance()
        self.scale_data()
        pass

    def check_compliance(self):
        if self.complier.check_compliance(self.data, self.parameters):
            return True
        else:
            wholecell_logger.debug(f"The data / parameters do not comply with standards.")
        return False

    def scale_data(self, log=False):
        if log:
            wholecell_logger.info("Scaling membrane voltage")
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            self.data[clamp] = ((self.data[clamp] - self.data[clamp][0])*1e-3)+self.parameters["Ess"][idx]
        return self.data
    
    def filter_membrane_potentials(self, log=False):
        if log:
            wholecell_logger.info("Filtering membrane potentials")
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        for inj in self.parameters["Iinj"]:
            self.data[inj] = self.filters["membrane_potentials"].propagate(self.data[inj], sampling_rate, log)
        return self.data
    
    def compute_activation_conductance_constants(self, log=False):
        if log:
            wholecell_logger.info("Computing activation constants alpha and beta")
        self.parameters["alpha"] = (1.0/self.parameters["Rin"])/(2.0*(self.parameters["Eact"] - self.parameters["Ess"]))
        self.parameters["beta"] = self.parameters["alpha"]*(self.parameters["Et"] - self.parameters["Ess"])
        self.parameters["alpha"] = self.parameters["alpha"]*self.parameters["xalpha"]
        self.parameters["beta"] = self.parameters["beta"]*self.parameters["xbeta"]
        if log:
            wholecell_logger.info(f"alpha = {self.parameters["alpha"].to_numpy().tolist()} & beta = {self.parameters["beta"].to_numpy().tolist()}")
        return self.parameters
    
    def compute_polarizations(self, log=False):
        if log:
            wholecell_logger.info("Computing polarizations")
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            self.data["depolarization_"+str(clamp)] = np.where(self.data[clamp] > self.parameters["Ess"][idx], self.data[clamp] - self.parameters["Ess"][idx], 0)
            self.data["hyperpolarization_"+str(clamp)] = np.where(self.data[clamp] < self.parameters["Ess"][idx], self.data[clamp] - self.parameters["Ess"][idx], 0)
        return self.data
    
    def compute_leakage_currents(self, log=False):
        if log:
            wholecell_logger.info("Computing leakage currents")
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            self.data["Ileakage_"+str(clamp)] = (1/self.parameters["Rin"][idx])*(self.data[clamp] - self.parameters["Er"][idx])
        return self.data
    
    def compute_activation_currents(self, log=False):
        if log:
            wholecell_logger.info("Computing activation currents")
        self.compute_activation_conductance_constants(log)
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            alpha_current = self.parameters["alpha"][idx]*(self.data[clamp] - self.parameters["Er"][idx])*(self.data[clamp] - self.parameters["Et"][idx])
            beta_current = self.parameters["beta"][idx]*(self.data[clamp] - self.parameters["Er"][idx])
            activation_current = alpha_current + beta_current
            activation_current[self.data[clamp] < self.parameters["Ess"][idx]] = 0.0
            activation_current[self.data[clamp] > self.parameters["Eact"][idx]] = 0.0
            self.data["Iactivation_"+str(clamp)] = activation_current
        return self.data
    
    def compute_membrane_currents(self, log=False):
        if log:
            wholecell_logger.info("Computing membrane currents")
        for idx, clamp in enumerate(self.parameters["Iinj"]):
            self.data["Imembrane_"+str(clamp)] = self.parameters["Cm"][idx]*(self.data[clamp].diff()/self.data["times"].diff())
            self.data.at[0, "Imembrane_"+str(clamp)] = 0.0
            self.data["Imembrane_"+str(clamp)] = self.data["Imembrane_"+str(clamp)] - self.data["Imembrane_"+str(clamp)][0]
        return self.data
    
    def filter_membrane_currents(self, log=False):
        if log:
            wholecell_logger.info("Filtering membrane currents")
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        for inj in self.parameters["Iinj"]:
            self.data["filtered_Imembrane_"+str(inj)] = self.filters["membrane_currents"].propagate(self.data["Imembrane_"+str(inj)], sampling_rate, log)
        return self.data
    
    def filter_activation_currents(self, log=False):
        if log:
            wholecell_logger.info("Filtering activation currents")
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        for inj in self.parameters["Iinj"]:
            activation_current = self.filters["activation_currents"].propagate(self.data["Iactivation_"+str(inj)], sampling_rate, log)
            self.data["filtered_Iactivation_"+str(inj)] = activation_current
        return self.data

    def squared_sum_of_negative_conductances(self, solution, log=False):
        if log:
            wholecell_logger.info("Computing objective function values")
        self.parameters["Eact"] = solution
        self.compute_activation_currents(log)
        self.filter_activation_currents(log)
        self.compute_passive_conductances(log)
        self.compute_stats(log)
        negative_going_excitation = self.data["excitation"][self.data["excitation"] < 0]
        negative_going_excitation = negative_going_excitation/np.amin(negative_going_excitation)
        negative_going_inhibition = self.data["inhibition"][self.data["inhibition"] < 0]
        negative_going_inhibition = negative_going_inhibition/np.amin(negative_going_inhibition)
        return np.square(np.sum(negative_going_excitation) + np.sum(negative_going_inhibition))
    
    def optimize(self, log=False):
        if log:
            wholecell_logger.info("Optimizing")
        self.estimate_conductances()
        def obj_wrapper(solution):
            return self.squared_sum_of_negative_conductances(solution)
        maximum_recorded_membrane_voltage = [max_val for max_val in self.data[[x for x in self.parameters["Iinj"]]].max()]
        steady_state_potential = np.asarray([x for x in self.parameters["Ess"]])
        activation_potential_bounds = Bounds(steady_state_potential, maximum_recorded_membrane_voltage)
        result = minimize(obj_wrapper, maximum_recorded_membrane_voltage, method='trust-constr', bounds=activation_potential_bounds)
        optimizer_logger.info(f"Optimization Success?{result.success} iterations:{result.niter}")
        self.parameters["Eact"] = result.x
        return result
    
    def compute_passive_conductances(self, log=False):
        if log:
            wholecell_logger.info("Estimating conductances")
        ntimesteps = self.data.shape[0]
        A = np.zeros((ntimesteps, 2, 2))
        B = np.zeros((ntimesteps, 2, 1))
        membrane_potential = self.data[[x for x in self.parameters["Iinj"]]].to_numpy()
        membrane_current = self.data[["filtered_Imembrane_"+str(x) for x in self.parameters["Iinj"]]].to_numpy()
        activation_current = self.data[["filtered_Iactivation_"+str(x) for x in self.parameters["Iinj"]]].to_numpy()
        leakage_current = self.data[["Ileakage_"+str(x) for x in self.parameters["Iinj"]]].to_numpy()
        excitatory_reversal_potential = self.parameters["Ee"].to_numpy()
        inhibitory_reversal_potentail = self.parameters["Ei"].to_numpy()
        injected_current = self.parameters["Iinj"].to_numpy()
        A[:, 0, 0] = np.sum(np.square(membrane_potential - excitatory_reversal_potential), axis=1)
        A[:, 0, 1] = np.sum((membrane_potential - excitatory_reversal_potential) * (membrane_potential- inhibitory_reversal_potentail), axis=1)
        A[:, 1, 0] = A[:, 0, 1]
        A[:, 1, 1] = np.sum(np.square(membrane_potential - inhibitory_reversal_potentail), axis=1)
        B[:, 0, 0] = -1.0*np.sum((membrane_current - activation_current - injected_current + leakage_current)*(membrane_potential - excitatory_reversal_potential), axis=1)
        B[:, 1, 0] = -1.0*np.sum((membrane_current - activation_current - injected_current + leakage_current)*(membrane_potential - inhibitory_reversal_potentail), axis=1)
        try:
            conductances = np.linalg.pinv(A) @ B
        except Exception as e:
            optimizer_logger.debug(f"pinv msg: {e}")
        self.data["excitation"] = conductances[:, 0, 0]
        self.data["inhibition"] = conductances[:, 1, 0]
        self.data["positive_excitation"] = self.data["excitation"]
        self.data["positive_inhibition"] = self.data["inhibition"]
        self.data.loc[self.data["positive_excitation"] < 0, "positive_excitation"] = 0.0
        self.data.loc[self.data["positive_inhibition"] < 0, "positive_inhibition"] = 0.0
        self.data["resultant_excitation"] = self.data["positive_excitation"] - self.data["positive_inhibition"]
        self.data["resultant_inhibition"] = self.data["positive_inhibition"] - self.data["positive_excitation"]
        self.data.loc[self.data["resultant_excitation"] < 0, "resultant_excitation"] = 0.0
        self.data.loc[self.data["resultant_inhibition"] < 0, "resultant_inhibition"] = 0.0
        return self.data
    
    def get_clamp_near_0(self, log=False):
        if log:
            wholecell_logger.info("computing the closest clamp to resting")
        index_of_minimum_injected_current = self.parameters["Iinj"].abs().argmin()
        minimum_injected_current = self.parameters["Iinj"][index_of_minimum_injected_current]
        return index_of_minimum_injected_current, minimum_injected_current
    
    def compute_stats(self, log=False):
        if log:
            wholecell_logger.info("Computing stats")
        index_of_minimum_injected_current, minimum_injected_current = self.get_clamp_near_0()
        stats = pd.DataFrame()
        paradigm_all_var_stats = self.data.mean(numeric_only=True).to_frame().T
        stats["depolarization"] = paradigm_all_var_stats[f"depolarization_{minimum_injected_current}"]
        stats["hyperpolarization"] = paradigm_all_var_stats[f"hyperpolarization_{minimum_injected_current}"]
        stats["Imembrane"] = paradigm_all_var_stats[f"filtered_Imembrane_{minimum_injected_current}"]
        stats["Ileakage"] = paradigm_all_var_stats[f"Ileakage_{minimum_injected_current}"]
        stats["Iactivation"] = paradigm_all_var_stats[f"filtered_Iactivation_{minimum_injected_current}"]
        stats["mean_excitation"] = paradigm_all_var_stats["positive_excitation"]
        stats["mean_inhibition"] = paradigm_all_var_stats["positive_inhibition"]
        stats["net_excitation"] = paradigm_all_var_stats["resultant_excitation"]
        stats["net_inhibition"] = paradigm_all_var_stats["resultant_inhibition"]
        stats["spikes_per_stimulus_repetition"] = self.parameters["sps"][index_of_minimum_injected_current]
        return stats

    def estimate_conductances(self, log=False):
        self.filter_membrane_potentials(log)
        self.compute_polarizations(log)
        self.compute_activation_currents(log)
        self.filter_activation_currents(log)
        self.compute_leakage_currents(log)
        self.compute_membrane_currents(log)
        self.filter_membrane_currents(log)
        self.compute_passive_conductances(log)
        self.stats = self.compute_stats(log)
        return self.data

class Analyzer:
    def __init__(self, filepaths: Path, output_path: Path, filter_configurations: dict):
        self.fileaddress_handler = FileAddressHandler()
        self.filepaths = self.fileaddress_handler.split_fileaddress_list(filepaths)
        self.output_path = output_path
        self.filter_configurations = filter_configurations
        pass

    def get_paradigm_to_optimize(self, recordings):
        analysis_logger.info("Finding the best paradigm to optimize")
        maximum_depolarizations = np.zeros(len(recordings))
        for idx, paradigm in enumerate(recordings):
            index_of_minimum_injected_current, minimum_injected_current = recordings[paradigm].get_clamp_near_0()
            membrane_potential_at_minimum_injected_current = recordings[paradigm].data[minimum_injected_current]
            depols = membrane_potential_at_minimum_injected_current - recordings[paradigm].parameters["Ess"][index_of_minimum_injected_current]
            maximum_depolarizations[idx] = np.max(depols)
        max_index = np.argmax(maximum_depolarizations)
        return list(recordings.keys())[max_index]
    
    def estimate_optimum_activation_potential(self, recordings):
        optim_paradigm = self.get_paradigm_to_optimize(recordings)
        analysis_logger.info(f"Optimizing activation potentials for current clamps in {optim_paradigm}")
        try:
            result = recordings[optim_paradigm].optimize()
        except Exception as e:
            analysis_logger.debug(f"{e}")
        dEact = recordings[optim_paradigm].parameters["Eact"] - recordings[optim_paradigm].parameters["Ess"]
        recordings[optim_paradigm].stats.insert(0, "paradigm", optim_paradigm)
        analysis_logger.info(f"Optimization results: convergence = {result.success}, iterations = {result.niter}")
        analysis_logger.info(f"Optimum activation potentials: {recordings[optim_paradigm].parameters["Eact"].to_numpy().tolist()}")
        overall_stats = recordings[optim_paradigm].stats.copy()
        for idx, paradigm in enumerate(recordings):
            if not paradigm == optim_paradigm:
                recordings[paradigm].parameters["Eact"] = recordings[paradigm].parameters["Ess"] + dEact
                recordings[paradigm].estimate_conductances(log=False)
                recordings[paradigm].stats.insert(0, "paradigm", paradigm)
                overall_stats = pd.concat([overall_stats, recordings[paradigm].stats], axis=0)
            analysis_logger.info(f"estimated conductances for {paradigm}")
        overall_stats = overall_stats.sort_values(by="paradigm", ascending=True)
        return recordings, overall_stats
    
    def estimate_optimum_activation_potential_each_paradigm(self, recordings):
        for idx, paradigm in enumerate(recordings):
            print(f"Optimizing Eact for {paradigm}")
            recordings[paradigm].optimize()
            recordings[paradigm].stats.insert(0, "paradigm", paradigm)
            if idx == 0:
                overall_stats = recordings[paradigm].stats.copy()
            else:
                overall_stats = pd.concat([overall_stats, recordings[paradigm].stats], axis=0)
            print(f"{paradigm} done")
        overall_stats = overall_stats.sort_values(by="paradigm", ascending=True)
        return recordings, overall_stats
    
    def estimation_without_optim_activation_potential(self, recordings):
        for idx, paradigm in enumerate(recordings):
            recordings[paradigm].estimate_conductances()
            recordings[paradigm].stats.insert(0, "paradigm", paradigm)
            if idx == 0:
                overall_stats = recordings[paradigm].stats.copy()
            else:
                overall_stats = pd.concat([overall_stats, recordings[paradigm].stats], axis=0)
            print(f"{paradigm} done")
        overall_stats = overall_stats.sort_values(by="paradigm", ascending=True)
        return recordings, overall_stats

    def analyze(self, fileaddress, optimize=1):
        filepath = fileaddress[0]
        filename = fileaddress[1]
        fileformat = fileaddress[2]
        analysis_logger.info(f"Reading: {filename}")
        reader = XLReader(Path(filepath) / (filename + fileformat))
        recordings = {}
        for _, paradigm in enumerate(reader.get_paradigms()):
            recordings[paradigm] = WholeCellRecording(
                reader.get_paradigm_data(paradigm), 
                reader.get_paradigm_parameters(paradigm),
                self.filter_configurations
            )
        if optimize == 0:
            analysis_logger.info(f"Level 0 optimization: No optimizing using user provided values.")
            recordings, overall_stats = self.estimation_without_optim_activation_potential(recordings)
            self.logger.info(f"Level 0 optimization: Completed.")
        elif optimize == 1:
            analysis_logger.info(f"Level 1 optimization: Activation potentials optimized using paradigm with maximum depolarization.")
            recordings, overall_stats = self.estimate_optimum_activation_potential(recordings)
            analysis_logger.info(f"Level 1 optimization: Completed.")
        elif optimize == 2:
            analysis_logger.info(f"Level 3 optimization: Activation potentials optimized for every paradigm.")
            recordings, overall_stats = self.estimate_optimum_activation_potential_each_paradigm(recordings)
            analysis_logger.info(f"Level 3 optimization: Complete.")
        result_filename = self.output_path / f"{filename}_analyzed"
        analysis_logger.info(f"Analysis of {filename} completed.")
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
        analysis_logger.info(f"Verbose plotting of conductance estimations for {self.filepaths[1]}")
        fig, axs = plt.subplots(nrows = 7, ncols = len(recordings), sharex="all", sharey="row", figsize=(15, 10), constrained_layout=True)
        for idx, paradigm in enumerate(recordings):
            if "representative" in recordings[paradigm].data:
                rep = recordings[paradigm].data["representative"].to_numpy()
            membrane_potential = recordings[paradigm].data[[x for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            membrane_current = recordings[paradigm].data[["filtered_Imembrane_"+str(x) for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            leakage_current = recordings[paradigm].data[["Ileakage_"+str(x) for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            activation_current = recordings[paradigm].data[["filtered_Iactivation_"+str(x) for x in recordings[paradigm].parameters["Iinj"]]].to_numpy()
            conductances = recordings[paradigm].data[["excitation", "inhibition"]].to_numpy()
            times = recordings[paradigm].data["times"].to_numpy()
            resting_potential = times*0 + recordings[paradigm].parameters["Er"][0]
            threshold_potential = times*0 + recordings[paradigm].parameters["Et"][0]
            activation_potential = membrane_potential*0 + recordings[paradigm].parameters["Eact"].to_numpy()
            if "stimulus" in recordings[paradigm].data:
                stim = recordings[paradigm].data["stimulus"].to_numpy()
            axs[0, idx].set_title(paradigm)
            if "representative" in recordings[paradigm].data:
                axs[0, idx].plot(times, rep)
                axs[0, idx].plot(times, resting_potential, '--k')
            axs[0, idx].set_ylabel("Rep. Vm (V)")
            axs[0, idx].set_title(paradigm)
            axs[0, idx].grid(True)
            curves = axs[1, idx].plot(times, membrane_potential)
            colors = [x.get_color() for x in curves]
            axs[1, idx].plot(times, resting_potential, '--k')
            axs[1, idx].plot(times, threshold_potential, linestyle='--', color=(0.5, 0.5, 0.5))
            for i in range(activation_potential.shape[-1]):
                axs[1, idx].plot(times, activation_potential[:, i], linestyle='--', color=colors[i])
            axs[1, idx].set_ylabel("Vm (V)")
            axs[1, idx].grid(True)
            axs[2, idx].plot(times, membrane_current)
            axs[2, idx].set_ylabel("Im (A)")
            axs[2, idx].grid(True)
            axs[3, idx].plot(times, leakage_current)
            axs[3, idx].set_ylabel("Ileak (A)")
            axs[3, idx].grid(True)
            axs[4, idx].plot(times, activation_current)
            axs[4, idx].set_ylabel("Iact (A)")
            axs[4, idx].grid(True)
            # axs[4, idx].set_title(", ".join([f"{val:.3f}" for val in recordings[paradigm].parameters["Eact"]]))
            axs[5, idx].plot(times, conductances[:, 0], c='r')
            axs[5, idx].plot(times, conductances[:, 1], c='b')
            axs[5, idx].plot(times, times*0, '--k')
            axs[5, idx].set_ylabel("G (S)")
            axs[5, idx].grid(True)
            if "stimulus" in recordings[paradigm].data:
                axs[6, idx].plot(times, stim)
            axs[6, idx].set_ylabel("Stimulus")
            axs[6, idx].set_xlabel("times(sec)")
            axs[6, idx].grid(True)
        analysis_logger.info(f"Verbose plotting of conductance estimations for {self.filepaths[1]} completed.")
        analysis_logger.info(f"Saving verbose plotting of conductance estimations for {self.filepaths[1]}")
        plt.savefig(str(filename)+f"_dev_traces.png")
        analysis_logger.info(f"Saving verbose plotting of conductance estimations for {self.filepaths[1]} completed.")
        pass

    def set_stats_scale(self, ax, scale_max, margin=0.1):
        scale_max = scale_max + margin*scale_max
        try:
            ax.set_ylim([-1*scale_max, scale_max])
        except Exception as e:
            plotter_logger.debug(f"{e}")
        pass

    def plot_stats_dev(self, recordings, filename: Path):
        analysis_logger.info(f"Verbose plotting of stats for {self.filepaths[1]}")
        fig, axs = plt.subplots(nrows = 5, ncols = 1, sharex="all", figsize=(15, 10), constrained_layout=True)
        mean_depolarizations = np.asarray([recordings[paradigm].stats["depolarization"][0] for paradigm in recordings])
        mean_hyperpolarizations = np.asarray([recordings[paradigm].stats["hyperpolarization"][0] for paradigm in recordings])
        mean_Im = np.asarray([recordings[paradigm].stats["Imembrane"][0] for paradigm in recordings])
        mean_Ileak = np.asarray([recordings[paradigm].stats["Ileakage"][0] for paradigm in recordings])
        mean_Iactivation = np.asarray([recordings[paradigm].stats["Iactivation"][0] for paradigm in recordings])
        mean_excitation = np.asarray([recordings[paradigm].stats["mean_excitation"][0] for paradigm in recordings])
        mean_inhibition = np.asarray([recordings[paradigm].stats["mean_inhibition"][0] for paradigm in recordings])
        net_excitation = np.asarray([recordings[paradigm].stats["net_excitation"][0] for paradigm in recordings])
        net_inhibition = np.asarray([recordings[paradigm].stats["net_inhibition"][0] for paradigm in recordings])
        spikes_per_stimulus_repetition = np.asarray([recordings[paradigm].stats["spikes_per_stimulus_repetition"][0] for paradigm in recordings])
        paradigms = [paradigm for paradigm in recordings]
        xlocations = np.asarray([x for x in range(len(paradigms))])
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
        ax.plot(xlocations, spikes_per_stimulus_repetition, color='k', marker = 'o')
        scale_max = np.amax(spikes_per_stimulus_repetition)
        self.set_stats_scale(ax, scale_max, 0.1)
        axs[4].set_ylabel("G (S)")
        analysis_logger.info(f"Verbose plotting of stats for {self.filepaths[1]} completed.")
        analysis_logger.info(f"Saving verbose plotting of stats for {self.filepaths[1]}")
        plt.savefig(str(filename)+f"_dev_stats.png")
        analysis_logger.info(f"Saving verbose plotting of stats for {self.filepaths[1]} completed.")
        pass

    def run(self, optimization_level=1):
        for filepath in self.filepaths:
            recordings, stats, result_filename = self.analyze(filepath, optimize=optimization_level)
            # self.write_to_excel(f"{result_filename}.xlsx", recordings, stats)
            self.plot_dev(recordings, result_filename)
            self.plot_stats_dev(recordings, result_filename)
        pass