from pathlib import Path
from pyhelpers.store import save_fig

from scipy.optimize import nnls

from libs.readers import XLReader, FilterCfg, AnalyzerCfg
from libs.utils import *

from typing import Dict, Optional, List, Tuple

class LowPassFilter:
    def __init__(self, passband: float, stopband: float, attenuation: float, ripple: float, name: str="generic") -> None:
        assert 0 < passband < stopband
        assert attenuation >= 0
        assert ripple >= 0

        self.passband: float = passband
        self.stopband: float = stopband
        self.attenuation: float = attenuation
        self.ripple: float = ripple
        self.name: str = name

    def compute_minimum_order(self, fs: float, log: bool=False) -> Tuple[int, float]:
        assert fs > 0

        nyquist_freq: float = fs / 2
        normalized_passband = self.passband / nyquist_freq
        normalized_stopband = self.stopband / nyquist_freq
        order, normalized_cutoff_frequency = filters.buttord(normalized_passband, normalized_stopband, self.ripple, self.attenuation)
        assert isinstance(normalized_cutoff_frequency, float)

        if log: 
            lpf_logger.info(f"{self.name} low pass filter ")
            lpf_logger.info(f"{self.name} low pass filter computed minimum order: {order} with transition gap: {self.stopband - self.passband}")

        return order, normalized_cutoff_frequency
    
    def append_samples(self, signal: np.ndarray) -> Tuple[np.ndarray, int]:
        nsamples: int = signal.shape[0]
        nappend: int = int(nsamples / 2)
        append_samples = np.zeros((nappend,), dtype=signal.dtype)+signal[0]
        return np.concatenate([append_samples, signal, append_samples], axis=0), nappend
    
    def deppend_samples(self, signal: np.ndarray, nappend: int) -> np.ndarray:
        return signal[nappend:-nappend, ...]
    
    def propagate(self, input: np.ndarray, sampling_rate: float, log=False):
        order, normalized_frequency = self.compute_minimum_order(sampling_rate, log)
        second_order_sections = filters.butter(order, normalized_frequency, output="sos")
        input_signal, nappend = self.append_samples(input)
        output_signal = filters.sosfiltfilt(second_order_sections, input_signal)
        return self.deppend_samples(output_signal, nappend)

class WholeCellRecording:
    def __init__(self, data: pd.DataFrame, parameters: pd.DataFrame, filter_parameters: Dict[str, 'FilterCfg'], current_clamps: Optional[List[float]]=None) -> None:
        """
        Args:
            data: pd.DataFrame
            parameters: pd.DataFrame
            filter_parameters: Dict[str, FilterCfg]
            current_clamps: Optional[List[float]]
        Returns:
            None
        """
        self.filters: Dict[str, LowPassFilter] = {key:LowPassFilter(**params.__dict__, name=key) for key, params in filter_parameters.items()}
        self.complier = Compliance()

        self.data: pd.DataFrame = data

        self.parameters: pd.DataFrame = parameters
        if current_clamps is not None:        
            if set(current_clamps) <= set(parameters["Iinj"]):
                self.parameters = parameters[parameters["Iinj"].isin(current_clamps)]
            else:
                wholecell_logger.warning("Invalid current_clamps argument in WholeCellRecording initialization: not a subset of clamps in excel file. Using all clamps present in file.")
        
        self.scale_data()

    def check_compliance(self):
        if self.complier.check_compliance(self.data, self.parameters):
            return True
        else:
            wholecell_logger.debug(f"The data / parameters do not comply with standards.")
        return False

    def scale_data(self, log=False):
        if log:
            wholecell_logger.info("Scaling membrane voltage")
        for idx, clamp in zip(self.parameters["Iinj"].keys(), self.parameters["Iinj"]):
            self.data[f"{clamp:.3e}"] = ((self.data[f"{clamp:.3e}"] - self.data[f"{clamp:.3e}"][0])*1e-3)+self.parameters["Ess"][idx]
        return self.data
    
    def filter_membrane_potentials(self, log=False):
        if log:
            wholecell_logger.info("Filtering membrane potentials")
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        for inj in self.parameters["Iinj"]:
            self.data[f"{inj:.3e}"] = self.filters["membrane_potentials"].propagate(self.data[f"{inj:.3e}"].to_numpy(), sampling_rate, log)
        return self.data
    
    def compute_activation_conductance_constants(self, log=False):
        if log:
            wholecell_logger.info("Computing activation constants alpha and beta")
        self.parameters["alpha"] = (1.0/self.parameters["Rin"])/(2.0*(self.parameters["Eact"] - self.parameters["Ess"]))
        self.parameters["beta"] = self.parameters["alpha"]*(self.parameters["Et"] - self.parameters["Ess"])
        self.parameters["alpha"] = self.parameters["alpha"]*self.parameters["xalpha"]
        self.parameters["beta"] = self.parameters["beta"]*self.parameters["xbeta"]
        if log:
            wholecell_logger.info(f"alpha = {self.parameters['alpha'].to_numpy().tolist()} & beta = {self.parameters['beta'].to_numpy().tolist()}")
        return self.parameters
    
    def compute_polarizations(self, log=False):
        if log:
            wholecell_logger.info("Computing polarizations")
        for idx, clamp in zip(self.parameters["Iinj"].keys(), self.parameters["Iinj"]):
            self.data[f"depolarization_{clamp:.3e}"] = np.where(self.data[f"{clamp:.3e}"] > self.parameters["Ess"][idx], self.data[f"{clamp:.3e}"] - self.parameters["Ess"][idx], 0)
            self.data[f"hyperpolarization_{clamp:.3e}"] = np.where(self.data[f"{clamp:.3e}"] < self.parameters["Ess"][idx], self.data[f"{clamp:.3e}"] - self.parameters["Ess"][idx], 0)
        return self.data
    
    def compute_leakage_currents(self, log=False):
        if log:
            wholecell_logger.info("Computing leakage currents")
        for idx, clamp in zip(self.parameters["Iinj"].keys(), self.parameters["Iinj"]):
            self.data[f"Ileakage_{clamp:.3e}"] = (1/self.parameters["Rin"][idx])*(self.data[f"{clamp:.3e}"] - self.parameters["Er"][idx])
        return self.data
    
    def compute_activation_currents(self, log=False):
        if log:
            wholecell_logger.info("Computing activation currents")
        self.compute_activation_conductance_constants(log)
        for idx, clamp in zip(self.parameters["Iinj"].keys(), self.parameters["Iinj"]):
            alpha_current = self.parameters["alpha"][idx]*(self.data[f"{clamp:.3e}"] - self.parameters["Ess"][idx])*(self.data[f"{clamp:.3e}"] - self.parameters["Et"][idx])
            beta_current = self.parameters["beta"][idx]*(self.data[f"{clamp:.3e}"] - self.parameters["Ess"][idx])
            activation_current = alpha_current + beta_current
            activation_current[self.data[f"{clamp:.3e}"] < self.parameters["Ess"][idx]] = 0.0
            activation_current[self.data[f"{clamp:.3e}"] > self.parameters["Et"][idx]] = 0.0
            self.data[f"Iactivation_{clamp:.3e}"] = activation_current
        return self.data
    
    def compute_membrane_currents(self, log=False):
        if log:
            wholecell_logger.info("Computing membrane currents")
        for idx, clamp in zip(self.parameters["Iinj"].keys(), self.parameters["Iinj"]):
            self.data[f"Imembrane_{clamp:.3e}"] = self.parameters["Cm"][idx]*(self.data[f"{clamp:.3e}"].diff()/self.data["times"].diff())
            self.data.at[0, f"Imembrane_{clamp:.3e}"] = 0.0
            self.data[f"Imembrane_{clamp:.3e}"] = self.data[f"Imembrane_{clamp:.3e}"] - self.data[f"Imembrane_{clamp:.3e}"][0]
        return self.data
    
    def filter_membrane_currents(self, log=False):
        if log:
            wholecell_logger.info("Filtering membrane currents")
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        for inj in list(self.parameters["Iinj"]):
            self.data[f"filtered_Imembrane_{inj:.3e}"] = self.filters["membrane_currents"].propagate(self.data[f"Imembrane_{inj:.3e}"], sampling_rate, log)
        return self.data
    
    def filter_activation_currents(self, log=False):
        if log:
            wholecell_logger.info("Filtering activation currents")
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        for inj in list(self.parameters["Iinj"]):
            activation_current = self.filters["activation_currents"].propagate(self.data[f"Iactivation_{inj:.3e}"], sampling_rate, log)
            self.data[f"filtered_Iactivation_{inj:.3e}"] = activation_current 
        return self.data
    
    def compute_passive_conductances(self, bin_ms: float=5.0, log: bool=False, include_beta: bool=False):
        """
        Level-1 conductance estimation:
        - Voltage-domain via bin-integrated membrane equation (no dv/dt)
        - Piecewise-constant (boxcar) ge, gi per bin
        - Nonnegativity via NNLS
        - No extra regularizers/priors

        Writes per-timestep ge/gi by expanding bin estimates back to the full timeline.
        """

        if log:
            wholecell_logger.info("Estimating conductances (Level 1: integrated + NNLS)")

        # --- pull parameters ---
        Iinj_levels = list(self.parameters["Iinj"])
        colnames_v = [f"{x:.3e}" for x in Iinj_levels]

        v = self.data[colnames_v].to_numpy()  # shape: (T, N)
        T, N = v.shape

        dt = self.data["times"][1] - self.data["times"][0]

        # constants
        C = float(self.parameters["Cm"][0])         # Farads
        gl = 1 / float(self.parameters["Rin"][0])   # Siemens
        Er = float(self.parameters["Er"][0])        # Volts
        Ee = float(self.parameters["Ee"][0])        # Volts
        Ei = float(self.parameters["Ei"][0])        # Volts

        beta = float(self.parameters["alpha"][0])
        Et = float(self.parameters["Et"][0])

        # injected currents per clamp (assume constant in time for now)
        Iinj = np.asarray(self.parameters["Iinj"].to_numpy(), dtype=float)  # shape: (N,)

        # binning
        bin_s = bin_ms / 1000.0
        bin_len = max(1, int(round(bin_s / dt)))
        nbins = int(np.ceil(T / bin_len))
        print(bin_s, bin_len, nbins)

        ge_bins = np.zeros(nbins)
        gi_bins = np.zeros(nbins)
        resnorm_bins = np.full(nbins, np.nan)
        cond_bins = np.full(nbins, np.nan)

        # precompute terms that get integrated
        # leak integrand: gl*(Er - v)
        leak = gl * (Er - v)  # (T, N)

        # injected current integrand: Iinj (constant over time)
        inj = np.tile(Iinj.reshape(1, -1), (T, 1))  # (T, N)

        if include_beta and beta != 0.0:
            # active current term: -beta(Er - v)(Et - v) appears on RHS in your original.
            # In our rearranged y we add +beta * integral((Er - v)(Et - v)) if using the form in earlier math.
            quad = (Er - v) * (Et - v)  # (T, N)
        else:
            quad = None

        # helper for trapezoid integral over [a:b] for each clamp column
        def trapz_segment(arr, a, b):
            # arr shape (T, N)
            # integrate per column over indices [a, b) using trapezoid rule
            seg = arr[a:b, :]
            if seg.shape[0] < 2:
                # fall back: rectangle
                return seg.sum(axis=0) * dt
            return np.trapz(seg, dx=dt, axis=0)

        # main loop: build y_k and X_k and solve NNLS
        for k in range(nbins):
            a = k * bin_len
            b = min((k + 1) * bin_len, T)
            if b - a < 1:
                continue

            # y_i,k = C*(v_i(b)-v_i(a)) - ∫[ gl(Er-v) + Iinj ] dt  (+ beta ∫ quad dt if include_beta)
            dv = v[b - 1, :] - v[a, :]  # (N,)
            yk = C * dv - (trapz_segment(leak + inj, a, b))

            if quad is not None:
                yk = yk + beta * trapz_segment(quad, a, b)

            # X columns: Ae_i,k = ∫(Ee - v_i)dt ; Ai_i,k = ∫(Ei - v_i)dt
            Ae = trapz_segment((Ee - v), a, b)  # (N,)
            Ai = trapz_segment((Ei - v), a, b)  # (N,)

            Xk: np.ndarray = np.column_stack([Ae, Ai])  # (N, 2)

            # conditioning diagnostic (optional but highly recommended)
            try:
                XtX = Xk.T @ Xk
                cond_bins[k] = np.linalg.cond(XtX)
            except Exception:
                cond_bins[k] = np.nan

            # NNLS solve (nonnegative ge, gi)
            # Note: if Xk is near-rank-deficient, NNLS will still return *a* solution;
            # that's why cond_bins matters.
            gk, rnorm = nnls(Xk, yk)
            ge_bins[k], gi_bins[k] = gk
            resnorm_bins[k] = rnorm

        # expand bins back to timesteps
        ge = np.repeat(ge_bins, bin_len)[:T]
        gi = np.repeat(gi_bins, bin_len)[:T]

        self.data["excitation"] = ge
        self.data["inhibition"] = gi

        # (Optional) store diagnostics
        self.data["bin_index"] = np.repeat(np.arange(nbins), bin_len)[:T]
        self.data["bin_resnorm"] = np.repeat(resnorm_bins, bin_len)[:T]
        self.data["bin_cond_XtX"] = np.repeat(cond_bins, bin_len)[:T]

        # remove all the post-hoc clipping logic; NNLS already enforces positivity
        # If you still want resultant E/I "dominance" signals, compute them cleanly:
        self.data["resultant_excitation"] = (self.data["excitation"] - self.data["inhibition"]).clip(lower=0.0)
        self.data["resultant_inhibition"] = (self.data["inhibition"] - self.data["excitation"]).clip(lower=0.0)

        return self.data
        
    def get_clamp_near_0(self, log=False) -> Tuple[int, float]:
        if log:
            wholecell_logger.info("computing the closest clamp to resting")
        index_of_minimum_injected_current: np.intp = np.argmin(np.abs(self.parameters["Iinj"]))
        minimum_injected_current: float = self.parameters["Iinj"][index_of_minimum_injected_current]
        return index_of_minimum_injected_current, minimum_injected_current
    
    def compute_stats(self, log=False):
        if log:
            wholecell_logger.info("Computing stats")
        index_of_minimum_injected_current, minimum_injected_current = self.get_clamp_near_0()
        stats = pd.DataFrame()
        paradigm_all_var_stats = self.data.mean(numeric_only=True).to_frame().T
        stats["depolarization"] = paradigm_all_var_stats[f"depolarization_{minimum_injected_current:.3e}"]
        stats["hyperpolarization"] = paradigm_all_var_stats[f"hyperpolarization_{minimum_injected_current:.3e}"]
        stats["Imembrane"] = paradigm_all_var_stats[f"filtered_Imembrane_{minimum_injected_current:.3e}"]
        stats["Ileakage"] = paradigm_all_var_stats[f"Ileakage_{minimum_injected_current:.3e}"]
        stats["Iactivation"] = paradigm_all_var_stats[f"filtered_Iactivation_{minimum_injected_current:.3e}"]
        stats["mean_excitation"] = paradigm_all_var_stats["excitation"]
        stats["mean_inhibition"] = paradigm_all_var_stats["inhibition"]
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
        self.compute_passive_conductances(log=log)
        self.stats = self.compute_stats(log)
        return self.data
    
    def squared_sum_of_negative_conductances(self, solution, log=False):
        if log:
            wholecell_logger.info("Computing objective function values")
        self.parameters["Eact"] = solution
        self.compute_activation_currents(log)
        self.filter_activation_currents(log)
        self.compute_passive_conductances(log=log)
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
        maximum_recorded_membrane_voltage = [max_val for max_val in self.data[list(self.parameters["Iinj"])].max()]
        steady_state_potential = np.asarray([x for x in self.parameters["Ess"]])
        activation_potential_bounds = Bounds(steady_state_potential, maximum_recorded_membrane_voltage)
        result = minimize(self.squared_sum_of_negative_conductances, maximum_recorded_membrane_voltage, method='trust-constr', bounds=activation_potential_bounds)
        optimizer_logger.info(f"Optimization Success?{result.success} iterations:{result.niter}")
        self.parameters["Eact"] = result.x
        return result

class Analyzer:
    def __init__(self, cfg: AnalyzerCfg):
        self.cfg: AnalyzerCfg = cfg

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
        analysis_logger.info(f"Optimum activation potentials: {recordings[optim_paradigm].parameters['Eact'].to_numpy().tolist()}")
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
    
    def estimation_without_optim_activation_potential(self, recordings: Dict[str, WholeCellRecording]):
        assert len(recordings) > 0

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

    def analyze(self, filepath: str, filter_configurations: Dict[str, 'FilterCfg'], optimize: int=0, current_clamps: Optional[List[float]]=None):
        basename: str = os.path.basename(filepath)
        analysis_logger.info(f"Reading: {basename}")
        reader = XLReader(filepath)
        recordings: Dict[str, WholeCellRecording] = {}
        for _, paradigm in enumerate(reader.get_paradigms()):
            recordings[paradigm] = WholeCellRecording(
                reader.get_paradigm_data(paradigm), 
                reader.get_paradigm_parameters(paradigm),
                filter_parameters=filter_configurations,
                current_clamps=current_clamps
            )
        if optimize == 0:
            analysis_logger.info(f"Level 0 optimization: No optimization--using user provided values.")
            recordings, overall_stats = self.estimation_without_optim_activation_potential(recordings)
            # self.logger.info(f"Level 0 optimization: Completed.")
        elif optimize == 1:
            analysis_logger.info(f"Level 1 optimization: Activation potentials optimized using paradigm with maximum depolarization.")
            recordings, overall_stats = self.estimate_optimum_activation_potential(recordings)
            analysis_logger.info(f"Level 1 optimization: Completed.")
        elif optimize == 2:
            analysis_logger.info(f"Level 2 optimization: Activation potentials optimized for every paradigm.")
            recordings, overall_stats = self.estimate_optimum_activation_potential_each_paradigm(recordings)
            analysis_logger.info(f"Level 2 optimization: Complete.")
        result_filename = os.path.join(self.cfg.image_save_dir, f"{os.path.splitext(basename)[0]}_analyzed")
        analysis_logger.info(f"Analysis of {basename} completed.")
        return recordings, overall_stats, result_filename

    def plot_dev(self, recordings, filename: Path, filetype: str="png", current_clamps: Optional[List[float]]=None):
        analysis_logger.info(f"Verbose plotting of conductance estimations for {filename}")
        fig, axs = plt.subplots(nrows = 8, ncols = len(recordings), sharex="all", sharey="row", figsize=(15, 10), constrained_layout=True)
        for idx, paradigm in enumerate(recordings):
            paradigm_iinj: List[float] = list(recordings[paradigm].parameters["Iinj"])
            if current_clamps is not None and set(current_clamps) <= set(paradigm_iinj):
                paradigm_iinj = list(set(paradigm_iinj).intersection(current_clamps))
                assert len(paradigm_iinj) > 0, "Cannot plot. Specified current clamps have no intersection with current clamps listed in parameters."

            if "representative" in recordings[paradigm].data:
                rep = recordings[paradigm].data["representative"].to_numpy()
            membrane_potential = recordings[paradigm].data[[f"{x:.3e}" for x in paradigm_iinj]].to_numpy()

            # predicted_membrane_potential = recordings[paradigm].data["predicted_membrane_potential"]

            membrane_current = recordings[paradigm].data[[f"filtered_Imembrane_{x:.3e}" for x in paradigm_iinj]].to_numpy()
            leakage_current = recordings[paradigm].data[[f"Ileakage_{x:.3e}" for x in paradigm_iinj]].to_numpy()
            activation_current = recordings[paradigm].data[[f"filtered_Iactivation_{x:.3e}" for x in paradigm_iinj]].to_numpy()
            conductances = recordings[paradigm].data[["excitation", "inhibition"]].to_numpy()
            # errors = recordings[paradigm].data["best_fit_squared_error"].to_numpy()
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
            # axs[1, idx].plot(times, predicted_membrane_potential, '--k')
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
            
            # axs[6, idx].plot(times, errors, '--k')
            axs[6, idx].set_ylabel("Err")
            axs[6, idx].grid(True)

            if "stimulus" in recordings[paradigm].data:
                axs[6, idx].plot(times, stim)
            axs[7, idx].set_ylabel("Stimulus")
            axs[7, idx].set_xlabel("times(sec)")
            axs[7, idx].grid(True)
        analysis_logger.info(f"Saving verbose plotting of conductance estimations for {filename}")
        filename: str = f"{str(filename)}_dev_traces."
        if filetype == "png":
            plt.savefig(f"{filename}{filetype}")
        elif filetype == "emf":
            save_fig(f"{filename}svg", dpi=300, conv_svg_to_emf=True, verbose=True)
        else:
            raise ValueError(f"The filetype requested ({filetype}) is not yet implemented :) Please consult James or Rishi.")
        plt.show()

    def set_stats_scale(self, ax, scale_max, margin=0.1):
        scale_max = scale_max + margin*scale_max
        try:
            ax.set_ylim([-1*scale_max, scale_max])
        except Exception as e:
            plotter_logger.debug(f"{e}")
        pass

    def plot_stats_dev(self, recordings, filename: Path):
        analysis_logger.info(f"Verbose plotting of stats for {filename}")
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
        analysis_logger.info(f"Saving verbose plotting of stats for {filename}")
        plt.savefig(str(filename)+f"_dev_stats.png")
        pass

    def run(self) -> None:
        for i in range(len(self.cfg.paths_to_all_analysis_spreadsheets)):
            iinj_clamps_to_use: Optional[List[List[float]]] = self.cfg.iinj_clamps_to_use if self.cfg.iinj_clamps_to_use is None else self.cfg.iinj_clamps_to_use[i]
            recordings, stats, result_filename = self.analyze(
                self.cfg.paths_to_all_analysis_spreadsheets[i], 
                optimize=self.cfg.optimization_level, 
                current_clamps=iinj_clamps_to_use,
                filter_configurations=self.cfg.timeseries_filter_cfgs
            )
            # self.write_to_excel(f"{result_filename}.xlsx", recordings, stats)
            self.plot_dev(recordings, result_filename, filetype=self.cfg.image_save_type, current_clamps=iinj_clamps_to_use)
            self.plot_stats_dev(recordings, result_filename)