import pandas as pd
import numpy as np

from scipy.signal import butter, buttord, sosfiltfilt

from pathlib import Path
from pyhelpers.store import save_fig

from libs.readers import XLReader

from typing import Dict, Optional, List, Tuple


class LowPassFilter:
    def __init__(self, passband: float, stopband: float, attenuation: float, ripple: float, fs: float, name: str="generic") -> None:
        assert passband < stopband, f"For {name} low pass filter stopband cannot be less than passband."
    
        self.name: str = name
        
        self.filter_design = butter(*buttord(passband, stopband, ripple, attenuation, fs=fs), output="sos", fs=fs)

    def propagate(self, raw_signal: np.ndarray) -> np.ndarray:
        return sosfiltfilt(self.filter_design, raw_signal)

class WholeCellRecording:
    def __init__(self, data: pd.DataFrame, parameters: pd.DataFrame, filters: Dict[str, LowPassFilter], current_clamps: Optional[List[float]]=None) -> None:
        self.filters: Dict[str, LowPassFilter] = filters

        self.data: pd.DataFrame = data

        self.parameters: pd.DataFrame = parameters

    def filter_membrane_potentials(self):
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        for inj in self.parameters["Iinj"]:
            self.data[f"{inj:.3e}"] = self.filters["membrane_potentials"].propagate(self.data[f"{inj:.3e}"], sampling_rate)
        return self.data
    
    def compute_activation_conductance_constants(self):
        self.parameters["alpha"] = (1.0/self.parameters["Rin"])/(2.0*(self.parameters["Eact"] - self.parameters["Ess"]))
        self.parameters["beta"] = self.parameters["alpha"]*(self.parameters["Et"] - self.parameters["Ess"])
        self.parameters["alpha"] = self.parameters["alpha"]*self.parameters["xalpha"]
        self.parameters["beta"] = self.parameters["beta"]*self.parameters["xbeta"]
        return self.parameters
    
    def compute_polarizations(self):
        for idx, clamp in zip(self.parameters["Iinj"].keys(), self.parameters["Iinj"]):
            self.data["depolarization_"+str(clamp)] = np.where(self.data[f"{clamp:.3e}"] > self.parameters["Ess"][idx], self.data[f"{clamp:.3e}"] - self.parameters["Ess"][idx], 0)
            self.data["hyperpolarization_"+str(clamp)] = np.where(self.data[f"{clamp:.3e}"] < self.parameters["Ess"][idx], self.data[f"{clamp:.3e}"] - self.parameters["Ess"][idx], 0)
        return self.data
    
    def compute_leakage_currents(self):
        for idx, clamp in zip(self.parameters["Iinj"].keys(), self.parameters["Iinj"]):
            self.data["Ileakage_"+str(clamp)] = (1/self.parameters["Rin"][idx])*(self.data[f"{clamp:.3e}"] - self.parameters["Er"][idx])
        return self.data
    
    def compute_activation_currents(self):
        self.compute_activation_conductance_constants(log)
        for idx, clamp in zip(self.parameters["Iinj"].keys(), self.parameters["Iinj"]):
            alpha_current = self.parameters["alpha"][idx]*(self.data[f"{clamp:.3e}"] - self.parameters["Ess"][idx])*(self.data[f"{clamp:.3e}"] - self.parameters["Et"][idx])
            beta_current = self.parameters["beta"][idx]*(self.data[f"{clamp:.3e}"] - self.parameters["Ess"][idx])
            activation_current = alpha_current + beta_current
            activation_current[self.data[f"{clamp:.3e}"] < self.parameters["Ess"][idx]] = 0.0
            activation_current[self.data[f"{clamp:.3e}"] > self.parameters["Et"][idx]] = 0.0
            self.data["Iactivation_"+str(clamp)] = activation_current
        return self.data
    
    def compute_membrane_currents(self):
        for idx, clamp in zip(self.parameters["Iinj"].keys(), self.parameters["Iinj"]):
            self.data["Imembrane_"+str(clamp)] = self.parameters["Cm"][idx]*(self.data[f"{clamp:.3e}"].diff()/self.data["times"].diff())
            self.data.at[0, "Imembrane_"+str(clamp)] = 0.0
            self.data["Imembrane_"+str(clamp)] = self.data["Imembrane_"+str(clamp)] - self.data["Imembrane_"+str(clamp)][0]
        return self.data
    
    def filter_membrane_currents(self):
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        for inj in list(self.parameters["Iinj"]):
            self.data["filtered_Imembrane_"+str(inj)] = self.filters["membrane_currents"].propagate(self.data["Imembrane_"+str(inj)], sampling_rate)
        return self.data
    
    def filter_activation_currents(self):
        sampling_rate = 1/(self.data["times"][1] - self.data["times"][0])
        for inj in list(self.parameters["Iinj"]):
            activation_current = self.filters["activation_currents"].propagate(self.data["Iactivation_"+str(inj)], sampling_rate)
            self.data["filtered_Iactivation_"+str(inj)] = activation_current 
        return self.data
    
    def compute_passive_conductances(self):
        ntimesteps = self.data.shape[0]
        A = np.zeros((ntimesteps, 2, 2))
        B = np.zeros((ntimesteps, 2, 1))
        membrane_potential = self.data[[f"{x:.3e}" for x in list(self.parameters["Iinj"])]].to_numpy()
        membrane_current = self.data[["filtered_Imembrane_"+str(x) for x in list(self.parameters["Iinj"])]].to_numpy()
        activation_current = self.data[["filtered_Iactivation_"+str(x) for x in list(self.parameters["Iinj"])]].to_numpy()
        leakage_current = self.data[["Ileakage_"+str(x) for x in list(self.parameters["Iinj"])]].to_numpy()
        excitatory_reversal_potential = self.parameters["Ee"].to_numpy()
        inhibitory_reversal_potentail = self.parameters["Ei"].to_numpy()
        injected_current = self.parameters["Iinj"].to_numpy()
        A[:, 0, 0] = np.sum(np.square(membrane_potential - excitatory_reversal_potential), axis=1)
        A[:, 0, 1] = np.sum((membrane_potential - excitatory_reversal_potential) * (membrane_potential- inhibitory_reversal_potentail), axis=1)
        A[:, 1, 0] = A[:, 0, 1]
        A[:, 1, 1] = np.sum(np.square(membrane_potential - inhibitory_reversal_potentail), axis=1)
        B[:, 0, 0] = -1.0*np.sum((membrane_current - activation_current - injected_current + leakage_current)*(membrane_potential - excitatory_reversal_potential), axis=1)
        B[:, 1, 0] = -1.0*np.sum((membrane_current - activation_current - injected_current + leakage_current)*(membrane_potential - inhibitory_reversal_potentail), axis=1)
        conductances = np.linalg.pinv(A) @ B
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
    
    def get_clamp_near_0(self) -> Tuple[int, float]:
        index_of_minimum_injected_current: int = np.argmin(np.abs(self.parameters["Iinj"]))
        minimum_injected_current: float = self.parameters["Iinj"][index_of_minimum_injected_current]
        return index_of_minimum_injected_current, minimum_injected_current
    
    def compute_stats(self):
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

    def estimate_conductances(self):
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

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from pyhelpers.store import save_fig
from pathlib import Path
from scipy.optimize import minimize, Bounds, lsq_linear
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
from scipy import sparse

from libs.readers import XLReader, AnalyzerCfg

from typing import Dict, List, Tuple

import numpy as np

def weighted_quantile(x, q, w=None, axis=-1):
    """
    Weighted quantiles of `x` at quantiles `q` in [0,1], along `axis`.

    Parameters
    ----------
    x : array_like
        Data.
    q : float or array_like
        Quantile(s) in [0,1].
    w : array_like or None
        Nonnegative weights, same shape as x, or broadcastable to x.
        If None, equivalent to np.quantile(x, q, axis=axis) but using sorting.
    axis : int
        Axis to compute along.

    Returns
    -------
    out : ndarray
        Shape is x.shape with `axis` removed, plus an extra quantile dimension (len(q)).
        If q is scalar, that quantile dimension is omitted.
    """
    x = np.asarray(x)
    q = np.asarray(q, dtype=float)
    if np.any((q < 0) | (q > 1)):
        raise ValueError("q must be in [0, 1].")

    if w is None:
        w = np.ones_like(x, dtype=float)
    else:
        w = np.asarray(w, dtype=float)
        w = np.broadcast_to(w, x.shape)

    if np.any(w < 0):
        raise ValueError("weights must be nonnegative.")

    # Move target axis to last for easier vectorization
    x = np.moveaxis(x, axis, -1)
    w = np.moveaxis(w, axis, -1)

    *batch, n = x.shape
    m = int(np.prod(batch)) if batch else 1
    x2 = x.reshape(m, n)
    w2 = w.reshape(m, n)

    # Sort each row
    idx = np.argsort(x2, axis=1)
    xs = np.take_along_axis(x2, idx, axis=1)
    ws = np.take_along_axis(w2, idx, axis=1)

    # CDF of weights
    cw = np.cumsum(ws, axis=1)
    total = cw[:, -1]
    if np.any(total <= 0):
        raise ValueError("each slice must have positive total weight.")

    cdf = cw / total[:, None]  # in (0,1]

    # Find first index where CDF >= q (vectorized)
    qv = q.ravel()
    mask = cdf[:, None, :] >= qv[None, :, None]          # (m, k, n)
    any_true = mask.any(axis=2)
    hi = mask.argmax(axis=2)                              # (m, k)
    hi = np.where(any_true, hi, n - 1)                    # safety (shouldn't trigger if q<=1)

    # Linear interpolation in (cdf, x) space
    lo = np.clip(hi - 1, 0, n - 1)

    x_hi = np.take_along_axis(xs, hi, axis=1)
    x_lo = np.take_along_axis(xs, lo, axis=1)

    c_hi = np.take_along_axis(cdf, hi, axis=1)
    c_lo = np.where(hi > 0, np.take_along_axis(cdf, lo, axis=1), 0.0)

    denom = (c_hi - c_lo)
    # If denom==0 (flat CDF step), fall back to x_hi
    t = np.where(denom > 0, (qv[None, :] - c_lo) / denom, 0.0)
    t = np.clip(t, 0.0, 1.0)

    out = x_lo + t * (x_hi - x_lo)

    # Reshape back: batch dims + (k,)
    k = qv.size
    out = out.reshape((*batch, k))

    # If q was scalar, drop the last dim
    if q.ndim == 0:
        out = out[..., 0]

    return out


class WholeCellRecording:
    def __init__(self, parameters: pd.DataFrame, bin_s: float, stimuli: Dict[str, pd.DataFrame]) -> None:
        assert len(stimuli) > 0
        
        self.Cm: float = parameters["Cm"][0]    # units: Farads
        self.Et: float = parameters["Et"][0]

        example_t: pd.Series = list(stimuli.values())[0]["times"] # units: Seconds, shape: (Nsamples,)
        self.dt: float = example_t[1] - example_t[0] # time assumed sampled at constant interval; units: Seconds

        self.bin_nsamples: int = round(bin_s / self.dt)
        self.bin_s: float = self.bin_nsamples * self.dt # units: Seconds

        self.Ee: float = np.nan # units: Volts
        self.Ei: float = np.nan # units: Volts
        
        self.LI_dvdt_vs_Vm: float = np.nan

        self.stimuli: Dict[str, WholeCellStimulus] = {name: WholeCellStimulus(self, data) for name, data in stimuli.items()}
        
        self.Er, self.Rin = self._estimate_Er_Rin()
        
        self.Eact: float = np.nan
        self.Eact = self._estimate_Eact()
        
        print(f"Parameter Estimates:\n\tEr: {self.Er*1e3:.1f} mV\n\tRin: {self.Rin*1e-9:.1f} Gohm\n\tEact: {self.Eact*1e3:.1f} mV")


    @property
    def gl(self) -> float:
        return 1 / self.Rin # units: Siemens

    @property
    def gact(self) -> float:
        return self.gl * (self.max_Eact - self.Er) / (self.max_Eact - self.Eact)

    @property
    def max_Vm(self) -> float:
        max_Vms: List[float] = []
        for stimulus in self.stimuli.values():
            max_Vms.append(stimulus.Vm.max())
        return np.max(max_Vms)
    
    @property
    def max_Eact(self) -> float:
        return max(self.Et, self.max_Vm)

    def _estimate_Er_Rin(self, bins: int=200, smooth_bins: float=4) -> Tuple[float, float]:
        # Step 1: Get noisy estimate of Vss: argmax of Vm density function.
        Vm_by_Iinj: Dict[float, List[np.ndarray]] = {}
        for stimulus in self.stimuli.values():
            for Iinj, Vm in zip(np.squeeze(stimulus.Iinj), stimulus.Vm):
                if Iinj in Vm_by_Iinj:
                    Vm_by_Iinj[Iinj].append(Vm)
                else:
                    Vm_by_Iinj[Iinj] = [Vm]
        
        fig, ax = plt.subplots(2, 1)

        Iinjs_list: List[float] = []
        ys_list: List[float] = []
        for Iinj, Vms in Vm_by_Iinj.items():
            counts, edges = np.histogram(np.concatenate(Vms), bins=bins)
            counts_s: np.ndarray = gaussian_filter1d(counts.astype(float), smooth_bins)
            centers: np.ndarray = 0.5 * (edges[:-1] + edges[1:])
            ax[0].plot(centers, counts_s)
            Vss_hat: float = centers[np.argmax(counts_s)]

            Iinjs_list.append(Iinj)
            ys_list.append(Vss_hat)

        # Step 2: Linear Regression Solution
        Iinjs: np.ndarray = np.array(Iinjs_list)
        Iinjs_mean: float = Iinjs.mean()
        Iinjs_centered: np.ndarray = Iinjs - Iinjs_mean

        ys: np.ndarray = np.array(ys_list)
        ys_mean: float = ys.mean() 
        ys_centered: np.ndarray = ys - ys_mean 
        
        Rin_hat: float = np.sum(Iinjs_centered * ys_centered) / np.sum(np.square(Iinjs_centered))
        Er_hat: float = ys_mean - Rin_hat * Iinjs_mean

        ax[1].scatter(Iinjs, ys)
        ax[1].plot(Iinjs, Er_hat + Rin_hat * Iinjs)

        return Er_hat, Rin_hat

    def estimate_Ee_Ei(self) -> Tuple[float, float]:
        Eeff_pool: List[float] = []
        for stimulus in self.stimuli.values():
            Eeff_pool.extend(stimulus.binned_timeseries["Eeff"])

        Ei_hat: float = float(np.nanquantile(Eeff_pool, 0.05)) # units: Volts
        Ee_hat: float = float(np.nanquantile(Eeff_pool, 0.95)) # units: Volts

        return Ee_hat, Ei_hat

    def estimate_Ee_Ei_from_Eeff(self, series_key: str, q_low: float = 0.05, q_high: float = 0.95) -> Tuple[float, float]:
        """
        Estimate Ee/Ei from pooled Eeff values stored in each stimulus' binned_timeseries[series_key].
        NaNs are ignored.
        Returns (Ee_hat, Ei_hat).
        """
        pool: List[float] = []
        for stimulus in self.stimuli.values():
            if series_key not in stimulus.binned_timeseries:
                continue
            vals = stimulus.binned_timeseries[series_key]
            # vals may be numpy array, list, or pandas Series
            pool.extend(list(np.asarray(vals, dtype=np.float64).ravel()))
        if len(pool) == 0:
            return np.nan, np.nan
        Ei_hat = float(np.nanquantile(pool, q_low))
        Ee_hat = float(np.nanquantile(pool, q_high))
        return Ee_hat, Ei_hat
        

    def run_analysis(self, verbose: bool=True):
        # --- Step 1: estimate reversal potentials from Eeff/gsyn stage ---
        if verbose: print("Estimating reversal potentials... ")

        for stimulus in self.stimuli.values():
            stimulus.calculate_target_Qsyn()
            stimulus.calculate_target_Qsyn_nonlinear()

            # Original (stepwise) Eeff/gsyn estimator (kept for backward compatibility)
            stimulus.estimate_Eeff()

            # New (piecewise-linear) Eeff/gsyn estimator (for side-by-side comparison)
            # stimulus.estimate_Eeff_gsyn_piecewise_linear(nonlinear=False, constrained=True)

        # Original reversal estimate (used by default downstream)
        self.Ee, self.Ei = self.estimate_Ee_Ei()

        # Alternate reversal estimate from piecewise-linear Eeff nodes (stored for comparison)
        self.Ee_pl, self.Ei_pl = self.estimate_Ee_Ei_from_Eeff("Eeff nodes piecewise linear constrained")

        if verbose:
            print(f"Estimated Reversals (stepwise Eeff): Ee = {self.Ee*1e3:.1f} mV, Ei = {self.Ei*1e3:.1f} mV")
            if np.isfinite(self.Ee_pl) and np.isfinite(self.Ei_pl):
                print(f"Estimated Reversals (piecewise-linear Eeff): Ee = {self.Ee_pl*1e3:.1f} mV, Ei = {self.Ei_pl*1e3:.1f} mV")

        # --- Step 2: estimate ge/gi and forward-predict Vm ---
        if verbose: print("Calculating synaptic conductances... ")

        for stimulus in self.stimuli.values():
            # Original stepwise ge/gi + analytic forward model (unchanged outputs)
            stimulus.estimate_ge_gi(nonlinear=True, constrained=True, fit_ge=True, fit_gi=True)
            stimulus.calculate_predicted_Vm()

            # New piecewise-linear ge/gi + Crank–Nicolson forward model (side-by-side)
            stimulus.estimate_ge_gi_piecewise_linear(nonlinear=True, constrained=True)
            stimulus.calculate_predicted_Vm_piecewise_linear(nonlinear=True, constrained=True)
    """ Eact Optimization """
    def predicted_Vm_SSE(self, x) -> float:
        self.Eact = x
        self.run_analysis(verbose=False)
        total = 0.0
        for stimulus in self.stimuli.values():
            sse = np.sum(np.square(stimulus.Vm - np.stack(stimulus.timeseries["predicted Vm"].to_numpy()).astype(np.float64).T)) #type: ignore
            total += sse
        return total
    
    def gegi_correlation(self, x) -> float:
        self.Eact = x
        self.run_analysis(verbose=False)
        total = 0.0
        for stimulus in self.stimuli.values():
            sse = np.sum(np.square(stimulus.Vm - np.stack(stimulus.timeseries["predicted Vm"].to_numpy()).astype(np.float64).T)) #type: ignore
            total += sse
        return total
    
    def _estimate_Eact(self) -> float:
        return self.max_Eact
        print(self.max_Vm, self.max_Eact, self.Er)
        eps: float = np.finfo(float).eps
        if self.max_Eact <= self.Er:
            return self.Er + eps
        else:
            Eact_initial_guess: float = (self.Er + self.max_Eact) / 2
            result = minimize(
                self.predicted_Vm_SSE, 
                Eact_initial_guess, 
                method='Powell', 
                bounds=Bounds(self.Er + eps, self.max_Eact - eps),
            )
            return result.x[0]


class WholeCellStimulus:
    def __init__(self, recording: WholeCellRecording, data: pd.DataFrame) -> None:
        self.recording: WholeCellRecording = recording

        self.times: np.ndarray = data["times"].to_numpy(dtype=np.float64)

        Iinj_colnames: List[str] = list(data.columns)
        
        for aux_key in ["times", "stimulus", "representative"]:
            if aux_key in Iinj_colnames:
                Iinj_colnames.remove(aux_key)

        # (* 1e-3 is to scale Vm from millivolts to volts)
        self.Vm: np.ndarray = data[Iinj_colnames].to_numpy(dtype=np.float64).T * 1e-3       # units: Volts, shape: [Nclamps, Nsamples]
        self.Iinj: np.ndarray = np.array(Iinj_colnames, dtype=np.float64)[:, np.newaxis]    # units: Amperes, shape: [Nclamps, 1]
                
        self.Nclamps: int = self.Iinj.size
        self.Nsamples: int = self.times.size

        self.timeseries: pd.DataFrame = pd.DataFrame({"times": self.times})

    def calculate_target_Qsyn_nonlinear(self) -> None:

        integral_Vm: np.ndarray = np.stack(self.binned_timeseries["integral Vm"].to_numpy()).astype(np.float64).T #type: ignore
        target_Qsyn: np.ndarray = np.stack(self.binned_timeseries["target Qsyn"].to_numpy()).astype(np.float64).T #type: ignore
        bin_duration: np.ndarray = self.binned_timeseries["bin duration"].to_numpy(np.float64)  # units: Seconds, shape: (Nbins,)

        integral_Iact: np.ndarray = np.maximum(0.0, self.recording.gact * (integral_Vm - self.recording.Eact * bin_duration[np.newaxis, :]))
        target_Qsyn_nonlinear: np.ndarray = target_Qsyn - integral_Iact

        self.binned_timeseries["integral Iact"] = integral_Iact.T.tolist()
        self.binned_timeseries["target Qsyn nonlinear"] = target_Qsyn_nonlinear.T.tolist()

    def calculate_target_Qsyn(self) -> None:

        l_bin_idxs: np.ndarray = self.binned_timeseries["left index"].to_numpy(np.int32)           # shape: (Nbins,)
        r_bin_idxs: np.ndarray = self.binned_timeseries["right index"].to_numpy(np.int32)           # shape: (Nbins,)
        bin_duration: np.ndarray = self.binned_timeseries["bin duration"].to_numpy(np.float64)  # units: Seconds, shape: (Nbins,)

        Vm_prefix_integral: np.ndarray = self._cumtrapz_prefix_integral(self.Vm)  # units: Webers, shape: (Nclamps, Nsamples + 1)
        integral_Vm: np.ndarray = Vm_prefix_integral[:, r_bin_idxs - 1] - Vm_prefix_integral[:, l_bin_idxs] # units : Webers, shape: (Nclamps, Nsamples)

        Delta_Vm: np.ndarray = self.Vm[:, r_bin_idxs - 1] - self.Vm[:, l_bin_idxs] # units: Volts, shape: (Nclamps, Nsamples - 1)

        integral_Il: np.ndarray = self.recording.gl * (self.recording.Er * bin_duration[np.newaxis, :] - integral_Vm) # units: Coulombs, shape: (Nclamps, Nbins)
        integral_Iinj: np.ndarray = bin_duration[np.newaxis, :] * self.Iinj    # units: Coulombs, shape: (Nclamps, Nbins)
        target_Qsyn: np.ndarray = self.recording.Cm * Delta_Vm - (integral_Il + integral_Iinj)      # units: Coulombs, shape: (Nclamps, Nbins)

        self.binned_timeseries["integral Vm"] = integral_Vm.T.tolist()
        self.binned_timeseries["integral Il"] = integral_Il.T.tolist()
        self.binned_timeseries["integral Iinj"] = integral_Iinj.T.tolist()
        self.binned_timeseries["target Qsyn"] = target_Qsyn.T.tolist()

    def estimate_Eeff(self) -> None:
        # Pull data (faster than .tolist() if these are arrays-in-cells, but keep if needed)
        integral_Vm = np.stack(self.binned_timeseries["integral Vm"].to_numpy()).astype(np.float64).T #type: ignore , (Nclamps, Nbins)
        target_Qsyn = np.stack(self.binned_timeseries["target Qsyn"].to_numpy()).astype(np.float64).T #type: ignore , (Nclamps, Nbins)

        # Means per bin
        integral_Vm_mean = integral_Vm.mean(axis=0) # (Nbins,)
        target_Qsyn_mean = target_Qsyn.mean(axis=0)

        # Centered
        centered_integral_Vm = integral_Vm - integral_Vm_mean[None, :]
        centered_target_Qsyn = target_Qsyn - target_Qsyn_mean[None, :]

        # Regression slope b and intercept a for each bin
        denom = np.sum(centered_integral_Vm * centered_integral_Vm, axis=0) # var * (Nclamps-1) up to scale
        numer = np.sum(centered_integral_Vm * centered_target_Qsyn, axis=0)

        # Handle degenerate bins where Phi has no variation across clamps
        eps = np.finfo(np.float64).tiny
        b = np.where(np.abs(denom) > eps, numer / denom, np.nan)   # slope (Nbins,)
        a = target_Qsyn_mean - b * integral_Vm_mean                                        # intercept

        # Your derived params
        gsyn = -b                                                  # Siemens
        # Avoid divide-by-zero when gsyn ~ 0
        gsyn_safe = np.where(np.abs(gsyn) > 0.1e-9, gsyn, np.nan)

        Eeff = a / (gsyn_safe * float(self.recording.bin_s))       # Volts

        # Diagnostics per bin
        # Your SSE formula: sum (Q - gsyn*(Eeff - Phi))^2
        # We can compute predicted Q directly from a + b*Phi (same fit)
        Q_hat = a[None, :] + b[None, :] * integral_Vm
        resid = target_Qsyn - Q_hat
        sse = np.sum(resid * resid, axis=0)

        # Proper per-bin R^2: 1 - SSE / SST, SST = sum (Q - mean(Q))^2 within the bin
        sst = np.sum((target_Qsyn - target_Qsyn_mean[None, :]) ** 2, axis=0)
        r2 = np.where(sst > 0, 1.0 - (sse / sst), np.nan)

        # Store
        self.binned_timeseries["Eeff"] = Eeff
        self.binned_timeseries["gsyn"] = gsyn_safe
        self.binned_timeseries["SSE least-squares Qsyn"] = sse
        self.binned_timeseries["r2 least-squares Qsyn"] = r2

    def estimate_ge_gi(self) -> None:
        # Pull arrays
        Vm = self.Vm
        Im = self.Im
        Ia = self.Ia
        Il = self.Il
        Iinj = self.Iinj

        Ee = self.Ee 
        Ei = self.Ei

        y = -(Im - Ia - Iinj + Ileak)

        X = np.stack([Vm - Ee, Vm - Ei], axis=-1)

        XtX = np.einsum("nki,nkj->nij", X, X)
        Xty = np.einsum("nki,nk->ni", X, y)

        try:
            conductances = np.linalg.solve(XtX, Xty)   # shape (Nsamples, 2)
        except np.linalg.LinAlgError:
            conductances = (np.linalg.pinv(XtX) @ Xty[..., None])[..., 0]

        self.data["excitation"] = conductances[:, 0]
        self.data["inhibition"] = conductances[:, 1]

        self.data["positive_excitation"] = self.data["excitation"].clip(lower=0.0)
        self.data["positive_inhibition"] = self.data["inhibition"].clip(lower=0.0)

        self.data["resultant_excitation"] = (
            self.data["positive_excitation"] - self.data["positive_inhibition"]
        ).clip(lower=0.0)

        self.data["resultant_inhibition"] = (
            self.data["positive_inhibition"] - self.data["positive_excitation"]
        ).clip(lower=0.0)

        return self.data
    
    def calculate_predicted_Vm(self) -> None:
        Vm: np.ndarray = self.Vm        # units: Volts, shape: (Nclamps, Nsamples)
        Iinj: np.ndarray = self.Iinj    # units: Amperes, shape: (Nclamps, 1)

        Ee: float = float(self.recording.Ee)    # units: Volts
        Ei: float = float(self.recording.Ei)    # units: Volts
        Er: float = float(self.recording.Er)    # units: Volts
        gl: float = float(self.recording.gl)    # units: Siemens
        Cm: float = float(self.recording.Cm)    # units: Farads
        dt: float = float(self.recording.dt)    # units: Seconds

        l_idx: np.ndarray = self.binned_timeseries["left index"].to_numpy(dtype=np.int64)   # shape: (Nbins,)
        r_idx: np.ndarray = self.binned_timeseries["right index"].to_numpy(dtype=np.int64)  # shape: (Nbins,)

        ge_bins: np.ndarray = self.binned_timeseries["ge nonlinear constrained"].to_numpy(dtype=np.float64) # units: Siemens, shape: (Nbins,)
        gi_bins: np.ndarray = self.binned_timeseries["gi nonlinear constrained"].to_numpy(dtype=np.float64) # units: Siemens, shape: (Nbins,)

        Vpred: np.ndarray = np.empty_like(Vm, dtype=np.float64) # units: Volts, shape: (Nclamps, Nsamples)
        Vpred[:, 0:1] = Er + Iinj / gl

        vsss = []

        for r, l, ge, gi in zip(r_idx, l_idx, ge_bins, gi_bins):
            v0: np.ndarray = Vpred[:, l].reshape(-1, 1)
            G: np.ndarray = ge + gi + gl
            vss: np.ndarray = (ge * Ee + gi * Ei + gl * Er + Iinj) / G
            t: np.ndarray = np.arange(r - l + 1) * dt
            Vpred[:, l:r + 1] = (v0 - vss) * np.exp(-G * t / Cm) + vss

            vsss.append(vss)

        self.timeseries[f"predicted Vm"] = Vpred.T.tolist()


class Analyzer:
    def __init__(self, cfg: AnalyzerCfg):
        self.cfg: AnalyzerCfg = cfg

    def plot_timeseries(
        self,
        recording: WholeCellRecording,
        filename: Path,
        filetype: str = "png",
        cond_warn: float = 1e6,
        cond_bad: float = 1e10,
        resnorm_factor_warn: float = 5.0,
        display: bool = True
    ):
        """
        Rows:
        0) Vm (solid) + Vpred (dotted), per Iinj color
        1) Im per Iinj color
        2) Il per Iinj color
        3) ge (red) and gi (blue)
        4) diagnostics: residual norm (warn line based on robust baseline)
        5) diagnostics: conditioning (log scale with warn/bad lines)
        """

        ncols = len(recording.stimuli)
        nrows = 7
        fig, axs = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            sharex="all",
            sharey="row",
            figsize=(15, 10),
            constrained_layout=True,
        )
        fig.suptitle(f"{filename.name}")

        if ncols == 1:
            axs = np.expand_dims(axs, axis=1)

        for idx, paradigm in enumerate(recording.stimuli):
            stimulus = recording.stimuli[paradigm]

            times: np.ndarray = stimulus.times
            bin_times: np.ndarray = (stimulus.binned_timeseries["left time"].to_numpy(np.float64) + stimulus.binned_timeseries["right time"].to_numpy(np.float64)) / 2
            node_times: np.ndarray = stimulus.node_timeseries["node time"].to_numpy(np.float64)

            bin_durations: np.ndarray = stimulus.binned_timeseries["bin duration"].to_numpy(np.float64)

            Vpred: np.ndarray = np.stack(stimulus.timeseries["predicted Vm"].to_numpy()).astype(np.float64).T #type: ignore
            Vpred_pw: np.ndarray = np.stack(stimulus.timeseries["predicted Vm piecewise linear nonlinear constrained"].to_numpy()).astype(np.float64).T #type: ignore
            integral_Iact: np.ndarray = np.stack(stimulus.binned_timeseries["integral Iact"].to_numpy()).astype(np.float64).T #type: ignore

            Eeff: np.ndarray = stimulus.binned_timeseries["Eeff"].to_numpy(np.float64) 
            # Eeff_pw: np.ndarray = stimulus.node_timeseries["Eeff nodes piecewise linear constrained"].to_numpy(np.float64)
            
            ge = stimulus.binned_timeseries["ge nonlinear constrained"].to_numpy(np.float64)
            gi = stimulus.binned_timeseries["gi nonlinear constrained"].to_numpy(np.float64)
            ge_pw = stimulus.node_timeseries["ge nodes piecewise linear nonlinear constrained"].to_numpy(np.float64)
            gi_pw = stimulus.node_timeseries["gi nodes piecewise linear nonlinear constrained"].to_numpy(np.float64)
            
            # diagnostics
            resnorm = stimulus.binned_timeseries["ge/gi residual norm nonlinear constrained"].to_numpy()
            cond = stimulus.binned_timeseries["ge/gi matrix conditioning nonlinear constrained"].to_numpy()

            axs[0, idx].set_title(paradigm)

            # ---- Row 0: Vm + Vpred (same colors) ----
            curves = axs[0, idx].plot(times, stimulus.Vm.T)  # solid Vm traces
            colors = [line.get_color() for line in curves]

            Ess = recording.Er + recording.Rin * stimulus.Iinj
            for j in range(stimulus.Nclamps):
                axs[0, idx].plot(times, Vpred[j, :], linestyle=":", color=colors[j])  # dotted Vpred
                axs[0, idx].plot(times, Vpred_pw[j, :], linestyle="--", color=colors[j])
                axs[0, idx].plot([times[0], times[-1]], [Ess[j]] * 2, color="grey", ls="--")
            axs[0, idx].plot([times[0], times[-1]], [recording.Eact] * 2, color="red", ls="--")

            axs[0, idx].set_ylabel("Vm (V)")
            axs[0, idx].grid(True)
            # if idx == 0:
            #     axs[0, idx].legend(loc="upper right")

            # ---- Row 1: Eeff ----
            axs[1, idx].grid(True)
            axs[1, idx].plot(bin_times, Eeff, c="black")
            # axs[1, idx].plot(node_times, Eeff_pw, c="black", ls="--")
            axs[1, idx].axhline(recording.Er, linestyle="--", color="k", linewidth=1, label="Er" if idx == 0 else None)
            axs[1, idx].axhline(recording.Ee, linestyle="--", color="r", linewidth=1, label="Ee" if idx == 0 else None)
            axs[1, idx].axhline(recording.Ei, linestyle="--", color="b", linewidth=1, label="Ei" if idx == 0 else None)

            axs[2, idx].grid(True)
            axs[2, idx].plot(bin_times, bin_times * 0, "k--")
            axs[2, idx].plot(bin_times, integral_Iact.T / bin_durations[:, np.newaxis])

            # ---- Row 3: conductances ----
            axs[3, idx].plot(bin_times, ge, c="r", label="ge")
            axs[3, idx].plot(bin_times, gi, c="b", label="gi")
            axs[3, idx].plot(node_times, ge_pw, c="r", ls=":")
            axs[3, idx].plot(node_times, gi_pw, c="b", ls=":")

            axs[3, idx].plot(bin_times, bin_times * 0, "--k", linewidth=1)
            axs[3, idx].set_ylabel("G (S)")
            axs[3, idx].grid(True)
            # if idx == 0:
            #     axs[3, idx].legend(loc="upper right")

            # ---- Row 4: residual norm + warning threshold ----
            ax_r = axs[4, idx]
            ax_r.grid(True)
            ax_r.set_ylabel("resnorm")

            ax_r.plot(bin_times, resnorm, linewidth=1)

            # robust baseline and warning line
            finite_r = resnorm[np.isfinite(resnorm)]
            if finite_r.size > 0:
                baseline = np.median(finite_r)
                warn_line = resnorm_factor_warn * baseline
                ax_r.axhline(warn_line, linestyle="--", linewidth=1)
                ax_r.text(
                    0.01, 0.95,
                    f"warn > {resnorm_factor_warn:g}×median",
                    transform=ax_r.transAxes,
                    va="top",
                )

            # ---- Row 5: conditioning + warning/bad thresholds ----
            ax_c = axs[5, idx]
            ax_c.grid(True)
            ax_c.set_ylabel("cond(XᵀX)")
            ax_c.set_yscale("log")

            cond_plot = np.where(np.isfinite(cond) & (cond > 0), cond, np.nan)
            ax_c.plot(bin_times, cond_plot, linewidth=1)

            # thresholds
            ax_c.axhline(cond_warn, linestyle="--", linewidth=1)
            ax_c.axhline(cond_bad, linestyle=":", linewidth=1)
            ax_c.text(
                0.01, 0.95,
                f"warn>{cond_warn:.0e}  bad>{cond_bad:.0e}",
                transform=ax_c.transAxes,
                va="top",
            )

            axs[5, idx].set_xlabel("time (s)")

        out = f"{str(filename)}_dev_level1."
        if filetype == "png":
            plt.savefig(f"{out}{filetype}")
        elif filetype == "emf":
            save_fig(f"{out}svg", dpi=300, conv_svg_to_emf=True, verbose=True)
        else:
            raise ValueError(f"Unsupported filetype: {filetype}")

        if display:
            plt.show()

    def run(self, display: bool) -> None:
        n_files: int = len(self.cfg.paths_to_spreadsheets)

        for i in range(n_files):
            path_to_spreadsheet: Path = self.cfg.paths_to_spreadsheets[i]
            print(f"Collecting data from {path_to_spreadsheet}")
            rdr: XLReader = XLReader(path_to_spreadsheet)
            
            stimuli: Dict[str, pd.DataFrame] = {paradigm:rdr.get_paradigm_data(paradigm) for paradigm in rdr.get_paradigms()}
            recording: WholeCellRecording = WholeCellRecording(rdr.get_paradigm_parameters(rdr.get_paradigms()[0]), 5e-3, stimuli)
            del rdr  # free excel file handle
            
            recording.run_analysis()

            # ---- plot ----
            if display:
                self.plot_timeseries(
                    recording,
                    self.cfg.image_save_dir / path_to_spreadsheet.stem,
                    filetype=self.cfg.image_save_type,
                    display=display,
                )

# class Analyzer:
#     def __init__(self, filepaths: List[str], output_path: Path):
#         self.filepaths: List[str] = filepaths
#         self.output_path: str = output_path

#     def get_paradigm_to_optimize(self, recordings):
#         analysis_logger.info("Finding the best paradigm to optimize")
#         maximum_depolarizations = np.zeros(len(recordings))
#         for idx, paradigm in enumerate(recordings):
#             index_of_minimum_injected_current, minimum_injected_current = recordings[paradigm].get_clamp_near_0()
#             membrane_potential_at_minimum_injected_current = recordings[paradigm].data[minimum_injected_current]
#             depols = membrane_potential_at_minimum_injected_current - recordings[paradigm].parameters["Ess"][index_of_minimum_injected_current]
#             maximum_depolarizations[idx] = np.max(depols)
#         max_index = np.argmax(maximum_depolarizations)
#         return list(recordings.keys())[max_index]

#     def estimation_without_optim_activation_potential(self, recordings):
#         for idx, paradigm in enumerate(recordings):
#             recordings[paradigm].estimate_conductances()
#             recordings[paradigm].stats.insert(0, "paradigm", paradigm)
#             if idx == 0:
#                 overall_stats = recordings[paradigm].stats.copy()
#             else:
#                 overall_stats = pd.concat([overall_stats, recordings[paradigm].stats], axis=0)
#             print(f"{paradigm} done")
#         overall_stats = overall_stats.sort_values(by="paradigm", ascending=True)
#         return recordings, overall_stats
    
#     def write_to_excel(self, filepath: Path, recordings, stats):
#         if not self.sysops.check_directory(filepath):
#             pd.DataFrame().to_excel(filepath)
#         with pd.ExcelWriter(filepath, mode='a', engine='openpyxl', if_sheet_exists='new') as writer:
#             for paradigm in recordings:
#                 recordings[paradigm].data.to_excel(writer, sheet_name=paradigm, index=False)
#                 recordings[paradigm].parameters.to_excel(writer, sheet_name="parameters_"+paradigm, index = False)
#             stats.to_excel(writer, sheet_name="stats", index=False)
#         pass

#     def write_analysis_to_excel(self, filepath: Path, paradigm: str, paradigm_data: pd.DataFrame):
#         if self.sysops.check_directory(filepath):
#             with pd.ExcelWriter(filepath, mode='a', engine='openpyxl', if_sheet_exists='replace') as writer:
#                 paradigm_data.to_excel(writer, sheet_name=paradigm, index=False)
#         else:
#             # File does not exist, create a new file
#             paradigm_data.to_excel(filepath, sheet_name=paradigm, index=False)
#         pass

#     def plot_dev(self, recordings, filename: Path, filetype: str="png", current_clamps: Optional[List[float]]=None):
#         analysis_logger.info(f"Verbose plotting of conductance estimations for {filename}")
#         fig, axs = plt.subplots(nrows = 7, ncols = len(recordings), sharex="all", sharey="row", figsize=(15, 10), constrained_layout=True)
#         for idx, paradigm in enumerate(recordings):
#             paradigm_iinj: List[float] = list(recordings[paradigm].parameters["Iinj"])
#             if current_clamps is not None and set(current_clamps) <= set(paradigm_iinj):
#                 paradigm_iinj = list(set(paradigm_iinj).intersection(current_clamps))
#                 assert len(paradigm_iinj) > 0, "Cannot plot. Specified current clamps have no intersection with current clamps listed in parameters."

#             if "representative" in recordings[paradigm].data:
#                 rep = recordings[paradigm].data["representative"].to_numpy()
#             membrane_potential = recordings[paradigm].data[[f"{x:.3e}" for x in paradigm_iinj]].to_numpy()
#             membrane_current = recordings[paradigm].data[["filtered_Imembrane_"+str(x) for x in paradigm_iinj]].to_numpy()
#             leakage_current = recordings[paradigm].data[["Ileakage_"+str(x) for x in paradigm_iinj]].to_numpy()
#             activation_current = recordings[paradigm].data[["filtered_Iactivation_"+str(x) for x in paradigm_iinj]].to_numpy()
#             conductances = recordings[paradigm].data[["excitation", "inhibition"]].to_numpy()
#             times = recordings[paradigm].data["times"].to_numpy()
#             resting_potential = times*0 + recordings[paradigm].parameters["Er"][0]
#             threshold_potential = times*0 + recordings[paradigm].parameters["Et"][0]
#             activation_potential = membrane_potential*0 + recordings[paradigm].parameters["Eact"].to_numpy()
#             if "stimulus" in recordings[paradigm].data:
#                 stim = recordings[paradigm].data["stimulus"].to_numpy()
#             axs[0, idx].set_title(paradigm)
#             if "representative" in recordings[paradigm].data:
#                 axs[0, idx].plot(times, rep)
#                 axs[0, idx].plot(times, resting_potential, '--k')
#             axs[0, idx].set_ylabel("Rep. Vm (V)")
#             axs[0, idx].set_title(paradigm)
#             axs[0, idx].grid(True)
#             curves = axs[1, idx].plot(times, membrane_potential)
#             colors = [x.get_color() for x in curves]
#             axs[1, idx].plot(times, resting_potential, '--k')
#             axs[1, idx].plot(times, threshold_potential, linestyle='--', color=(0.5, 0.5, 0.5))
#             for i in range(activation_potential.shape[-1]):
#                 axs[1, idx].plot(times, activation_potential[:, i], linestyle='--', color=colors[i])
#             axs[1, idx].set_ylabel("Vm (V)")
#             axs[1, idx].grid(True)
#             axs[2, idx].plot(times, membrane_current)
#             axs[2, idx].set_ylabel("Im (A)")
#             axs[2, idx].grid(True)
#             axs[3, idx].plot(times, leakage_current)
#             axs[3, idx].set_ylabel("Ileak (A)")
#             axs[3, idx].grid(True)
#             axs[4, idx].plot(times, activation_current)
#             axs[4, idx].set_ylabel("Iact (A)")
#             axs[4, idx].grid(True)
#             # axs[4, idx].set_title(", ".join([f"{val:.3f}" for val in recordings[paradigm].parameters["Eact"]]))
#             axs[5, idx].plot(times, conductances[:, 0], c='r')
#             axs[5, idx].plot(times, conductances[:, 1], c='b')
#             axs[5, idx].plot(times, times*0, '--k')
#             axs[5, idx].set_ylabel("G (S)")
#             axs[5, idx].grid(True)
#             if "stimulus" in recordings[paradigm].data:
#                 axs[6, idx].plot(times, stim)
#             axs[6, idx].set_ylabel("Stimulus")
#             axs[6, idx].set_xlabel("times(sec)")
#             axs[6, idx].grid(True)
#         analysis_logger.info(f"Saving verbose plotting of conductance estimations for {filename}")
#         filename: str = f"{str(filename)}_dev_traces."
#         if filetype == "png":
#             plt.savefig(f"{filename}{filetype}")
#         elif filetype == "emf":
#             save_fig(f"{filename}svg", dpi=300, conv_svg_to_emf=True, verbose=True)
#         else:
#             raise ValueError(f"The filetype requested ({filetype}) is not yet implemented :) Please consult James or Rishi.")
#         plt.show()

#     def set_stats_scale(self, ax, scale_max, margin=0.1):
#         scale_max = scale_max + margin*scale_max
#         try:
#             ax.set_ylim([-1*scale_max, scale_max])
#         except Exception as e:
#             plotter_logger.debug(f"{e}")
#         pass

#     def plot_stats_dev(self, recordings, filename: Path):
#         analysis_logger.info(f"Verbose plotting of stats for {filename}")
#         fig, axs = plt.subplots(nrows = 5, ncols = 1, sharex="all", figsize=(15, 10), constrained_layout=True)
#         mean_depolarizations = np.asarray([recordings[paradigm].stats["depolarization"][0] for paradigm in recordings])
#         mean_hyperpolarizations = np.asarray([recordings[paradigm].stats["hyperpolarization"][0] for paradigm in recordings])
#         mean_Im = np.asarray([recordings[paradigm].stats["Imembrane"][0] for paradigm in recordings])
#         mean_Ileak = np.asarray([recordings[paradigm].stats["Ileakage"][0] for paradigm in recordings])
#         mean_Iactivation = np.asarray([recordings[paradigm].stats["Iactivation"][0] for paradigm in recordings])
#         mean_excitation = np.asarray([recordings[paradigm].stats["mean_excitation"][0] for paradigm in recordings])
#         mean_inhibition = np.asarray([recordings[paradigm].stats["mean_inhibition"][0] for paradigm in recordings])
#         net_excitation = np.asarray([recordings[paradigm].stats["net_excitation"][0] for paradigm in recordings])
#         net_inhibition = np.asarray([recordings[paradigm].stats["net_inhibition"][0] for paradigm in recordings])
#         spikes_per_stimulus_repetition = np.asarray([recordings[paradigm].stats["spikes_per_stimulus_repetition"][0] for paradigm in recordings])
#         paradigms = [paradigm for paradigm in recordings]
#         xlocations = np.asarray([x for x in range(len(paradigms))])
#         axs[0].bar(xlocations, mean_depolarizations, align='center', color="red")
#         axs[0].bar(xlocations, mean_hyperpolarizations, align='center', color="blue")
#         axs[0].axhline(0, color='grey', linewidth=0.8)
#         axs[0].set_ylabel("polarizations (V)")
#         scale_max = np.amax([np.amax(mean_depolarizations), np.amax(mean_hyperpolarizations)])
#         self.set_stats_scale(axs[0], scale_max, 0.1)
#         axs[1].bar(xlocations, mean_Im, align='center', color="black")
#         axs[1].axhline(0, color='grey', linewidth=0.8)
#         axs[1].set_ylabel("Im (A)")
#         scale_max = np.amax(mean_Im)
#         self.set_stats_scale(axs[1], scale_max, 0.1)
#         axs[2].bar(xlocations, mean_Iactivation, align='center', color="black")
#         axs[2].axhline(0, color='grey', linewidth=0.8)
#         axs[2].set_ylabel("Iact (A)")
#         scale_max = np.amax(mean_Iactivation)
#         self.set_stats_scale(axs[2], scale_max, 0.1)
#         axs[3].bar(xlocations, mean_Ileak, align='center', color="black")
#         axs[3].axhline(0, color='grey', linewidth=0.8)
#         axs[3].set_ylabel("Ileak (A)")
#         scale_max = np.amax(mean_Ileak)
#         self.set_stats_scale(axs[3], scale_max, 0.1)
#         axs[4].bar(xlocations, mean_excitation, align='center', color="red")
#         axs[4].bar(xlocations, -1*mean_inhibition, align='center', color="blue")
#         axs[4].axhline(0, color='grey', linewidth=0.8)
#         scale_max = np.amax([np.amax(mean_excitation), np.amax(mean_inhibition)])
#         self.set_stats_scale(axs[4], scale_max, 0.1)
#         ax = axs[4].twinx()
#         ax.plot(xlocations, spikes_per_stimulus_repetition, color='k', marker = 'o')
#         scale_max = np.amax(spikes_per_stimulus_repetition)
#         self.set_stats_scale(ax, scale_max, 0.1)
#         axs[4].set_ylabel("G (S)")
#         analysis_logger.info(f"Saving verbose plotting of stats for {filename}")
#         plt.savefig(str(filename)+f"_dev_stats.png")
#         pass

#     def run(self, current_clamps: Optional[List[float] | List[List[float]]]=None, filetype: str="png") -> None:
        for i in range(len(self.filepaths)):
            ccs = (current_clamps if current_clamps is None or isinstance(current_clamps[i], float) else current_clamps[i])
            reader = XLReader(filepath)
            recordings = {}
            for _, paradigm in enumerate(reader.get_paradigms()):
                recordings[paradigm] = WholeCellRecording(
                    reader.get_paradigm_data(paradigm), 
                    reader.get_paradigm_parameters(paradigm),
                    filters=lpfs,
                    current_clamps=current_clamps
                )
            recordings, overall_stats = self.estimation_without_optim_activation_potential(recordings)
            result_filename = os.path.join(self.output_path, f"{os.path.splitext(os.path.basename(filepath))[0]}_analyzed")
            # self.write_to_excel(f"{result_filename}.xlsx", recordings, stats)
            self.plot_dev(recordings, result_filename, filetype=filetype, current_clamps=ccs)
            self.plot_stats_dev(recordings, result_filename)