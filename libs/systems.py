# Array Handling
import pandas as pd
import numpy as np

# Optimization
from numba import njit

from scipy.optimize import minimize, Bounds

# Signal Processing
from scipy.signal import butter, buttord, sosfiltfilt
from scipy.ndimage import gaussian_filter1d

# Plotting & Graphics
import matplotlib.pyplot as plt
from pyhelpers.store import save_fig

# Local
from libs.readers import XLReader, AnalyzerCfg

# OS
from pathlib import Path

# Annotation
from typing import Dict, List, Tuple


class LowPassFilter:
    def __init__(self, passband: float, stopband: float, attenuation: float, ripple: float) -> None:
        assert passband < stopband, f"For low pass filter stopband cannot be less than passband."
        self.filter_design = lambda fs: butter(*buttord(passband, stopband, ripple, attenuation, fs=fs), output="sos", fs=fs)

    def propagate(self, raw_signal: np.ndarray, fs: float) -> np.ndarray:
        return sosfiltfilt(self.filter_design(fs), raw_signal)
    

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
    def __init__(self, parameters: pd.DataFrame, stimuli: Dict[str, pd.DataFrame], filters: Dict[str, LowPassFilter]) -> None:
        assert len(stimuli) > 0
        assert set(filters.keys()) == {"Vm", "Im", "g", "Eeff"}

        self.filters: Dict[str, LowPassFilter] = filters
        
        self.Cm: float = parameters["Cm"][0]    # units: Farads
        self.Et_measured: float = parameters["Et"][0]

        example_t: pd.Series = list(stimuli.values())[0]["times"] # units: Seconds, shape: (Nsamples,)
        self.dt: float = example_t[1] - example_t[0] # time assumed sampled at constant interval; units: Seconds

        self.Ee: float = np.nan # units: Volts
        self.Ei: float = np.nan # units: Volts
        
        self.LI_dvdt_vs_Vm: float = np.nan

        self.stimuli: Dict[str, WholeCellStimulus] = {name: WholeCellStimulus(self, data) for name, data in stimuli.items()}
        
        self.Er, self.Rin = self._estimate_Er_Rin()
        self._gact_override: float = np.nan
        
        self.Eact: float = np.nan
        self.Eact, self._gact_override = self._estimate_Eact_and_gact()
        
        print(f"Parameter Estimates:\n\tEr: {self.Er*1e3:.1f} mV\n\tRin: {self.Rin*1e-9:.1f} Gohm\n\tEact: {self.Eact*1e3:.1f} mV\n\tgact: {self.gact*1e9:.1f} nS\n\tEt: {self.Et*1e3:.1f} mV")
        print(f"\tgl: {self.gl*1e9:.1f} nS")

    @property
    def Et(self) -> float:
        return (self.Eact * self.gact - self.Er * self.gl) / (self.gact - self.gl)

    @property
    def gl(self) -> float:
        return 1 / self.Rin # units: Siemens

    @property
    def gact(self) -> float:
        if not np.isnan(self._gact_override):
            return self._gact_override
        return self.gl * (self.max_Vm - self.Er) / (self.max_Vm - self.Eact)

    @property
    def max_Vm(self) -> float:
        max_Vms: List[float] = []
        for stimulus in self.stimuli.values():
            max_Vms.append(stimulus.Vm.max())
        return np.max(max_Vms)

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
            Eeff_pool.extend(stimulus.timeseries["Eeff filtered"])

        Ei_hat: float = float(np.nanquantile(Eeff_pool, 0.05)) # units: Volts
        Ee_hat: float = float(np.nanquantile(Eeff_pool, 0.95)) # units: Volts

        Ee_hat = max(Ee_hat, self.Eact)

        return Ee_hat, Ei_hat

    def run_analysis(self, verbose: bool=True):
        # --- Step 1: estimate reversal potentials from Eeff/gsyn stage ---
        if verbose: print("Estimating reversal potentials... ")

        for stimulus in self.stimuli.values():
            stimulus.calculate_target_Isyn()
            stimulus.estimate_Eeff()

        # Original reversal estimate (used by default downstream)
        self.Ee, self.Ei = self.estimate_Ee_Ei()

        if verbose: print(f"Estimated Reversals: Ee = {self.Ee*1e3:.1f} mV, Ei = {self.Ei*1e3:.1f} mV")

        # --- Step 2: estimate ge/gi and forward-predict Vm ---
        if verbose: print("Calculating synaptic conductances... ")

        for stimulus in self.stimuli.values():
            stimulus.estimate_ge_gi()
            stimulus.calculate_predicted_Im()

    """ Eact Optimization """
    def sos_negative_g(self, x) -> float:
        self.Eact = x
        self.run_analysis(verbose=False)
        total = 0.0
        for stimulus in self.stimuli.values():
            ge = stimulus.timeseries["ge filtered"].to_numpy()
            sse_ge = np.sum(np.square(np.clip(ge, -np.inf, 0)))
            gi = stimulus.timeseries["gi filtered"].to_numpy()
            sse_gi = np.sum(np.square(np.clip(gi, -np.inf, 0)))
            total += sse_gi + sse_ge
        print(x, total)
        return 1e10 * total

    def _estimate_Eact_and_gact(self) -> Tuple[float, float]:
        eps_v: float = 2e-3

        lo_E_init: float = self.Er + eps_v
        hi_E_init: float = self.max_Vm
        max_V_safe: float = self.max_Vm + eps_v

        def gact_upper(E: float) -> float:
            return float(self.gl) * (max_V_safe - self.Er) / (max_V_safe - E)

        def unpack(x: np.ndarray) -> Tuple[float, float]:
            E, alpha = float(x[0]), float(x[1])
            g = float(self.gl) + alpha * (gact_upper(E) - float(self.gl))
            return E, g

        def objective(x: np.ndarray) -> float:
            E, g = unpack(x)
            self.Eact = E
            self._gact_override = g
            self.run_analysis(verbose=False)
            total = 0.0
            for stimulus in self.stimuli.values():
                ge = stimulus.timeseries["ge filtered"].to_numpy()
                gi = stimulus.timeseries["gi filtered"].to_numpy()
                total += np.sum(np.square(np.minimum(ge, 0.0)))
                total += np.sum(np.square(np.minimum(gi, 0.0)))
            return total

        if self.max_Vm <= lo_E_init:
            return lo_E_init, float(self.gl) * 1.1

        n_epochs: int = 3    # number of refinement passes
        depth: int = 5       # grid points per dimension per pass
        shrink: float = 0.5  # how much to shrink the search window each epoch

        lo_E: float = lo_E_init
        hi_E: float = hi_E_init
        lo_a: float = 0.1
        hi_a: float = 1.0

        best_x: np.ndarray = np.array([0.5 * (lo_E + hi_E), 0.5])
        best_val: float = np.inf

        for epoch in range(n_epochs):
            Eact_grid: np.ndarray = np.linspace(lo_E, hi_E, depth)
            alpha_grid: np.ndarray = np.linspace(lo_a, hi_a, depth)

            for E in Eact_grid:
                for a in alpha_grid:
                    v = objective(np.array([E, a]))
                    if v < best_val:
                        best_val = v
                        best_x = np.array([E, a])

            # Shrink the search window around the current best, keeping it inside the original bounds
            E_half: float = 0.5 * shrink * (hi_E - lo_E)
            a_half: float = 0.5 * shrink * (hi_a - lo_a)
            lo_E = max(lo_E_init, best_x[0] - E_half)
            hi_E = min(hi_E_init, best_x[0] + E_half)
            lo_a = max(0.0, best_x[1] - a_half)
            hi_a = min(1.0, best_x[1] + a_half)

        E_opt, g_opt = unpack(best_x)
        return E_opt, g_opt


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

    def calculate_target_Isyn(self) -> None:
        Vm: np.ndarray = self.Vm        # units: Volts, shape: [Nclamps, Nsamples]
        
        Iinj: np.ndarray = self.Iinj    # units: Amperes, shape: [Nclamps, 1]

        gl: float = self.recording.gl   # units: Siemens
        Cm: float = self.recording.Cm   # units: Faradays
        Er: float = self.recording.Er   # units: Volts

        dt: float = self.recording.dt # units: Seconds
        fs: float = 1 / dt # units: Herz

        Vm_filtered: np.ndarray = self.recording.filters["Vm"].propagate(Vm, fs) # units: Volts, shape: [Nclamps, Nsamples]

        Im: np.ndarray = Cm * np.gradient(Vm, dt, axis=-1) # units: Amperes, shape: [Nclamps, Nsamples]
        Im_filtered: np.ndarray = self.recording.filters["Im"].propagate(Im, fs) # units: Amperes, shape: [Nclamps, Nsamples]

        Il: np.ndarray = gl * (Er - Vm_filtered) # units: Amperes, shape: [Nclamps, Nsamples]

        Iact: np.ndarray = np.where(
            Vm_filtered > self.recording.Eact,
            self.recording.gact * (self.recording.Eact - Vm_filtered),
            0.0,
        )

        target_Isyn = -Im + Iinj + Il - Iact # units: Amperes, shape: [Nclamps, Nsamples]
        
        self.timeseries["Vm filtered"] = Vm_filtered.T.tolist()
        self.timeseries["Im"] = Im.T.tolist()
        self.timeseries["Im filtered"] = Im_filtered.T.tolist()
        self.timeseries["Il"] = Il.T.tolist()
        self.timeseries["Iact"] = Iact.T.tolist()
        self.timeseries["target Isyn"] = target_Isyn.T.tolist()

    def estimate_Eeff(self) -> None:
        # Pull data (faster than .tolist() if these are arrays-in-cells, but keep if needed)
        Vm_filtered: np.ndarray = np.stack(self.timeseries["Vm filtered"]).astype(np.float64).T  # type: ignore , units: Volts, shape: [Nclamps, Nsamples]
        target_Isyn: np.ndarray = np.stack(self.timeseries["target Isyn"]).astype(np.float64).T  # type: ignore , units: Amperes, shape: [Nclamps, Nsamples]

        # Means per bin
        Vm_filtered_mean: np.ndarray = Vm_filtered.mean(axis=0)  # units: Volts, shape: [Nsamples,]
        target_Isyn_mean: np.ndarray = target_Isyn.mean(axis=0)  # units: Amperes, shape: [Nsamples,]

        # Centered
        centered_integral_Vm: np.ndarray = Vm_filtered - Vm_filtered_mean[None, :]  # units: Volts, shape: [Nclamps, Nsamples]
        centered_target_Isyn: np.ndarray = target_Isyn - target_Isyn_mean[None, :]  # units: Amperes, shape: [Nclamps, Nsamples]

        # Regression slope b and intercept a for each bin
        denom = np.sum(centered_integral_Vm * centered_integral_Vm, axis=0)  # var * (Nclamps-1) up to scale
        numer = np.sum(centered_integral_Vm * centered_target_Isyn, axis=0)

        # Handle degenerate bins where Phi has no variation across clamps
        eps = np.finfo(np.float64).tiny
        b = np.where(np.abs(denom) > eps, numer / denom, np.nan)  # slope (Nbins,)
        a = target_Isyn_mean - b * Vm_filtered_mean  # intercept

        # Your derived params
        gsyn = -b  # Siemens
        # Avoid divide-by-zero when gsyn ~ 0
        gsyn_safe = np.where(np.abs(gsyn) > 1e-10, gsyn, np.nan)

        Eeff = a / gsyn_safe  # Volts

        # Diagnostics per bin
        # Your SSE formula: sum (Q - gsyn*(Eeff - Phi))^2
        # We can compute predicted Q directly from a + b*Phi (same fit)
        I_hat = a[None, :] + b[None, :] * Vm_filtered
        resid = target_Isyn - I_hat
        sse = np.sum(resid * resid, axis=0)

        # Proper per-bin R^2: 1 - SSE / SST, SST = sum (Q - mean(Q))^2 within the bin
        sst = np.sum((target_Isyn - target_Isyn_mean[None, :]) ** 2, axis=0)
        r2 = np.where(sst > 0, 1.0 - (sse / sst), np.nan)

        # Filter Eeff, handling NaNs by linear interpolation before filtering
        # and restoring NaNs in originally undefined regions afterward.
        nan_mask = np.isnan(Eeff)
        if nan_mask.all():
            Eeff_filtered = np.full_like(Eeff, np.nan)
        else:
            if nan_mask.any():
                idx = np.arange(Eeff.size)
                Eeff_clean = Eeff.copy()
                Eeff_clean[nan_mask] = np.interp(idx[nan_mask], idx[~nan_mask], Eeff[~nan_mask])
            else:
                Eeff_clean = Eeff

            fs = 1 / self.recording.dt
            Eeff_filtered = self.recording.filters["Eeff"].propagate(Eeff_clean[None, :], fs)[0]
            # Restore NaNs where Eeff was originally undefined so the filter output doesn't
            # invent values in regions with no meaningful synaptic conductance.
            Eeff_filtered[nan_mask] = np.nan

        # Store
        self.timeseries["Eeff"] = Eeff
        self.timeseries["Eeff filtered"] = Eeff_filtered
        self.timeseries["gsyn"] = gsyn_safe
        self.timeseries["SSE least-squares Qsyn"] = sse
        self.timeseries["r2 least-squares Qsyn"] = r2

    def estimate_ge_gi(self) -> None:
        Vm_filtered: np.ndarray = np.stack(self.timeseries["Vm filtered"]).astype(np.float64).T  # type: ignore , units: Volts, shape: [Nclamps, Nsamples]
        target_Isyn: np.ndarray = np.stack(self.timeseries["target Isyn"]).astype(np.float64).T  # type: ignore , units: Amperes, shape: [Nclamps, Nsamples]

        Ee: float = self.recording.Ee  # units: Volts
        Ei: float = self.recording.Ei  # units: Volts

        driving_force: np.ndarray = np.stack([Vm_filtered - Ee, Vm_filtered - Ei], axis=0)  # units: Volts, shape: [2, Nclamps, Nsamples]

        XtX: np.ndarray = np.einsum("ijn,kjn->nik", driving_force, driving_force)  # shape: [Nsamples, 2, 2]
        Xty: np.ndarray = np.einsum("ijn,jn->ni", driving_force, target_Isyn)  # shape: [Nsamples, 2]

        try:
            conductances = np.linalg.solve(XtX, Xty[..., None])[..., 0]  # shape: [Nsamples, 2]
        except np.linalg.LinAlgError:
            print("Cannot use least-squares, attempting pseudoinverse...")
            conductances = (np.linalg.pinv(XtX) @ Xty[..., None])[..., 0]  # shape: [Nsamples, 2]

        ge: np.ndarray = conductances[:, 0]  # units: Siemens, shape: [Nsamples]
        gi: np.ndarray = conductances[:, 1]  # units: Siemens, shape: [Nsamples]

        fs: float = 1 / self.recording.dt  # units: Hertz
        ge_filtered: np.ndarray = self.recording.filters["g"].propagate(ge[None, :], fs)[0]  # units: Siemens, shape: [Nsamples]
        gi_filtered: np.ndarray = self.recording.filters["g"].propagate(gi[None, :], fs)[0]  # units: Siemens, shape: [Nsamples]

        self.timeseries["ge"] = ge.tolist()
        self.timeseries["gi"] = gi.tolist()
        self.timeseries["ge filtered"] = ge_filtered.tolist()
        self.timeseries["gi filtered"] = gi_filtered.tolist()
    
    def calculate_predicted_Im(self) -> None:
        """
        Compute predicted membrane current from the fitted conductance model,
        evaluated at the measured Vm. If the model is correct and fits well,
        predicted Im should match measured Im up to noise.

        Current balance:
            Cm·dVm/dt = Iinj - gl·(Vm - Er) - ge·(Vm - Ee) - gi·(Vm - Ei) - Iact(Vm)
        so
            Im_predicted = Cm·dVm/dt = Iinj - Il_model - Isyn_model - Iact_model
        where all RHS terms are evaluated at the measured Vm.
        """
        Iinj: np.ndarray = self.Iinj  # units: Amperes, shape: [Nclamps, 1]

        Ee: float = float(self.recording.Ee)        # units: Volts
        Ei: float = float(self.recording.Ei)        # units: Volts
        Er: float = float(self.recording.Er)        # units: Volts
        Eact: float = float(self.recording.Eact)    # units: Volts
        gl: float = float(self.recording.gl)        # units: Siemens
        gact: float = float(self.recording.gact)    # units: Siemens

        Vm_filtered: np.ndarray = np.stack(self.timeseries["Vm filtered"]).astype(np.float64).T  # type: ignore , units: Volts, shape: [Nclamps, Nsamples]

        ge: np.ndarray = np.stack(self.timeseries["ge"]).astype(np.float64).T  # type: ignore , units: Siemens, shape: [Nsamples,]
        gi: np.ndarray = np.stack(self.timeseries["gi"]).astype(np.float64).T  # type: ignore , units: Siemens, shape: [Nsamples,]

        # Model currents evaluated at measured Vm
        Il_model: np.ndarray = gl * (Vm_filtered - Er)                              # units: Amperes, shape: [Nclamps, Nsamples]
        Ie_model: np.ndarray = ge[None, :] * (Vm_filtered - Ee)                     # units: Amperes, shape: [Nclamps, Nsamples]
        Ii_model: np.ndarray = gi[None, :] * (Vm_filtered - Ei)                     # units: Amperes, shape: [Nclamps, Nsamples]
        Iact_model: np.ndarray = np.where(
            Vm_filtered > Eact,
            gact * (Eact - Vm_filtered),
            0.0,
        )                                                                            # units: Amperes, shape: [Nclamps, Nsamples]

        # Predicted membrane current from current balance
        Im_predicted: np.ndarray = Iinj - Il_model - Ie_model - Ii_model - Iact_model  # units: Amperes, shape: [Nclamps, Nsamples]
        Im_predicted_filtered: np.ndarray = self.recording.filters["Im"].propagate(Im_predicted, 1 / self.recording.dt)

        self.timeseries["Iact predicted"] = Iact_model.T.tolist()
        self.timeseries["Ie predicted"] = Ie_model.T.tolist()
        self.timeseries["Ii predicted"] = Ii_model.T.tolist()
        self.timeseries["Im predicted"] = Im_predicted.T.tolist()
        self.timeseries["Im predicted filtered"] = Im_predicted_filtered.T.tolist()

    def calculate_predicted_Vm(self) -> None:
        Iinj: np.ndarray = self.Iinj    # units: Amperes, shape: [Nclamps, 1]
        V0: np.ndarray = self.Vm[:, 0]  # units: Volts, shape: [Nclamps,]

        Ee: float = float(self.recording.Ee)        # units: Volts
        Ei: float = float(self.recording.Ei)        # units: Volts
        Er: float = float(self.recording.Er)        # units: Volts
        Eact: float = float(self.recording.Eact)    # units: Volts
        gl: float = float(self.recording.gl)        # units: Siemens
        gact: float = float(self.recording.gact)    # units: Siemens
        Cm: float = float(self.recording.Cm)        # units: Farads
        dt: float = float(self.recording.dt)        # units: Seconds

        ge: np.ndarray = np.stack(self.timeseries["ge"]).astype(np.float64).T  # type: ignore , units: Siemens, shape: [Nsamples,]
        gi: np.ndarray = np.stack(self.timeseries["gi"]).astype(np.float64).T  # type: ignore , units: Siemens, shape: [Nsamples,]

        # Passive components (do not depend on Vm)
        g_passive: np.ndarray = gl + ge + gi                                    # units: Siemens, shape: [Nsamples,]
        num_passive: np.ndarray = gl * Er + ge * Ee + gi * Ei + Iinj            # units: Amperes, shape: [Nclamps, Nsamples]

        pred_Vm: np.ndarray = self._vm_loop_with_iact(
            num_passive.astype(np.float32),
            g_passive.astype(np.float32),
            np.float32(gact),
            np.float32(Eact),
            np.float32(Cm),
            np.float32(dt),
            V0.astype(np.float32),
        )

        self.timeseries["predicted Vm"] = pred_Vm.T.tolist()

    @staticmethod
    @njit(cache=True)
    def _vm_loop_with_iact(
        num_passive: np.ndarray,    # gl*Er + ge*Ee + gi*Ei + Iinj, shape: [Nclamps, Nsamples]
        g_passive: np.ndarray,      # gl + ge + gi, shape: [Nsamples,]
        gact: np.float32,
        Eact: np.float32,
        Cm: np.float32,
        dt: np.float32,
        v0: np.ndarray,             # shape: [Nclamps,]
    ) -> np.ndarray:
        Nclamps, Nsamples = num_passive.shape
        Vm = np.empty_like(num_passive, dtype=np.float32)
        Vm[:, 0] = v0
        for i in range(1, Nsamples):
            for c in range(Nclamps):
                v_prev = Vm[c, i - 1]
                # Decide whether active current is on, based on previous Vm
                if v_prev > Eact:
                    g_tot = g_passive[i] + gact
                    v_inf = (num_passive[c, i] + gact * Eact) / g_tot
                else:
                    g_tot = g_passive[i]
                    v_inf = num_passive[c, i] / g_tot
                decay = np.exp(-g_tot * dt / Cm)
                Vm[c, i] = v_inf + (v_prev - v_inf) * decay
        return Vm


class Analyzer:
    def __init__(self, cfg: AnalyzerCfg):
        self.cfg: AnalyzerCfg = cfg

    def plot_timeseries(
        self,
        recording: WholeCellRecording,
        filename: Path,
        filetype: str = "png",
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
        nrows = 5
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

            filtered_Vm: np.ndarray = np.stack(stimulus.timeseries["Vm filtered"]).astype(np.float64).T #type: ignore
            Eeff: np.ndarray = stimulus.timeseries["Eeff filtered"].to_numpy(np.float64) 
            Im: np.ndarray = np.stack(stimulus.timeseries["Im filtered"]).astype(np.float64).T  # type: ignore
            ge: np.ndarray = stimulus.timeseries["ge filtered"].to_numpy(np.float64)
            gi: np.ndarray = stimulus.timeseries["gi filtered"].to_numpy(np.float64)
            
            axs[0, idx].set_title(paradigm)

            # ---- Row 0: Vm, Vss, Eact ----
            curves = axs[0, idx].plot(times, filtered_Vm.T)  # solid Vm traces
            colors = [line.get_color() for line in curves]

            Ess = recording.Er + recording.Rin * stimulus.Iinj
            for j in range(stimulus.Nclamps):
                axs[0, idx].plot([times[0], times[-1]], [Ess[j]] * 2, color="grey", ls="--")
            axs[0, idx].plot([times[0], times[-1]], [recording.Eact] * 2, color="red", ls="--")

            axs[0, idx].set_ylabel("Vm (V)")
            axs[0, idx].grid(True)

            # ---- Row 1: Eeff, Ee, Ei ----
            axs[1, idx].grid(True)
            axs[1, idx].plot(times, Eeff, c="black")
            axs[1, idx].axhline(recording.Er, linestyle="--", color="k", linewidth=1, label="Er" if idx == 0 else None)
            axs[1, idx].axhline(recording.Ee, linestyle="--", color="r", linewidth=1, label="Ee" if idx == 0 else None)
            axs[1, idx].axhline(recording.Ei, linestyle="--", color="b", linewidth=1, label="Ei" if idx == 0 else None)
            
            # ---- Row 2: Im ----
            for j in range(stimulus.Nclamps):
                axs[2, idx].plot(times, Im[j, :], color=colors[j])
            axs[2, idx].plot([times[0], times[-1]], [0, 0], color="grey", ls="--")
            axs[2, idx].set_ylabel("Im (A)")
            axs[2, idx].grid(True)

            # ---- Row 3: ge, gi ----
            axs[3, idx].plot(times, ge, c="r", label="ge")
            axs[3, idx].plot(times, gi, c="b", label="gi")

            axs[3, idx].plot(times, times * 0, "--k", linewidth=1)
            axs[3, idx].set_ylabel("G (S)")
            axs[3, idx].grid(True)

            # ---- Row 4: Iact predicted ----
            Iact_pred: np.ndarray = np.stack(stimulus.timeseries["Iact predicted"]).astype(np.float64).T  # type: ignore
            for j in range(stimulus.Nclamps):
                axs[4, idx].plot(times, Iact_pred[j, :], color=colors[j])
            axs[4, idx].plot(times, times * 0, "k--")
            axs[4, idx].set_ylabel("Iact (A)")
            axs[4, idx].grid(True)

            axs[4, idx].set_xlabel("time (s)")

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
        filters: Dict[str, LowPassFilter] = {
            "Im": LowPassFilter(30.0, 80.0, 80.0, 3),
            "Vm": LowPassFilter(30.0, 80.0, 80.0, 3),
            "g": LowPassFilter(30.0, 80.0, 80.0, 3),
            "Eeff": LowPassFilter(30.0, 80.0, 80.0, 3)
        }

        n_files: int = len(self.cfg.paths_to_spreadsheets)

        for i in range(n_files):
            path_to_spreadsheet: Path = self.cfg.paths_to_spreadsheets[i]
            print(f"Collecting data from {path_to_spreadsheet}")
            rdr: XLReader = XLReader(path_to_spreadsheet)
            
            stimuli: Dict[str, pd.DataFrame] = {paradigm:rdr.get_paradigm_data(paradigm) for paradigm in rdr.get_paradigms()}
            recording: WholeCellRecording = WholeCellRecording(rdr.get_paradigm_parameters(rdr.get_paradigms()[0]), stimuli, filters)
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