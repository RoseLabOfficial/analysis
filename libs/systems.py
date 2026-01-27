import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from pyhelpers.store import save_fig
from pathlib import Path
from scipy.optimize import minimize, Bounds
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter

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
        
    def run_analysis(self, verbose: bool=True):
        if verbose: print("Estimating reversal potentials... ")
        for stimulus in self.stimuli.values():
            stimulus.calculate_target_Qsyn()
            stimulus.calculate_target_Qsyn_nonlinear()
            stimulus.estimate_Eeff()
        self.Ee, self.Ei = self.estimate_Ee_Ei()
        if verbose: print(f"Estimated Reversals: Ee = {self.Ee*1e3:.1f} mV, Ei = {self.Ei*1e3:.1f} mV")

        if verbose: print("Calculating synaptic conductances... ")
        for stimulus in self.stimuli.values():
            stimulus.estimate_ge_gi(nonlinear=True, constrained=True)
            stimulus.calculate_predicted_Vm()
        
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
        self.Nbins: int = np.ceil(self.Nsamples / self.recording.bin_nsamples)

        self.timeseries: pd.DataFrame = pd.DataFrame({"times": self.times})

        l_index: np.ndarray = np.arange(0, self.Nsamples, self.recording.bin_nsamples)
        r_index: np.ndarray = np.minimum(l_index + self.recording.bin_nsamples, self.Nsamples - 1)
        l_time: np.ndarray = l_index * self.recording.dt
        r_time: np.ndarray = r_index * self.recording.dt
        bin_duration: np.ndarray = r_time - l_time
        self.binned_timeseries: pd.DataFrame = pd.DataFrame({
            "left index": l_index,
            "right index": r_index, 
            "left time": l_time, 
            "right time": r_time,
            "bin duration": bin_duration
        })
    
    def _cumtrapz_prefix_integral(self, arr: np.ndarray) -> np.ndarray:
        pref: np.ndarray = np.zeros((self.Nclamps, self.Nsamples + 1), dtype=arr.dtype) # shape: (Nclamps, Nsamples + 1)
        area: np.ndarray = 0.5 * (arr[:, 1:] + arr[:, :-1]) * self.recording.dt         # shape: (Nclamps, Nsamples)
        pref[:, 2:] = np.cumsum(area, axis=1)
        return pref

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

    def estimate_ge_gi(
        self,
        nonlinear: bool = False,
        constrained: bool = True,
        fit_ge: bool = False,
        fit_gi: bool = True,
    ) -> None:
        """
        Estimate ge and gi per bin.

        fit_ge / fit_gi control whether each conductance is fit.
        - fit_ge=True,  fit_gi=True  : fit both (your current behavior)
        - fit_ge=True,  fit_gi=False : excitation-only (gi == 0 for all bins)
        - fit_ge=False, fit_gi=True  : inhibition-only (ge == 0 for all bins)
        - fit_ge=False, fit_gi=False : both forced to 0 (degenerate)
        """

        bin_duration: np.ndarray = self.binned_timeseries["bin duration"].to_numpy(dtype=np.float64)
        integral_Vm: np.ndarray = np.stack(self.binned_timeseries["integral Vm"].to_numpy()).astype(np.float64).T  # type: ignore

        if nonlinear:
            target_Qsyn = np.stack(self.binned_timeseries["target Qsyn nonlinear"].to_numpy()).astype(np.float64).T  # type: ignore
        else:
            target_Qsyn = np.stack(self.binned_timeseries["target Qsyn"].to_numpy()).astype(np.float64).T  # type: ignore

        Ee: float = self.recording.Ee
        Ei: float = self.recording.Ei

        Phie: np.ndarray = Ee * bin_duration[np.newaxis, :] - integral_Vm
        Phii: np.ndarray = Ei * bin_duration[np.newaxis, :] - integral_Vm

        # dot products per bin (sum over clamps axis=0)
        eTe: np.ndarray = np.sum(Phie * Phie, axis=0)
        iTi: np.ndarray = np.sum(Phii * Phii, axis=0)
        eTi: np.ndarray = np.sum(Phie * Phii, axis=0)
        eTq: np.ndarray = np.sum(Phie * target_Qsyn, axis=0)
        iTq: np.ndarray = np.sum(Phii * target_Qsyn, axis=0)
        qTq: np.ndarray = np.sum(target_Qsyn * target_Qsyn, axis=0)

        # --- handle "excitation-only" / "inhibition-only" modes first ---
        if fit_ge and not fit_gi:
            # gi forced to 0
            ge = np.where(eTe > 0, eTq / eTe, 0.0)
            if constrained:
                ge = np.maximum(ge, 0.0)
            gi = np.zeros_like(ge)
            r2 = qTq - 2.0 * ge * eTq + (ge * ge) * eTe

        elif fit_gi and not fit_ge:
            # ge forced to 0
            gi = np.where(iTi > 0, iTq / iTi, 0.0)
            if constrained:
                gi = np.maximum(gi, 0.0)
            ge = np.zeros_like(gi)
            r2 = qTq - 2.0 * gi * iTq + (gi * gi) * iTi

        elif (not fit_ge) and (not fit_gi):
            # both forced to 0 (degenerate)
            ge = np.zeros_like(eTe)
            gi = np.zeros_like(eTe)
            r2 = qTq.copy()

        else:
            # --- fit both (your original logic) ---
            det: np.ndarray = eTe * iTi - eTi * eTi
            safe_det: np.ndarray = np.where(np.abs(det) > np.finfo(np.float64).tiny, det, np.nan)

            ge_u: np.ndarray = ( iTi * eTq - eTi * iTq) / safe_det
            gi_u: np.ndarray = (-eTi * eTq + eTe * iTq) / safe_det

            r2_u = (
                qTq
                - 2.0 * ge_u * eTq
                - 2.0 * gi_u * iTq
                + (ge_u * ge_u) * eTe
                + 2.0 * ge_u * gi_u * eTi
                + (gi_u * gi_u) * iTi
            )

            if constrained:
                # boundary candidates (nonnegative)
                ge_a = np.where(eTe > 0, eTq / eTe, 0.0)
                ge_a = np.maximum(ge_a, 0.0)
                r2_a = qTq - 2.0 * ge_a * eTq + (ge_a * ge_a) * eTe

                gi_b = np.where(iTi > 0, iTq / iTi, 0.0)
                gi_b = np.maximum(gi_b, 0.0)
                r2_b = qTq - 2.0 * gi_b * iTq + (gi_b * gi_b) * iTi

                feasible_u = (
                    (ge_u >= 0.0) & (gi_u >= 0.0) &
                    np.isfinite(ge_u) & np.isfinite(gi_u) &
                    np.isfinite(r2_u)
                )

                # pick best
                ge = ge_a.copy()
                gi = np.zeros_like(ge)
                r2 = r2_a.copy()

                pick_b = r2_b < r2
                ge[pick_b] = 0.0
                gi[pick_b] = gi_b[pick_b]
                r2[pick_b] = r2_b[pick_b]

                pick_u = feasible_u & (r2_u < r2)
                ge[pick_u] = ge_u[pick_u]
                gi[pick_u] = gi_u[pick_u]
                r2[pick_u] = r2_u[pick_u]

            else:
                ge = ge_u.copy()
                gi = gi_u.copy()
                r2 = r2_u.copy()

        resnorm = np.sqrt(np.maximum(r2, 0.0))

        # cond(XtX) from eigenvalues of 2x2 gram matrix
        tr = eTe + iTi
        disc = np.sqrt((eTe - iTi) ** 2 + 4.0 * (eTi ** 2))
        lam1 = 0.5 * (tr + disc)
        lam2 = 0.5 * (tr - disc)
        cond = np.where(lam2 > 0, lam1 / lam2, np.inf)

        suffix = f"{' nonlinear' if nonlinear else ''}{' constrained' if constrained else ''}"
        self.binned_timeseries[f"ge{suffix}"] = ge
        self.binned_timeseries[f"gi{suffix}"] = gi
        self.binned_timeseries[f"ge/gi residual norm{suffix}"] = resnorm
        self.binned_timeseries[f"ge/gi matrix conditioning{suffix}"] = cond
    
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
            bin_times: np.ndarray = stimulus.binned_timeseries["left time"].to_numpy(np.float64)

            bin_durations: np.ndarray = stimulus.binned_timeseries["bin duration"].to_numpy(np.float64)

            Vpred: np.ndarray = np.stack(stimulus.timeseries["predicted Vm"].to_numpy()).astype(np.float64).T #type: ignore
            integral_Iact: np.ndarray = np.stack(stimulus.binned_timeseries["integral Iact"].to_numpy()).astype(np.float64).T #type: ignore

            Eeff: np.ndarray = stimulus.binned_timeseries["Eeff"].to_numpy(np.float64) 
            
            ge = stimulus.binned_timeseries["ge nonlinear constrained"].to_numpy(np.float64)
            gi = stimulus.binned_timeseries["gi nonlinear constrained"].to_numpy(np.float64)
            
            # diagnostics
            resnorm = stimulus.binned_timeseries["ge/gi residual norm nonlinear constrained"].to_numpy()
            cond = stimulus.binned_timeseries["ge/gi matrix conditioning nonlinear constrained"].to_numpy()

            axs[0, idx].set_title(paradigm)

            r2_iv = stimulus.timeseries.get("r2_dvdt_vs_Vm", None)
            LI_iv = stimulus.timeseries.get("LI_dvdt_vs_Vm", None)
            if r2_iv is not None:
                r2_iv = np.asarray(r2_iv, dtype=np.float64)
                axs[6, idx].plot(times, r2_iv, linewidth=1)  # or pick a new row
                if LI_iv is not None:
                    LI_val = float(np.asarray(LI_iv, dtype=np.float64)[0])
                    axs[0, idx].set_title(f"{paradigm}  LI={LI_val:.3f}")

            # ---- Row 0: Vm + Vpred (same colors) ----
            curves = axs[0, idx].plot(times, stimulus.Vm.T)  # solid Vm traces
            colors = [line.get_color() for line in curves]

            Ess = recording.Er + recording.Rin * stimulus.Iinj
            for j in range(stimulus.Nclamps):
                axs[0, idx].plot(times, Vpred[j, :], linestyle=":", color=colors[j])  # dotted Vpred
                axs[0, idx].plot([times[0], times[-1]], [Ess[j]] * 2, color="grey", ls="--")
            axs[0, idx].plot([times[0], times[-1]], [recording.Eact] * 2, color="red", ls="--")

            axs[0, idx].set_ylabel("Vm (V)")
            axs[0, idx].grid(True)
            # if idx == 0:
            #     axs[0, idx].legend(loc="upper right")

            # ---- Row 1: Eeff ----
            axs[1, idx].grid(True)
            axs[1, idx].plot(bin_times, Eeff, c="black")
            axs[1, idx].axhline(recording.Er, linestyle="--", color="k", linewidth=1, label="Er" if idx == 0 else None)
            axs[1, idx].axhline(recording.Ee, linestyle="--", color="r", linewidth=1, label="Ee" if idx == 0 else None)
            axs[1, idx].axhline(recording.Ei, linestyle="--", color="b", linewidth=1, label="Ei" if idx == 0 else None)

            axs[2, idx].grid(True)
            axs[2, idx].plot(bin_times, bin_times * 0, "k--")
            axs[2, idx].plot(bin_times, integral_Iact.T / bin_durations[:, np.newaxis])

            # ---- Row 3: conductances ----
            axs[3, idx].plot(bin_times, ge, c="r", label="ge")
            axs[3, idx].plot(bin_times, gi, c="b", label="gi")

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
