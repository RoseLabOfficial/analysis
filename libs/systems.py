# Array Handling
import os
import json
import pandas as pd
import numpy as np

# Optimization
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy

# Signal Processing
from scipy.signal import butter, buttord, sosfiltfilt
from scipy.ndimage import gaussian_filter1d
from scipy.stats import spearmanr

from sklearn.mixture import GaussianMixture

# Plotting & Graphics
import matplotlib.pyplot as plt
from pyhelpers.store import save_fig

# Local
from libs.readers import (
    XLReader,
    AnalyzerCfg,
    FilterCfg,
    RecordingCfg,
    ErRinCfg,
    EeEiCfg,
    EtOptCfg,
    NumericsCfg,
)

# OS
from pathlib import Path

# Annotation
from typing import Dict, List, Tuple, Optional, Any


class LowPassFilter:
    def __init__(self, passband: float, stopband: float, attenuation: float, ripple: float) -> None:
        assert passband < stopband, f"For low pass filter stopband cannot be less than passband."
        self.passband = passband
        self.stopband = stopband
        self.attenuation = attenuation
        self.ripple = ripple

    @classmethod
    def from_cfg(cls, cfg: "FilterCfg") -> "LowPassFilter":
        return cls(
            passband=cfg.passband_hz,
            stopband=cfg.stopband_hz,
            attenuation=cfg.attenuation_db,
            ripple=cfg.ripple_db,
        )

    def _design(self, fs: float):
        return butter(
            *buttord(self.passband, self.stopband, self.ripple, self.attenuation, fs=fs),
            output="sos",
            fs=fs,
        )

    def propagate(self, raw_signal: np.ndarray, fs: float) -> np.ndarray:
        return sosfiltfilt(self._design(fs), raw_signal)

def weighted_quantile(x, q, w=None, axis=-1):
    """
    Weighted quantiles of `x` at quantiles `q` in [0,1], along `axis`.
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

    x = np.moveaxis(x, axis, -1)
    w = np.moveaxis(w, axis, -1)

    *batch, n = x.shape
    m = int(np.prod(batch)) if batch else 1
    x2 = x.reshape(m, n)
    w2 = w.reshape(m, n)

    idx = np.argsort(x2, axis=1)
    xs = np.take_along_axis(x2, idx, axis=1)
    ws = np.take_along_axis(w2, idx, axis=1)

    cw = np.cumsum(ws, axis=1)
    total = cw[:, -1]
    if np.any(total <= 0):
        raise ValueError("each slice must have positive total weight.")

    cdf = cw / total[:, None]

    qv = q.ravel()
    mask = cdf[:, None, :] >= qv[None, :, None]
    any_true = mask.any(axis=2)
    hi = mask.argmax(axis=2)
    hi = np.where(any_true, hi, n - 1)

    lo = np.clip(hi - 1, 0, n - 1)

    x_hi = np.take_along_axis(xs, hi, axis=1)
    x_lo = np.take_along_axis(xs, lo, axis=1)

    c_hi = np.take_along_axis(cdf, hi, axis=1)
    c_lo = np.where(hi > 0, np.take_along_axis(cdf, lo, axis=1), 0.0)

    denom = (c_hi - c_lo)
    t = np.where(denom > 0, (qv[None, :] - c_lo) / denom, 0.0)
    t = np.clip(t, 0.0, 1.0)

    out = x_lo + t * (x_hi - x_lo)

    k = qv.size
    out = out.reshape((*batch, k))

    if q.ndim == 0:
        out = out[..., 0]

    return out


# ---------------------------------------------------------------------------
# Worker function for grid-search parallelism. Module-level so it pickles.
# ---------------------------------------------------------------------------
# Parameter vector layout: theta = (Et_base, delta_1, ..., delta_{N-1}, x_beta)
# where Iinj-sorted Et is reconstructed as
#     Et_sorted[k] = Et_base + sum_{j<k} delta_j,    delta_j >= 0
# enforcing monotonic non-decreasing Et with Iinj (more depolarizing -> higher Et).

def _nan_to_none(x: float) -> Optional[float]:
    """JSON has no NaN; encode as None on save."""
    return None if (x is None or (isinstance(x, float) and np.isnan(x))) else float(x)


def _none_to_nan(x: Optional[float]) -> float:
    """Inverse of _nan_to_none."""
    return float("nan") if x is None else float(x)


def _theta_to_Et_dict(theta: np.ndarray, Iinj_sorted: np.ndarray) -> Tuple[Dict[float, float], float]:
    """Decode a parameter vector into (Et_by_Iinj, x_beta)."""
    N: int = Iinj_sorted.size
    assert theta.size == N + 1, f"theta has {theta.size} elems, expected {N+1}"

    Et_base: float = float(theta[0])
    deltas: np.ndarray = np.maximum(theta[1:N], 0.0)        # enforce >=0 even if optimizer drifts
    x_beta: float = float(theta[N])

    Et_sorted: np.ndarray = Et_base + np.concatenate([[0.0], np.cumsum(deltas)])
    Et_by_Iinj: Dict[float, float] = {float(I): float(Et_sorted[k]) for k, I in enumerate(Iinj_sorted)}
    return Et_by_Iinj, x_beta


def _evaluate_active_params(args: Tuple[np.ndarray, "WholeCellRecording"]) -> Tuple[np.ndarray, float]:
    """
    Worker: takes (theta, recording_copy), mutates the copy, runs analysis,
    returns the SOS-negative-Δgsyn loss summed over all clamps and stimuli.

    NOTE (TODO item 2): This loss penalizes any negative Δgsyn excursion.
    Under the corrected Δg framing, negative Δg has two sources:
      (a) model misspecification (active currents, Cm error, etc.) -- the
          α/β correction *should* absorb these
      (b) real disinhibition / withdrawal of tonic input -- the α/β correction
          should NOT absorb these
    The current loss conflates (a) and (b). Replacement options under
    discussion: voltage-dependence prior, timescale prior, or pharmacology-
    anchored joint fit. See TODO.md item 2.
    """
    theta, rec = args
    Iinj_sorted: np.ndarray = rec._Iinj_sorted_for_opt
    Et_by_Iinj, x_beta = _theta_to_Et_dict(theta, Iinj_sorted)

    rec.Et = Et_by_Iinj
    rec.x_beta = x_beta
    rec.run_analysis(verbose=False, complete=False)

    total: float = 0.0
    for stim in rec.stimuli.values():
        dgsyn: np.ndarray = stim.timeseries["dgsyn filtered"].to_numpy()
        total += float(np.sum(np.square(np.minimum(dgsyn, 0.0))))

    return theta, total


class WholeCellRecording:
    def __init__(
        self,
        parameters: pd.DataFrame,
        stimuli: Dict[str, pd.DataFrame],
        filters: Dict[str, LowPassFilter],
        n_workers: int,
        recording_cfg: RecordingCfg,
        Er_Rin_cfg: ErRinCfg,
        Ee_Ei_cfg: EeEiCfg,
        Et_opt_cfg: EtOptCfg,
        numerics_cfg: NumericsCfg,
        cached_params: Optional[Dict[str, Any]] = None,
    ) -> None:
        assert len(stimuli) > 0
        assert set(filters.keys()) == {"Vm", "Im", "dg", "Eeff", "dgsyn"}

        self.filters: Dict[str, LowPassFilter] = filters

        # Cfg objects -- stored so workers and stimuli can read them after deepcopy.
        # recording_cfg is read by WholeCellStimulus to apply the LJP correction
        # at ingestion (the only voltage-shifting step in the pipeline).
        self.recording_cfg: RecordingCfg = recording_cfg
        self.Er_Rin_cfg: ErRinCfg = Er_Rin_cfg
        self.Ee_Ei_cfg: EeEiCfg = Ee_Ei_cfg
        self.Et_opt_cfg: EtOptCfg = Et_opt_cfg
        self.numerics_cfg: NumericsCfg = numerics_cfg

        self.Cm: float = parameters["Cm"][0]
        # User-supplied Et from spreadsheet is in raw (uncorrected) volts;
        # apply LJP correction at ingestion so all internal voltages share
        # the same reference frame as the corrected Vm traces.
        # NOTE: Eact and Ess in the parameter sheet are NOT used downstream
        # (only sanity-checked at read in XLReader.get_paradigm_parameters);
        # if they ever start being used, they need the same correction.
        V_LJP: float = recording_cfg.liquid_junction_potential_volts
        self.Et_measured: float = parameters["Et"][0] - V_LJP

        example_t: pd.Series = list(stimuli.values())[0]["times"]
        self.dt: float = example_t[1] - example_t[0]

        self.Ee: float = np.nan
        self.Ei: float = np.nan
        self.LI_dvdt_vs_Vm: float = np.nan

        self.stimuli: Dict[str, WholeCellStimulus] = {
            name: WholeCellStimulus(self, data, paradigm=name) for name, data in stimuli.items()
        }

        # Step 1: precompute Vm_filtered, Im, Im_filtered for every stimulus.
        # These are stable across the optimization loop.
        for stim in self.stimuli.values():
            stim.precompute_filtered_traces()

        # Per-stimulus raw mode estimates (paradigm -> Iinj -> Vss_hat).
        # Filled by _estimate_Vss_per_stimulus.
        self.Vss_per_stimulus: Dict[str, Dict[float, float]] = {}
        # Per-stimulus single-stimulus (Er, Rin) regression fits, used as input
        # to drift clustering. paradigm -> (Er_hat, Rin_hat).
        self.Er_Rin_per_stimulus: Dict[str, Tuple[float, float]] = {}
        # Cluster assignment: paradigm -> cluster index (0-based).
        self.cluster_assignment: Dict[str, int] = {}
        # Cluster-pooled (Er, Rin) estimates: cluster index -> value.
        self.Er_by_cluster: Dict[int, float] = {}
        self.Rin_by_cluster: Dict[int, float] = {}
        # Cluster-pooled Vss: paradigm -> Iinj -> Vss_hat (each paradigm gets its
        # cluster's pooled estimate at each of its Iinj levels).
        self.Vss: Dict[str, Dict[float, float]] = {}

        # Per-clamp Et (one value per unique Iinj), shared scalar x_beta in [0, x_beta_max]
        self.Et: Dict[float, float] = {}
        self.x_beta: float = np.nan
        self.opt_objective: float = float("nan")

        # Canonical Iinj ordering (ascending: most-hyperpolarizing first).
        # Used by the optimizer constraint and per-clamp max_Vm bounds.
        self._Iinj_sorted_for_opt: np.ndarray = self._collect_unique_Iinj_sorted()
        self._max_Vm_by_Iinj: Dict[float, float] = self._collect_max_Vm_by_Iinj()

        if cached_params is None:
            # Step 2: resting-state pipeline (uses Vm_filtered + Im_filtered for the
            # weighted-mode estimator; produces Vss, Er_by_cluster, Rin_by_cluster).
            self._estimate_resting_state()

            # Step 3: precompute Il for every stimulus, now that Vss and gl(paradigm)
            # are known. Also stable across the optimization loop.
            for stim in self.stimuli.values():
                stim.precompute_Il()

            # Step 4: optimization loop. Only Iact + target_Isyn are recomputed per
            # iteration; everything else above is reused.
            self.Et, self.x_beta, self.opt_objective = self._estimate_Et_beta(n_workers=n_workers)
        else:
            # Cache hit: skip the expensive steps; load Vss / Er_by_cluster /
            # Rin_by_cluster / Et / x_beta from the JSON. Then precompute Il so
            # downstream analysis (run_analysis -> calculate_target_Isyn -> ...)
            # has what it needs.
            self._load_cached_params(cached_params)
            for stim in self.stimuli.values():
                stim.precompute_Il()

        # Pretty-print summary
        source_label = "loaded from cache" if cached_params is not None else "computed"
        Et_str = ", ".join(f"{self.Et[float(I)]*1e3:.1f}" for I in self._Iinj_sorted_for_opt)
        cluster_lines: List[str] = []
        for k in sorted(self.Er_by_cluster.keys()):
            paradigms_in_k = [p for p, c in self.cluster_assignment.items() if c == k]
            cluster_lines.append(
                f"\n\t  cluster {k}: Er={self.Er_by_cluster[k]*1e3:.1f} mV, "
                f"Rin={self.Rin_by_cluster[k]*1e-9:.2f} GOhm, "
                f"gl={1.0/self.Rin_by_cluster[k]*1e9:.1f} nS, "
                f"paradigms={paradigms_in_k}"
            )
        print(
            f"Parameter Estimates ({source_label}):"
            f"\n\tdrift clusters: {len(self.Er_by_cluster)}"
            + "".join(cluster_lines) +
            f"\n\tEt (per clamp, mV): [{Et_str}]"
            f"\n\txbeta: {self.x_beta:.3f}"
            f"\n\topt objective: {self.opt_objective:.3e}"
        )

    def _collect_unique_Iinj_sorted(self) -> np.ndarray:
        """Union of Iinj values across all stimuli, ascending."""
        all_Iinj: List[float] = []
        for stim in self.stimuli.values():
            all_Iinj.extend(float(i) for i in stim.Iinj.flatten())
        return np.array(sorted(set(all_Iinj)), dtype=np.float64)

    def _collect_max_Vm_by_Iinj(self) -> Dict[float, float]:
        """Per-clamp max(Vm) across all stimuli that contain that clamp."""
        out: Dict[float, List[float]] = {}
        for stim in self.stimuli.values():
            for k, I in enumerate(stim.Iinj.flatten()):
                out.setdefault(float(I), []).append(float(stim.Vm[k].max()))
        return {I: float(np.max(v)) for I, v in out.items()}

    # ---- Parameter cache: serialize / deserialize the expensive results ----
    def to_cache_dict(self) -> Dict[str, Any]:
        """
        Serializable dict of all computed parameters worth caching. Excludes
        timeseries (cheap to recompute, expensive to store) and quantities that
        are trivially recomputed (max_Vm, dt). Float keys in the dicts are
        converted to strings (JSON requirement); load reverses this.
        """
        return {
            "Cm":                   float(self.Cm),
            "Ee":                   _nan_to_none(self.Ee),
            "Ei":                   _nan_to_none(self.Ei),
            "Vss_per_stimulus":     {p: {str(k): v for k, v in d.items()} for p, d in self.Vss_per_stimulus.items()},
            "Er_Rin_per_stimulus":  {p: list(t) for p, t in self.Er_Rin_per_stimulus.items()},
            "cluster_assignment":   dict(self.cluster_assignment),
            "Er_by_cluster":        {str(k): float(v) for k, v in self.Er_by_cluster.items()},
            "Rin_by_cluster":       {str(k): float(v) for k, v in self.Rin_by_cluster.items()},
            "Vss":                  {p: {str(k): v for k, v in d.items()} for p, d in self.Vss.items()},
            "Et":                   {str(k): float(v) for k, v in self.Et.items()},
            "x_beta":               float(self.x_beta),
            "opt_objective":        float(self.opt_objective),
        }

    def _load_cached_params(self, d: Dict[str, Any]) -> None:
        """Inverse of to_cache_dict. Populates self with cached values."""
        # Cm comes from the parameters dataframe; sanity-check consistency
        if abs(float(d["Cm"]) - float(self.Cm)) > 1e-15:
            print(f"Warning: cached Cm ({d['Cm']}) differs from parameters Cm ({self.Cm}); using parameters value.")

        self.Ee = _none_to_nan(d.get("Ee"))
        self.Ei = _none_to_nan(d.get("Ei"))
        self.Vss_per_stimulus = {
            p: {float(k): float(v) for k, v in inner.items()}
            for p, inner in d["Vss_per_stimulus"].items()
        }
        self.Er_Rin_per_stimulus = {
            p: (float(t[0]), float(t[1])) for p, t in d["Er_Rin_per_stimulus"].items()
        }
        self.cluster_assignment = {p: int(c) for p, c in d["cluster_assignment"].items()}
        self.Er_by_cluster = {int(k): float(v) for k, v in d["Er_by_cluster"].items()}
        self.Rin_by_cluster = {int(k): float(v) for k, v in d["Rin_by_cluster"].items()}
        self.Vss = {
            p: {float(k): float(v) for k, v in inner.items()}
            for p, inner in d["Vss"].items()
        }
        self.Et = {float(k): float(v) for k, v in d["Et"].items()}
        self.x_beta = float(d["x_beta"])
        self.opt_objective = float(d["opt_objective"])

    def beta(self, paradigm: str, Iinj: np.ndarray) -> np.ndarray:
        """
        Per-clamp beta for a given paradigm, using its cluster's Vss and the
        global per-clamp Et. Returns shape [Nclamps, 1] aligned with `Iinj`.
        """
        Iinj_flat = Iinj.flatten()
        Et_arr: np.ndarray = np.array([[self.Et[float(i)]] for i in Iinj_flat])
        Vss_arr: np.ndarray = np.array([[self.Vss[paradigm][float(i)]] for i in Iinj_flat])
        return self.x_beta * self.gl(paradigm) / (Et_arr - Vss_arr)

    def gl(self, paradigm: str) -> float:
        """Leak conductance for a paradigm = 1 / Rin of its cluster."""
        cluster = self.cluster_assignment[paradigm]
        return 1.0 / self.Rin_by_cluster[cluster]

    def gact(self) -> Dict[str, Dict[float, float]]:
        """
        Per-paradigm, per-clamp gact: max slope of Iact at Vm = max_Vm_k.
        gact[paradigm][Iinj] = x_beta * gl(paradigm) * (max_Vm[Iinj] - Vss[paradigm][Iinj])
                                                     / (Et[Iinj]     - Vss[paradigm][Iinj]).
        Diagnostic only.
        """
        out: Dict[str, Dict[float, float]] = {}
        for paradigm, stim in self.stimuli.items():
            out[paradigm] = {}
            gl_p = self.gl(paradigm)
            for Iinj in np.squeeze(stim.Iinj, axis=-1):
                I_f = float(Iinj)
                Et_k = self.Et[I_f]
                Vss_k = self.Vss[paradigm][I_f]
                mV_k = self._max_Vm_by_Iinj[I_f]
                denom = Et_k - Vss_k
                if denom <= 0:
                    out[paradigm][I_f] = float("nan")
                else:
                    out[paradigm][I_f] = gl_p * (mV_k - Vss_k) / denom * self.x_beta
        return out

    @property
    def max_Vm(self) -> float:
        """Global max Vm across all clamps and stimuli (kept for backward compat)."""
        return float(max(self._max_Vm_by_Iinj.values()))

    # ---- Resting-state estimation with drift clustering ----
    def _estimate_resting_state(self) -> None:
        """
        Estimate (Er, Rin) and Vss for each stimulus, allowing for drift across
        stimuli (e.g. drug wash-in, dialysis). Pipeline:

            1. Per stimulus: histogram-mode Vss for each (paradigm, Iinj).
               -> self.Vss_per_stimulus[paradigm][Iinj]

            2. Per stimulus: 2-parameter OLS regression of (Iinj, Vss) within
               that stimulus.  Stimuli with <2 clamps fall back to NaN here.
               -> self.Er_Rin_per_stimulus[paradigm]

            3. Cluster paradigms by (Er, Rin) using a Gaussian mixture, with K
               selected by BIC over K in {1, ..., min(N-1, 4)}. Z-score the two
               axes before fitting so they're commensurable.
               -> self.cluster_assignment[paradigm]

            4. Per cluster: re-fit (Er, Rin) by pooling raw Vm traces across all
               stimuli in that cluster and re-running histogram-mode + OLS.
               -> self.Er_by_cluster, self.Rin_by_cluster

            5. Per paradigm: store the cluster's pooled Vss values, sampled at
               that paradigm's Iinj levels.
               -> self.Vss[paradigm][Iinj]
        """
        self._estimate_Vss_per_stimulus()
        self._fit_Er_Rin_per_stimulus()
        self._cluster_drift_states()
        self._refit_clusters_pooled()

    def _estimate_Vss_per_stimulus(self) -> None:
        """
        Weighted-mode Vss for each (paradigm, Iinj), using the filtered Vm trace
        and weighting samples by exp(-|Im_filtered| / Vss_Im_scale) so quiet
        baseline samples dominate the mode.
        """
        for paradigm, stim in self.stimuli.items():
            self.Vss_per_stimulus[paradigm] = {}
            for k, Iinj in enumerate(np.squeeze(stim.Iinj, axis=-1)):
                Vm_clamp = stim.Vm_filtered[k]
                Im_clamp = stim.Im_filtered[k]
                weights  = np.exp(-np.abs(Im_clamp) / self.Er_Rin_cfg.Vss_Im_scale)
                Vss_hat = self._weighted_smoothed_mode(Vm_clamp, weights)
                self.Vss_per_stimulus[paradigm][float(Iinj)] = Vss_hat

    def _weighted_smoothed_mode(self, samples: np.ndarray, weights: np.ndarray) -> float:
        """
        Volts-based weighted-histogram mode.
          - Bin edges spaced by Vss_bin_width (volts) over [min(samples), max(samples)].
          - Counts accumulate `weights` instead of unit increments.
          - Smoothed by a Gaussian of sigma Vss_smooth_sigma (volts), converted
            to bin units internally.
          - Returns the bin-center of the maximum smoothed count.
        """
        bin_w:    float = self.Er_Rin_cfg.Vss_bin_width
        sigma_v:  float = self.Er_Rin_cfg.Vss_smooth_sigma

        smin: float = float(samples.min())
        smax: float = float(samples.max())
        # Pad slightly to avoid edge-of-range artifacts.
        pad: float = 2.0 * sigma_v
        edges: np.ndarray = np.arange(smin - pad, smax + pad + bin_w, bin_w)
        if edges.size < 2:
            # Fall back to simple weighted mean if range is degenerate.
            return float(np.average(samples, weights=weights)) if weights.sum() > 0 else float(samples.mean())

        counts, _ = np.histogram(samples, bins=edges, weights=weights)
        sigma_bins: float = sigma_v / bin_w
        counts_s: np.ndarray = gaussian_filter1d(counts.astype(float), sigma_bins)
        centers: np.ndarray = 0.5 * (edges[:-1] + edges[1:])
        return float(centers[np.argmax(counts_s)])

    def _fit_Er_Rin_per_stimulus(self) -> None:
        """
        OLS fit of Vss = Er + Rin * Iinj within each stimulus.
        Stimuli with <2 distinct Iinj clamps cannot fit and get NaN; they will
        be assigned to the dominant cluster downstream by majority vote.
        """
        for paradigm, vss_dict in self.Vss_per_stimulus.items():
            Iinjs = np.array(list(vss_dict.keys()), dtype=np.float64)
            Vsss  = np.array(list(vss_dict.values()), dtype=np.float64)
            if Iinjs.size < 2 or np.std(Iinjs) == 0.0:
                self.Er_Rin_per_stimulus[paradigm] = (float("nan"), float("nan"))
                continue
            Er_hat, Rin_hat = self._ols_Er_Rin(Iinjs, Vsss)
            self.Er_Rin_per_stimulus[paradigm] = (Er_hat, Rin_hat)

    @staticmethod
    def _ols_Er_Rin(Iinjs: np.ndarray, Vsss: np.ndarray) -> Tuple[float, float]:
        """OLS regression Vss = Er + Rin * Iinj. Returns (Er, Rin)."""
        Imean = float(Iinjs.mean())
        Vmean = float(Vsss.mean())
        Ic = Iinjs - Imean
        Vc = Vsss - Vmean
        denom = float(np.sum(Ic * Ic))
        if denom <= 0:
            return float("nan"), float("nan")
        Rin_hat = float(np.sum(Ic * Vc) / denom)
        Er_hat = float(Vmean - Rin_hat * Imean)
        return Er_hat, Rin_hat

    def _cluster_drift_states(self) -> None:
        """
        Cluster paradigms by their per-stimulus (Er, Rin) using a Gaussian
        mixture; pick K by BIC over K in {1, ..., min(N_fittable - 1, 4)}.
        Stimuli that couldn't be fit (NaN) are assigned to the largest cluster.
        Sets self.cluster_assignment.
        """
        # Separate fittable from non-fittable stimuli
        fittable: List[str] = []
        feats: List[Tuple[float, float]] = []
        for paradigm, (Er, Rin) in self.Er_Rin_per_stimulus.items():
            if np.isfinite(Er) and np.isfinite(Rin):
                fittable.append(paradigm)
                feats.append((Er, Rin))

        if len(fittable) == 0:
            # No paradigm could be fit; everyone goes in cluster 0
            for p in self.stimuli:
                self.cluster_assignment[p] = 0
            return

        if len(fittable) == 1:
            # Trivially one cluster
            for p in self.stimuli:
                self.cluster_assignment[p] = 0
            return

        X = np.array(feats, dtype=np.float64)               # shape [N_fittable, 2]

        # Normalize each axis by a physical scale rather than sample SD. This is
        # critical with small N: z-scoring inflates within-noise variation to
        # unit scale, so even truly-identical stimuli look like they span unit
        # variance, and GMM K=2 fits this spurious structure with high log-
        # likelihood gain. Normalizing by a physical noise scale keeps
        # noise-only data at ~1 unit and real drift at >> 1 unit, so BIC
        # naturally prefers K=1 in the noise case.
        scales = np.array([self.Er_Rin_cfg.cluster_scale_Er,
                           self.Er_Rin_cfg.cluster_scale_Rin], dtype=np.float64)
        Xn = (X - X.mean(axis=0, keepdims=True)) / scales

        # Penalized BIC: BIC_alpha = -2 log L + alpha * k * log N. alpha = 1
        # reproduces standard BIC. alpha > 1 makes adding clusters harder.
        K_max = len(fittable) - 1
        alpha: float = self.Er_Rin_cfg.cluster_penalty_alpha
        N_fit: int = len(fittable)
        log_N: float = float(np.log(N_fit))
        best_K, best_bic, best_labels = 1, np.inf, np.zeros(N_fit, dtype=int)

        for K in range(1, K_max + 1):
            # Spherical covariance: each cluster has a single shared variance.
            # reg_covar floors cluster variance at cluster_noise_floor (in
            # normalized units). Without a meaningful floor, GMM can drive
            # variance to ~0 and log-likelihood to +inf, defeating BIC. The floor
            # should correspond to the typical per-stimulus noise on (Er, Rin)
            # estimates, expressed as a fraction of cluster_scale_*. E.g. with
            # cluster_scale_Er = 5 mV and typical Er-estimate noise of ~0.5 mV,
            # cluster_noise_floor = (0.5/5)^2 = 0.01.
            gm = GaussianMixture(
                n_components=K,
                covariance_type="spherical",
                n_init=5,
                random_state=self.Et_opt_cfg.rng_seed,
                reg_covar=self.Er_Rin_cfg.cluster_noise_floor,
            )
            try:
                gm.fit(Xn)
            except Exception:
                continue
            # gm.score(X) is mean log-likelihood per sample; multiply by N for total.
            log_L: float = float(gm.score(Xn) * N_fit)
            k_params: int = int(gm._n_parameters())
            bic: float = -2.0 * log_L + alpha * k_params * log_N
            if bic < best_bic:
                best_bic = bic
                best_K = K
                best_labels = gm.predict(Xn).astype(int)

        # Re-label so cluster ids are dense 0..best_K-1 in order of first appearance.
        # (sklearn already does this in practice but be explicit.)
        remap: Dict[int, int] = {}
        next_id = 0
        clean_labels: List[int] = []
        for lbl in best_labels:
            if int(lbl) not in remap:
                remap[int(lbl)] = next_id
                next_id += 1
            clean_labels.append(remap[int(lbl)])

        for paradigm, lbl in zip(fittable, clean_labels):
            self.cluster_assignment[paradigm] = int(lbl)

        # Non-fittable stimuli -> largest cluster
        if len(fittable) < len(self.stimuli):
            counts: Dict[int, int] = {}
            for lbl in clean_labels:
                counts[lbl] = counts.get(lbl, 0) + 1
            majority = max(counts, key=lambda k: counts[k])
            for p in self.stimuli:
                if p not in self.cluster_assignment:
                    self.cluster_assignment[p] = majority

    def _refit_clusters_pooled(self) -> None:
        """
        For each cluster, pool filtered Vm traces (with weights from filtered
        Im) across all member paradigms grouped by Iinj, and run weighted-mode +
        OLS on the pooled data. This gives the cluster's authoritative (Er, Rin)
        and Vss(Iinj). Each paradigm then inherits its cluster's Vss values at
        its own Iinj levels.
        """
        cluster_ids = sorted(set(self.cluster_assignment.values()))
        # Per cluster: Iinj -> pooled Vss_hat
        Vss_by_cluster: Dict[int, Dict[float, float]] = {}

        for k in cluster_ids:
            paradigms_in_k = [p for p, c in self.cluster_assignment.items() if c == k]

            # Pool filtered Vm samples + filtered Im (for weights) by Iinj across
            # this cluster's paradigms.
            Vm_by_Iinj: Dict[float, List[np.ndarray]] = {}
            Im_by_Iinj: Dict[float, List[np.ndarray]] = {}
            for p in paradigms_in_k:
                stim = self.stimuli[p]
                for clamp_idx, Iinj in enumerate(np.squeeze(stim.Iinj, axis=-1)):
                    Vm_by_Iinj.setdefault(float(Iinj), []).append(stim.Vm_filtered[clamp_idx])
                    Im_by_Iinj.setdefault(float(Iinj), []).append(stim.Im_filtered[clamp_idx])

            Vss_by_cluster[k] = {}
            for Iinj in Vm_by_Iinj:
                pooled_Vm = np.concatenate(Vm_by_Iinj[Iinj])
                pooled_Im = np.concatenate(Im_by_Iinj[Iinj])
                weights = np.exp(-np.abs(pooled_Im) / self.Er_Rin_cfg.Vss_Im_scale)
                Vss_by_cluster[k][Iinj] = self._weighted_smoothed_mode(pooled_Vm, weights)

            # OLS on the cluster-pooled (Iinj, Vss) pairs
            Iinjs = np.array(list(Vss_by_cluster[k].keys()), dtype=np.float64)
            Vsss  = np.array(list(Vss_by_cluster[k].values()), dtype=np.float64)
            if Iinjs.size >= 2 and np.std(Iinjs) > 0:
                Er_hat, Rin_hat = self._ols_Er_Rin(Iinjs, Vsss)
            else:
                # Single Iinj level in this cluster -- can't get Rin from it alone.
                # Fall back to the per-stimulus median for paradigms in this cluster.
                fallback_Rin = float(np.nanmedian([
                    self.Er_Rin_per_stimulus[p][1] for p in paradigms_in_k
                    if np.isfinite(self.Er_Rin_per_stimulus[p][1])
                ]))
                fallback_Er = float(np.nanmedian([
                    self.Er_Rin_per_stimulus[p][0] for p in paradigms_in_k
                    if np.isfinite(self.Er_Rin_per_stimulus[p][0])
                ]))
                Er_hat, Rin_hat = fallback_Er, fallback_Rin

            self.Er_by_cluster[k] = Er_hat
            self.Rin_by_cluster[k] = Rin_hat

        # Map each paradigm to its cluster's Vss dict, restricted to that
        # paradigm's actual Iinj levels.
        for paradigm, stim in self.stimuli.items():
            k = self.cluster_assignment[paradigm]
            self.Vss[paradigm] = {}
            for Iinj in np.squeeze(stim.Iinj, axis=-1):
                Iinj_f = float(Iinj)
                # Should always be in Vss_by_cluster[k] (we built it from the same stims)
                self.Vss[paradigm][Iinj_f] = Vss_by_cluster[k][Iinj_f]

    def estimate_Ee_Ei(self) -> Tuple[float, float]:
        """
        Estimate (Ee, Ei) by minimizing physical-bound violations.

        Background. The membrane equation gives target_Isyn(t) = Δgsyn(t)·V(t)
        - A(t), where Δgsyn(t) = Δge + Δgi and A(t) = Δge·Ee + Δgi·Ei. The
        cross-clamp regression yields Δgsyn(t) and A(t) directly, *without*
        needing Ee or Ei. For any (Ee, Ei) with Ee ≠ Ei:
            Δge(t) = (A(t) - Δgsyn(t)·Ei) / (Ee - Ei)
            Δgi(t) = (Δgsyn(t)·Ee - A(t)) / (Ee - Ei)

        Per-cluster physical bound on Δg. For cluster j with effective leak
        gl'_j = 1/Rin_j, the cell's pure leak gl satisfies gl ≤ gl'_j (because
        gl'_j absorbs tonic ge0_j + gi0_j ≥ 0). Using a global floor gl_min
        on the pure leak as a config-supplied prior, the tightest guaranteed
        lower bound on each component is
            Δge_j(t) ≥ -(gl'_j - gl_min)
            Δgi_j(t) ≥ -(gl'_j - gl_min)
        (Same number for both components: this is the loosest bound that
        respects ge0 + gi0 = gl'_j - gl, and ge0, gi0 ≥ 0.)

        Loss. Sum of squared violations across stimuli and timepoints:
            L(Ee, Ei) = Σ_j Σ_t [max(0, B_j - Δge_j(t))² + max(0, B_j - Δgi_j(t))²]
        where B_j = -(gl'_j - gl_min) is the per-cluster floor.

        Search. Grid-search over the box. The data fundamentally underdetermines
        (Ee, Ei) per cell (one cell, two unknowns, infinite-dimensional data
        but only one constraint per timepoint), so L typically has a flat
        minimum region. Tiebreak among the zero-loss (or near-min-loss) solutions
        by Euclidean distance to (Ee_prior_center, Ei_prior_center). The prior
        center is an honest fallback for the underdetermination -- not a
        statement that those values are correct, but a default when the data
        cannot pin them down.

        Returns (Ee, Ei) in volts.
        """
        cfg = self.Ee_Ei_cfg

        # --- Step 1: collect per-stimulus (Δgsyn, A, lower_bound) tuples ---
        # Each stimulus knows its paradigm and hence its cluster, so we can
        # apply the right per-cluster lower bound per timepoint.
        per_stim: List[Tuple[np.ndarray, np.ndarray, float]] = []
        for paradigm, stim in self.stimuli.items():
            dgsyn_arr = np.asarray(stim.timeseries["dgsyn filtered"], dtype=np.float64)
            A_arr     = np.asarray(stim.timeseries["A filtered"],     dtype=np.float64)
            valid = np.isfinite(dgsyn_arr) & np.isfinite(A_arr)
            if not np.any(valid):
                continue
            dgsyn_v = dgsyn_arr[valid]
            A_v     = A_arr[valid]
            # Per-cluster lower bound for THIS paradigm.
            gl_prime = self.gl(paradigm)               # = 1/Rin_by_cluster[cluster]
            B_j      = -(gl_prime - cfg.gl_min)        # negative number; the floor on Δg
            per_stim.append((dgsyn_v, A_v, B_j))

        if not per_stim:
            raise RuntimeError(
                "estimate_Ee_Ei: no stimuli have valid (Δgsyn, A) samples; "
                "cannot estimate (Ee, Ei)."
            )

        # --- Step 2: bound-violation loss as a function of (Ee, Ei) ---
        # We compute on a 2D grid for clarity; the box is small and grid eval
        # is fast. The loss is non-smooth (max(0, ·)²) but well-behaved for
        # gradient-free search.
        def violation_loss(Ee: float, Ei: float) -> float:
            if abs(Ee - Ei) < 1e-9:
                return float("inf")
            denom = Ee - Ei
            total = 0.0
            for dgsyn_v, A_v, B_j in per_stim:
                dge_v = (A_v - dgsyn_v * Ei) / denom
                dgi_v = (dgsyn_v * Ee - A_v) / denom
                # max(0, B_j - x) = max(0, x_below_floor) is the violation
                vio_e = np.maximum(0.0, B_j - dge_v)
                vio_i = np.maximum(0.0, B_j - dgi_v)
                total += float(np.sum(vio_e * vio_e) + np.sum(vio_i * vio_i))
            return total

        # --- Step 3: grid search over the box ---
        # Resolution: 0.5 mV per axis -> 41 x 91 = 3731 grid points for the
        # default box. Still fast even with ~1e6 timepoints across stimuli.
        grid_step: float = 0.5e-3
        Ee_grid: np.ndarray = np.arange(cfg.Ee_min, cfg.Ee_max + grid_step / 2, grid_step)
        Ei_grid: np.ndarray = np.arange(cfg.Ei_min, cfg.Ei_max + grid_step / 2, grid_step)

        loss_grid: np.ndarray = np.empty((Ee_grid.size, Ei_grid.size), dtype=np.float64)
        for i, Ee in enumerate(Ee_grid):
            for j, Ei in enumerate(Ei_grid):
                loss_grid[i, j] = violation_loss(float(Ee), float(Ei))

        # --- Step 4: tiebreak ---
        # Find all grid points within a small tolerance of the minimum loss;
        # among them, pick the one closest to the prior center.
        L_min: float = float(np.nanmin(loss_grid))
        # Tolerance handling: when L_min is exactly 0 (data is consistent with
        # bounds for many (Ee, Ei)) we still want to absorb floating-point
        # roundoff at the plateau boundary. Set absolute tolerance using a
        # data-driven scale: typical Δg² magnitude ~ var(Δgsyn). Floating-point
        # error in the violation computation scales with that, times machine
        # epsilon, times the number of timepoints summed.
        all_dgsyn = np.concatenate([d[0] for d in per_stim])
        scale = float(np.var(all_dgsyn)) if all_dgsyn.size > 0 else 1.0
        n_total = sum(d[0].size for d in per_stim)
        abs_tol_floor = scale * n_total * np.finfo(np.float64).eps * 100.0
        if L_min > 0:
            tol = max(L_min * 1e-9, abs_tol_floor)
        else:
            tol = abs_tol_floor
        feasible_mask: np.ndarray = loss_grid <= L_min + tol

        if not np.any(feasible_mask):
            # Should never happen since loss_grid contains L_min, but guard.
            return cfg.Ee_prior_center, cfg.Ei_prior_center

        # Distance from prior center for each (Ee_grid, Ei_grid) pair.
        Ee_mesh, Ei_mesh = np.meshgrid(Ee_grid, Ei_grid, indexing="ij")
        d2: np.ndarray = (Ee_mesh - cfg.Ee_prior_center) ** 2 \
                       + (Ei_mesh - cfg.Ei_prior_center) ** 2
        # Mask out infeasible points by setting their distance to +inf.
        d2_feasible: np.ndarray = np.where(feasible_mask, d2, np.inf)
        flat_idx: int = int(np.argmin(d2_feasible))
        i_best, j_best = np.unravel_index(flat_idx, d2_feasible.shape)
        return float(Ee_grid[i_best]), float(Ei_grid[j_best])

    def run_analysis(self, verbose: bool = True, complete: bool = True) -> None:
        if verbose:
            print("Estimating reversal potentials... ")

        for stimulus in self.stimuli.values():
            stimulus.calculate_target_Isyn()
            stimulus.estimate_Eeff_dgsyn()

        if complete:
            self.Ee, self.Ei = self.estimate_Ee_Ei()

            if verbose:
                print(f"Estimated Reversals: Ee = {self.Ee*1e3:.1f} mV, Ei = {self.Ei*1e3:.1f} mV")
                print("Calculating synaptic conductance deviations Δge, Δgi... ")

            for stimulus in self.stimuli.values():
                stimulus.estimate_dge_dgi()

    def evaluate_Cm_quality(self) -> Dict[str, float]:
        """
        Diagnostic for whether the supplied Cm value is consistent with the data.

        Theory: a Cm error contributes -dC * Cov_j(Vm, dVm/dt) / Var_j(Vm) to the
        inferred Δgsyn at each timepoint, where the cross-clamp covariance is taken
        across the stimulus's clamps at that timepoint. If Cm is biased, the
        magnitude of this term should correlate with the magnitude of negative-Δgsyn
        excursions: bigger |cross-clamp Vm/dVm covariance| -> bigger negative-Δgsyn
        events. If Cm is correct, the two should be uncorrelated (any remaining
        negatives come from other model error).

        We compute, per stimulus, the Spearman rho between
            x(t) = | Cov_j(Vm[j,t], dVm/dt[j,t]) / Var_j(Vm[j,t]) |
                   (the per-timepoint Cm-error sensitivity, in units of S/F)
            y(t) = min(Δgsyn(t), 0)^2
                   (the squared negative-Δgsyn excursion at t)

        and report two flavors:
            "rho_all":      using all timepoints in the analysis window
            "rho_negatives": using only timepoints where Δgsyn(t) < 0

        Each is aggregated across stimuli via Fisher-z averaging. NaN if too few
        valid timepoints.

        Interpretation: rho_negatives is the more discriminating of the two.
        Synthetic ground-truth tests show |rho_negatives| ~ 0.04 when Cm is
        correct and ~0.8 when Cm is off by a factor of 2, with the rest of the
        model held fixed. rho_all tends to light up even with correct Cm because
        Vm-dVm structure correlates with Δgsyn-magnitude in benign ways too.
        Use rho_negatives as the primary signal; treat rho_all as supportive.
        Thresholds are empirical -- calibrate on a recording where you trust Cm.

        CAVEAT (TODO item 3): under the corrected Δg framing, negative Δgsyn has
        two sources: (a) Cm/active-current misspecification, which is what this
        diagnostic targets; (b) real disinhibition / withdrawal of tonic input,
        which is biophysically valid signal. The current diagnostic cannot
        distinguish them. Real disinhibition is slower and not |dV/dt|-correlated,
        so it would not light up rho_all but might wash out rho_negatives. Plan
        to add a timescale separation step (high-pass before correlation) so the
        diagnostic is specific to PSP-timescale model misspecification. Until
        then, treat the diagnostic as detecting "PSP-timescale misspecification"
        rather than Cm-specific error.

        NOTE: This diagnostic reports magnitude only. A signed direction-of-error
        version was attempted (predicting sign(Delta_C) from sign(s_signed) at
        negative-Δgsyn timepoints), but synthetic tests showed the sign of s_signed
        at the negative-Δgsyn timepoints is determined more by *which phase of the
        PSP transient the negatives land in* than by sign(Delta_C). Negatives tend
        to cluster at the same PSP phase regardless of Cm bias, which means the
        sign signal cancels out. Use the magnitude diagnostic to flag a problem;
        determine the direction empirically by perturbing Cm and re-running.

        Notes:
        - Uses Vm and Δgsyn _after_ filtering, consistent with how target_Isyn was
          built.
        - dVm/dt is computed via central differences on the filtered Vm.
        - Stimuli with fewer than 2 clamps are skipped (covariance undefined).
        """
        rhos_all: List[float] = []
        rhos_neg: List[float] = []
        n_total: int = 0
        n_neg: int = 0

        for stim in self.stimuli.values():
            if stim.Nclamps < 2:
                continue

            Vm_filt: np.ndarray = np.stack(stim.timeseries["Vm filtered"]).astype(np.float64).T  # [Nclamps, Nsamples]
            dgsyn: np.ndarray = stim.timeseries["dgsyn filtered"].to_numpy(np.float64)           # [Nsamples]

            dVm: np.ndarray = np.gradient(Vm_filt, self.dt, axis=-1)

            # Per-timepoint cross-clamp moments (axis=0 averages across clamps).
            Vm_mean: np.ndarray = Vm_filt.mean(axis=0, keepdims=True)
            dVm_mean: np.ndarray = dVm.mean(axis=0, keepdims=True)
            Vm_c: np.ndarray = Vm_filt - Vm_mean
            dVm_c: np.ndarray = dVm - dVm_mean

            cov_VdV: np.ndarray = np.mean(Vm_c * dVm_c, axis=0)   # [Nsamples]
            var_V:   np.ndarray = np.mean(Vm_c * Vm_c, axis=0)    # [Nsamples]

            # Cm-error sensitivity per timepoint (units: 1/time = S/F).
            # Skip samples where var_V is tiny (clamps degenerate to a single voltage).
            eps: float = self.numerics_cfg.box_volume_eps
            valid: np.ndarray = var_V > eps
            sensitivity: np.ndarray = np.where(valid, np.abs(cov_VdV / np.where(valid, var_V, 1.0)), np.nan)

            neg_sq: np.ndarray = np.minimum(dgsyn, 0.0) ** 2

            # rho_all: over all valid timepoints with finite sensitivity and Δgsyn
            mask_all: np.ndarray = np.isfinite(sensitivity) & np.isfinite(dgsyn)
            if mask_all.sum() >= 3:
                rho, _ = spearmanr(sensitivity[mask_all], neg_sq[mask_all])
                if np.isfinite(rho):
                    rhos_all.append(float(rho))
                    n_total += int(mask_all.sum())

            # rho_neg: restricted to timepoints where Δgsyn < 0
            mask_neg: np.ndarray = mask_all & (dgsyn < 0.0)
            if mask_neg.sum() >= 3:
                rho, _ = spearmanr(sensitivity[mask_neg], neg_sq[mask_neg])
                if np.isfinite(rho):
                    rhos_neg.append(float(rho))
                    n_neg += int(mask_neg.sum())

        def fisher_avg(rs: List[float]) -> float:
            if not rs:
                return float("nan")
            # Clip to avoid arctanh blowup at exact +/-1
            zs = np.arctanh(np.clip(rs, -0.999999, 0.999999))
            return float(np.tanh(np.mean(zs)))

        return {
            "rho_all":         fisher_avg(rhos_all),
            "rho_negatives":   fisher_avg(rhos_neg),
            "n_stimuli_all":   float(len(rhos_all)),
            "n_stimuli_neg":   float(len(rhos_neg)),
            "n_samples_all":   float(n_total),
            "n_samples_neg":   float(n_neg),
        }

    """ Et / x_beta optimization """
    def _estimate_Et_beta(self, n_workers: int) -> Tuple[Dict[float, float], float, float]:
        """
        Optimize over theta = (Et_base, delta_1, ..., delta_{N-1}, x_beta).

        Constraints (all enforced by construction — no rejection sampling):
          - delta_j >= 0                        (monotonic non-decreasing Et with Iinj)
          - Et_sorted[k] >= max_Vm[k] + eps_v   (per-clamp lower bound)
          - x_beta in [0, x_beta_max]

        Strategy: random search in the feasible box, refining around the best point.
        Cheap and parallel; easy to swap for L-BFGS-B / CMA-ES later.
        """
        eps_v: float = self.Et_opt_cfg.eps_v

        Iinj_sorted: np.ndarray = self._Iinj_sorted_for_opt
        N: int = Iinj_sorted.size

        # Per-clamp lower bound for Et_sorted[k]
        max_Vm_sorted: np.ndarray = np.array(
            [self._max_Vm_by_Iinj[float(I)] for I in Iinj_sorted], dtype=np.float64
        )
        Et_lb_per_clamp: np.ndarray = max_Vm_sorted + eps_v
        Et_ub_global: float = 0.0

        if Et_ub_global <= Et_lb_per_clamp[0]:
            # Degenerate (cell sat above 0 mV?): bail out gracefully
            Et_dict: Dict[float, float] = {float(I): float(Et_lb_per_clamp[0]) for I in Iinj_sorted}
            return Et_dict, 0.0, float("nan")

        x_beta_max: float = self.Et_opt_cfg.x_beta_max

        n_epochs: int = self.Et_opt_cfg.n_epochs
        n_samples: int = self.Et_opt_cfg.n_samples_per_epoch
        shrink: float = self.Et_opt_cfg.shrink
        rng = np.random.default_rng(self.Et_opt_cfg.rng_seed)

        box_eps: float = self.numerics_cfg.box_volume_eps
        box_floor: float = self.numerics_cfg.box_volume_floor

        # We parameterize the search box directly in terms of the cumulative values
        # Et_sorted[k] (k = 0, ..., N-1) and x_beta. This makes the per-clamp lower
        # bound Et_sorted[k] >= Et_lb_per_clamp[k] a simple per-coordinate constraint
        # and lets us sample feasibly with 100% yield (no rejection): independent
        # uniforms per coordinate, then sort the Et_sorted columns to enforce
        # monotonicity, then clip up to the per-clamp lower bound (which itself is
        # non-decreasing in k, so clipping preserves monotonicity).
        #
        # box_lo[:N] / box_hi[:N] are the search bounds on Et_sorted[0..N-1];
        # box_lo[-1] / box_hi[-1] are bounds on x_beta.
        box_lo: np.ndarray = np.concatenate([Et_lb_per_clamp.copy(),    [0.0]])
        box_hi: np.ndarray = np.concatenate([np.full(N, Et_ub_global),  [x_beta_max]])

        # Always include the "no active current" baseline (x_beta = 0, Et at lower bound)
        theta_zero: np.ndarray = self._cum_to_theta(Et_lb_per_clamp, 0.0)
        best_theta: np.ndarray = theta_zero.copy()
        best_val: float = np.inf

        def sample_feasible(n: int) -> np.ndarray:
            """Draw n feasible thetas with 100% yield."""
            raw_Et: np.ndarray = rng.uniform(box_lo[:N], box_hi[:N], size=(n, N))
            sorted_Et: np.ndarray = np.sort(raw_Et, axis=1)
            sorted_Et = np.maximum(sorted_Et, Et_lb_per_clamp[None, :])
            xbeta_col: np.ndarray = rng.uniform(box_lo[-1], box_hi[-1], size=(n, 1))
            cum: np.ndarray = np.concatenate([sorted_Et, xbeta_col], axis=1)
            return np.stack([self._cum_to_theta(row[:N], row[N]) for row in cum], axis=0)

        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            for epoch in range(n_epochs):
                samples = sample_feasible(n_samples)

                # Always include current best + zero baseline
                always_include = np.stack([best_theta, theta_zero], axis=0)
                candidates = np.concatenate([always_include, samples], axis=0)

                # Pretty-print the current box bounds for Et_sorted[0] and x_beta
                print(
                    f"epoch {epoch}: Et_sorted[0] in [{box_lo[0]*1e3:.1f}, {box_hi[0]*1e3:.1f}] mV, "
                    f"Et_sorted[-1] in [{box_lo[N-1]*1e3:.1f}, {box_hi[N-1]*1e3:.1f}] mV, "
                    f"x_beta in [{box_lo[-1]:.2f}, {box_hi[-1]:.2f}] "
                    f"({len(candidates)} evaluations on {n_workers} workers)"
                )

                tasks = [(theta, deepcopy(self)) for theta in candidates]
                results = list(pool.map(_evaluate_active_params, tasks))

                for theta, val in results:
                    if val < best_val:
                        best_val = val
                        best_theta = theta

                Et_dict_now, xb_now = _theta_to_Et_dict(best_theta, Iinj_sorted)
                Et_str = ", ".join(f"{Et_dict_now[float(I)]*1e3:.1f}" for I in Iinj_sorted)
                print(f"  best so far: Et=[{Et_str}] mV, x_beta={xb_now:.3f}, obj={best_val:.3e}")

                # Shrink box around best, in the cumulative (Et_sorted, x_beta) parameterization.
                # First convert best_theta -> cumulative form, then shrink each coordinate
                # symmetrically around its best value, clipped to the original feasible bounds.
                best_Et_sorted: np.ndarray = self._theta_to_Et_sorted(best_theta)
                best_cum: np.ndarray = np.concatenate([best_Et_sorted, [best_theta[-1]]])

                half: np.ndarray = 0.5 * shrink * (box_hi - box_lo)
                new_lo: np.ndarray = np.maximum(box_lo, best_cum - half)
                new_hi: np.ndarray = np.minimum(box_hi, best_cum + half)

                # Floor each Et_sorted lower bound by its per-clamp physical floor;
                # x_beta lower bound by 0.
                new_lo[:N] = np.maximum(new_lo[:N], Et_lb_per_clamp)
                new_lo[-1] = max(new_lo[-1], 0.0)

                # Make sure box still has volume on every axis
                degenerate = new_hi <= new_lo + box_eps
                if np.any(degenerate):
                    new_hi = np.where(degenerate, new_lo + box_floor, new_hi)
                box_lo, box_hi = new_lo, new_hi

        Et_dict, x_beta = _theta_to_Et_dict(best_theta, Iinj_sorted)
        return Et_dict, x_beta, float(best_val)

    @staticmethod
    def _cum_to_theta(Et_sorted: np.ndarray, x_beta: float) -> np.ndarray:
        """Convert (Et_sorted[0..N-1], x_beta) -> theta = (Et_base, deltas, x_beta)."""
        Et_base: float = float(Et_sorted[0])
        deltas: np.ndarray = np.diff(Et_sorted)
        return np.concatenate([[Et_base], deltas, [x_beta]])

    @staticmethod
    def _theta_to_Et_sorted(theta: np.ndarray) -> np.ndarray:
        """Convert theta -> Et_sorted[0..N-1] (drops x_beta)."""
        N: int = theta.size - 1
        Et_base: float = float(theta[0])
        deltas: np.ndarray = np.maximum(theta[1:N], 0.0)
        return Et_base + np.concatenate([[0.0], np.cumsum(deltas)])



class WholeCellStimulus:
    def __init__(self, recording: WholeCellRecording, data: pd.DataFrame, paradigm: str) -> None:
        self.recording: WholeCellRecording = recording
        self.paradigm: str = paradigm

        self.times: np.ndarray = data["times"].to_numpy(dtype=np.float64)

        Iinj_colnames: List[str] = list(data.columns)
        
        for aux_key in ["times", "stimulus", "representative"]:
            if aux_key in Iinj_colnames:
                Iinj_colnames.remove(aux_key)

        # Vm ingestion. Two transformations:
        #   1. Scale from millivolts (spreadsheet convention) to volts.
        #   2. Subtract the liquid junction potential.
        #
        # LJP convention (Marino et al. 2014, LJPcalc): V_LJP is the bath
        # potential relative to the pipette. The amplifier records V_measured =
        # V_true + V_LJP because the amplifier was zeroed in bath. To recover
        # the true membrane potential, subtract V_LJP from every reading.
        # This is the ONLY place this correction is applied in the pipeline;
        # all downstream code sees LJP-corrected voltages.
        V_LJP: float = recording.recording_cfg.liquid_junction_potential_volts
        self.Vm: np.ndarray = data[Iinj_colnames].to_numpy(dtype=np.float64).T * 1e-3 - V_LJP   # units: Volts, shape: [Nclamps, Nsamples]
        print(self.Vm[0,0], V_LJP, data[Iinj_colnames].to_numpy(dtype=np.float64)[0,0], "!!!!!")
        self.Iinj: np.ndarray = np.array(Iinj_colnames, dtype=np.float64)[:, np.newaxis]    # units: Amperes, shape: [Nclamps, 1]
                
        self.Nclamps: int = self.Iinj.size
        self.Nsamples: int = self.times.size

        # Precomputed traces. Filled by _precompute_filtered_traces (called once
        # before resting-state estimation), then Il by _precompute_Il (after
        # resting-state estimation has filled self.recording.Vss). Read by both
        # the Vss-estimation pipeline and calculate_target_Isyn.
        self.Vm_filtered: np.ndarray = np.empty((0,))     # [Nclamps, Nsamples]
        self.Im:          np.ndarray = np.empty((0,))     # [Nclamps, Nsamples]
        self.Im_filtered: np.ndarray = np.empty((0,))     # [Nclamps, Nsamples]
        self.Il:          np.ndarray = np.empty((0,))     # [Nclamps, Nsamples]
        self._Vss_clamp_arr: np.ndarray = np.empty((0,))  # [Nclamps, 1]

        self.timeseries: pd.DataFrame = pd.DataFrame({"times": self.times})

    def precompute_filtered_traces(self) -> None:
        """Compute Vm_filtered, Im, Im_filtered. Stable across the optimization."""
        Cm: float = self.recording.Cm
        dt: float = self.recording.dt
        fs: float = 1 / dt

        self.Vm_filtered = self.recording.filters["Vm"].propagate(self.Vm, fs)
        self.Im          = Cm * np.gradient(self.Vm, dt, axis=-1)
        self.Im_filtered = self.recording.filters["Im"].propagate(self.Im, fs)

    def precompute_Il(self) -> None:
        """
        Compute Il (and the per-clamp Vss array used in target-Isyn) using the
        cluster-pooled Vss for this paradigm. Must be called after the resting-
        state pipeline has populated self.recording.Vss.
        """
        Iinj_flat = self.Iinj.flatten()
        self._Vss_clamp_arr = np.array([[self.recording.Vss[self.paradigm][float(i)]] for i in Iinj_flat])
        gl: float = self.recording.gl(self.paradigm)
        self.Il = gl * (self._Vss_clamp_arr - self.Vm_filtered)

    def calculate_target_Isyn(self) -> None:
        """
        Compute Iact and target_Isyn from the current (Et, x_beta) opt params.
        Reads precomputed Vm_filtered, Im, Il, _Vss_clamp_arr from self
        (filled once at init).
        """
        Iinj: np.ndarray = self.Iinj
        Vss = self._Vss_clamp_arr
        Et: np.ndarray = np.array([[self.recording.Et[float(i)]] for i in Iinj.flatten()])

        beta: np.ndarray = self.recording.beta(self.paradigm, Iinj)

        Vm_filtered = self.Vm_filtered
        Im = self.Im
        Im_filtered = self.Im_filtered
        Il = self.Il

        # Per-clamp Et and beta broadcast over samples; gating still per-clamp via Vss
        Iact: np.ndarray = np.where(Vm_filtered > Vss, beta * (Vm_filtered - Vss) * (Et - Vm_filtered), 0)

        target_Isyn = -Im + Il - Iact

        self.timeseries["Vm filtered"] = Vm_filtered.T.tolist()
        self.timeseries["Im"] = Im.T.tolist()
        self.timeseries["Im filtered"] = Im_filtered.T.tolist()
        self.timeseries["Il"] = Il.T.tolist()
        self.timeseries["Iact"] = Iact.T.tolist()
        self.timeseries["target Isyn"] = target_Isyn.T.tolist()

    def estimate_Eeff_dgsyn(self) -> None:
        """
        Estimate effective synaptic reversal Eeff and total synaptic
        deviation-from-baseline conductance Δgsyn from cross-clamp regression
        of target_Isyn against Vm.

        Δgsyn = Δge + Δgi is signed (estimator places no sign constraint).
        It represents the *change* in total synaptic conductance from
        prestimulus baseline; negative values are biophysically valid.
        """
        Vm_filtered: np.ndarray = np.stack(self.timeseries["Vm filtered"]).astype(np.float64).T  # type: ignore
        target_Isyn: np.ndarray = np.stack(self.timeseries["target Isyn"]).astype(np.float64).T  # type: ignore

        Vm_filtered_mean: np.ndarray = Vm_filtered.mean(axis=0)
        target_Isyn_mean: np.ndarray = target_Isyn.mean(axis=0)

        centered_integral_Vm: np.ndarray = Vm_filtered - Vm_filtered_mean[None, :]
        centered_target_Isyn: np.ndarray = target_Isyn - target_Isyn_mean[None, :]

        denom = np.sum(centered_integral_Vm * centered_integral_Vm, axis=0)
        numer = np.sum(centered_integral_Vm * centered_target_Isyn, axis=0)

        eps = np.finfo(np.float64).tiny
        b = np.where(np.abs(denom) > eps, numer / denom, np.nan)
        a = target_Isyn_mean - b * Vm_filtered_mean

        # dgsyn_zero_threshold is a numerical floor for division (Eeff = -a/Δgsyn);
        # it gates on |Δgsyn|, NOT on sign(Δgsyn). Negative Δgsyn with
        # |Δgsyn| > threshold is valid signal and passes through unchanged.
        dgsyn = b
        dgsyn_safe = np.where(np.abs(dgsyn) > self.recording.numerics_cfg.dgsyn_zero_threshold, dgsyn, np.nan)

        Eeff = -a / dgsyn_safe

        # A(t) = Δge·Ee + Δgi·Ei = -a, used by the analytical (Ee, Ei) estimator
        # in WholeCellRecording.estimate_Ee_Ei. Stored alongside Δgsyn so the
        # downstream regression can run without re-deriving from target_Isyn.
        A_data = -a

        I_hat = a[None, :] + b[None, :] * Vm_filtered
        resid = target_Isyn - I_hat
        sse = np.sum(resid * resid, axis=0)

        sst = np.sum((target_Isyn - target_Isyn_mean[None, :]) ** 2, axis=0)
        r2 = np.where(sst > 0, 1.0 - (sse / sst), np.nan)

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
            Eeff_filtered[nan_mask] = np.nan

        dgsyn_filtered = self.recording.filters["dgsyn"].propagate(dgsyn[None, :], fs)[0]
        # A(t) shares the same low-pass filter as Δgsyn since they enter the
        # downstream Ee/Ei estimator together.
        A_filtered = self.recording.filters["dgsyn"].propagate(A_data[None, :], fs)[0]

        self.timeseries["Eeff"] = Eeff
        self.timeseries["Eeff filtered"] = Eeff_filtered
        self.timeseries["dgsyn"] = dgsyn_safe
        self.timeseries["dgsyn filtered"] = dgsyn_filtered
        self.timeseries["A"] = A_data
        self.timeseries["A filtered"] = A_filtered
        self.timeseries["SSE least-squares Qsyn"] = sse
        self.timeseries["r2 least-squares Qsyn"] = r2

    def estimate_dge_dgi(self) -> None:
        """
        Solve for Δge(t) and Δgi(t): the *signed* deviations of excitatory and
        inhibitory conductances from their prestimulus baseline values
        (g_e0 and g_i0, which are absorbed into gleak per the model).

        Sign convention: Δge, Δgi ∈ [−g_e0, ∞) and [−g_i0, ∞) respectively.
        Negative values are biophysically valid (disinhibition / withdrawal of
        tonic input) and should NOT be assumed to be artifacts. See TODO.md
        items 2, 3 for principled handling of misspecification vs. real
        negative excursions.
        """
        Vm_filtered: np.ndarray = np.stack(self.timeseries["Vm filtered"]).astype(np.float64).T  # type: ignore
        target_Isyn: np.ndarray = np.stack(self.timeseries["target Isyn"]).astype(np.float64).T  # type: ignore

        Ee: float = self.recording.Ee
        Ei: float = self.recording.Ei

        driving_force: np.ndarray = np.stack([Vm_filtered - Ee, Vm_filtered - Ei], axis=0)

        XtX: np.ndarray = np.einsum("ijn,kjn->nik", driving_force, driving_force)
        Xty: np.ndarray = np.einsum("ijn,jn->ni", driving_force, target_Isyn)

        try:
            conductances = np.linalg.solve(XtX, Xty[..., None])[..., 0]
        except np.linalg.LinAlgError:
            print("Cannot use least-squares, attempting pseudoinverse...")
            conductances = (np.linalg.pinv(XtX) @ Xty[..., None])[..., 0]

        dge: np.ndarray = conductances[:, 0]
        dgi: np.ndarray = conductances[:, 1]

        fs: float = 1 / self.recording.dt
        dge_filtered: np.ndarray = self.recording.filters["dg"].propagate(dge[None, :], fs)[0]
        dgi_filtered: np.ndarray = self.recording.filters["dg"].propagate(dgi[None, :], fs)[0]

        self.timeseries["dge"] = dge.tolist()
        self.timeseries["dgi"] = dgi.tolist()
        self.timeseries["dge filtered"] = dge_filtered.tolist()
        self.timeseries["dgi filtered"] = dgi_filtered.tolist()


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

            filtered_Vm: np.ndarray = np.stack(stimulus.timeseries["Vm filtered"]).astype(np.float64).T  # type: ignore
            Eeff: np.ndarray = stimulus.timeseries["Eeff filtered"].to_numpy(np.float64)
            Im: np.ndarray = np.stack(stimulus.timeseries["Im filtered"]).astype(np.float64).T  # type: ignore
            dge: np.ndarray = stimulus.timeseries["dge filtered"].to_numpy(np.float64)
            dgi: np.ndarray = stimulus.timeseries["dgi filtered"].to_numpy(np.float64)
            dgsyn: np.ndarray = stimulus.timeseries["dgsyn filtered"].to_numpy(np.float64)

            axs[0, idx].set_title(paradigm)

            curves = axs[0, idx].plot(times, filtered_Vm.T)
            colors = [line.get_color() for line in curves]

            Vss = [recording.Vss[paradigm][float(i)] for i in stimulus.Iinj.flatten()]
            for j in range(stimulus.Nclamps):
                axs[0, idx].plot([times[0], times[-1]], [Vss[j]] * 2, color="grey", ls="--")

            axs[0, idx].set_ylabel("Vm (V)")
            axs[0, idx].grid(True)

            axs[1, idx].grid(True)
            axs[1, idx].plot(times, Eeff, c="black")
            Er_paradigm = recording.Er_by_cluster[recording.cluster_assignment[paradigm]]
            axs[1, idx].axhline(Er_paradigm, linestyle="--", color="k", linewidth=1, label="Er" if idx == 0 else None)
            axs[1, idx].axhline(recording.Ee, linestyle="--", color="r", linewidth=1, label="Ee" if idx == 0 else None)
            axs[1, idx].axhline(recording.Ei, linestyle="--", color="b", linewidth=1, label="Ei" if idx == 0 else None)

            for j in range(stimulus.Nclamps):
                axs[2, idx].plot(times, Im[j, :], color=colors[j])
            axs[2, idx].plot([times[0], times[-1]], [0, 0], color="grey", ls="--")
            axs[2, idx].set_ylabel("Im (A)")
            axs[2, idx].grid(True)

            axs[3, idx].plot(times, dge, c="r", label="Δge")
            axs[3, idx].plot(times, dgi, c="b", label="Δgi")
            axs[3, idx].plot(times, dgsyn, c="k", label="Δgsyn", linestyle="--")
            axs[3, idx].plot(times, times * 0, "--k", linewidth=1)
            axs[3, idx].set_ylabel("ΔG (S)")
            axs[3, idx].grid(True)

            Iact_pred: np.ndarray = np.stack(stimulus.timeseries["Iact"]).astype(np.float64).T  # type: ignore
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
            "Vm":    LowPassFilter.from_cfg(self.cfg.filters.Vm),
            "Im":    LowPassFilter.from_cfg(self.cfg.filters.Im),
            "dg":    LowPassFilter.from_cfg(self.cfg.filters.dg),
            "Eeff":  LowPassFilter.from_cfg(self.cfg.filters.Eeff),
            "dgsyn": LowPassFilter.from_cfg(self.cfg.filters.dgsyn),
        }

        paths = self.cfg.paths.paths_to_spreadsheets
        cache_dir = self.cfg.paths.cache_dir
        cfg_hash = self.cfg.hash_for_cache()

        for path_to_spreadsheet in paths:
            print(f"Collecting data from {path_to_spreadsheet}")

            # Cache lookup
            cached_params: Optional[Dict[str, Any]] = None
            cache_path: Optional[Path] = None
            if cache_dir is not None:
                cache_path = cache_dir / f"{path_to_spreadsheet.stem}.json"
                cached_params = self._try_load_cache(
                    cache_path, path_to_spreadsheet, cfg_hash
                )

            rdr: XLReader = XLReader(path_to_spreadsheet)
            stimuli: Dict[str, pd.DataFrame] = {
                paradigm: rdr.get_paradigm_data(paradigm) for paradigm in rdr.get_paradigms()
            }
            recording: WholeCellRecording = WholeCellRecording(
                parameters=rdr.get_paradigm_parameters(rdr.get_paradigms()[0]),
                stimuli=stimuli,
                filters=filters,
                n_workers=self.cfg.compute.n_workers,
                recording_cfg=self.cfg.recording,
                Er_Rin_cfg=self.cfg.Er_Rin_estimation,
                Ee_Ei_cfg=self.cfg.Ee_Ei_estimation,
                Et_opt_cfg=self.cfg.Et_optimization,
                numerics_cfg=self.cfg.numerics,
                cached_params=cached_params,
            )
            del rdr

            recording.run_analysis()

            # Save cache after a successful run (only when we just computed it).
            if cache_dir is not None and cached_params is None and cache_path is not None:
                self._save_cache(cache_path, recording, path_to_spreadsheet, cfg_hash)

            if display:
                self.plot_timeseries(
                    recording,
                    self.cfg.paths.image_save_dir / path_to_spreadsheet.stem,
                    filetype=self.cfg.paths.image_save_type,
                    display=display,
                )

    def _try_load_cache(
        self,
        cache_path: Path,
        spreadsheet_path: Path,
        cfg_hash: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Returns the cached params dict if the cache is valid (file exists, cfg
        hash matches, spreadsheet mtime matches, and cache_invalidate is False).
        Otherwise returns None.
        """
        if self.cfg.paths.cache_invalidate:
            return None
        if not cache_path.is_file():
            return None
        try:
            with open(cache_path, "r") as f:
                payload = json.load(f)
        except Exception as e:
            print(f"  cache: failed to read {cache_path}: {e}; will recompute.")
            return None

        cached_hash = payload.get("cfg_hash")
        cached_mtime = payload.get("spreadsheet_mtime")
        actual_mtime = spreadsheet_path.stat().st_mtime

        if cached_hash != cfg_hash:
            print(f"  cache: cfg_hash mismatch (cached={cached_hash}, current={cfg_hash}); will recompute.")
            return None
        if cached_mtime is None or abs(float(cached_mtime) - actual_mtime) > 1e-6:
            print(f"  cache: spreadsheet mtime changed; will recompute.")
            return None

        print(f"  cache: HIT -> {cache_path}")
        return payload["params"]

    def _save_cache(
        self,
        cache_path: Path,
        recording: WholeCellRecording,
        spreadsheet_path: Path,
        cfg_hash: str,
    ) -> None:
        """Write recording params + cache metadata to disk."""
        payload = {
            "cfg_hash":          cfg_hash,
            "spreadsheet_mtime": spreadsheet_path.stat().st_mtime,
            "spreadsheet_name":  spreadsheet_path.name,
            "params":            recording.to_cache_dict(),
        }
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_path, "w") as f:
                json.dump(payload, f, indent=2)
            print(f"  cache: saved -> {cache_path}")
        except Exception as e:
            print(f"  cache: failed to write {cache_path}: {e}")