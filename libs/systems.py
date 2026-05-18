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

# Plotting & Graphics
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pyhelpers.store import save_fig

# Local
from libs.readers import (
    XLReader,
    AnalyzerCfg,
    FilterCfg,
    ErRinCfg,
    EeEiCfg,
    EtOptCfg,
    NumericsCfg,
    ManualClusterDataCfg,
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
        """
        Zero-phase low-pass filter the signal with DC-replicate boundary padding.

        Standard sosfiltfilt with default ('odd') padding rings on the boundary
        of a steady-state signal: it reflects the signal with sign flip around
        the endpoint, creating a discontinuity in derivative that the filter
        rings on. The ring decays over ~filter_order / passband_corner seconds
        and can have peak deviation of ~1 mV for typical 10-order Butterworth
        filters at 30 Hz — large enough to corrupt the cross-clamp Δg
        regression downstream.

        Fix: pad the signal with a long stretch of constant DC (the mean of
        a window at each endpoint), filter, then crop the pad back off. With
        ~100 ms of pad, the filter's transient ring fully decays inside the
        artificial padding and the cropped output stays within ~50 μV of true
        DC at the boundary -- ~20× better than the default.

        The DC pad is anchored to the mean of the first/last ~100 samples
        (10 ms at typical 10 kHz fs), which is robust to single-sample
        outliers at the very edge.
        """
        sos = self._design(fs)
        # Pad length: ~100 ms is safe for any reasonable filter we'd use here.
        # If fs varies (it shouldn't, but be defensive), scale with fs.
        pad_n: int = int(round(0.100 * fs))
        # Use the mean of a small endpoint window for the pad value, to avoid
        # anchoring the pad on a single noisy sample. 10 ms window.
        anchor_n: int = min(int(round(0.010 * fs)), max(1, raw_signal.shape[-1] // 4))

        # Operate along the last axis (matches sosfiltfilt's default).
        # Build endpoint anchors:
        head_anchor = np.mean(raw_signal[..., :anchor_n], axis=-1, keepdims=True)
        tail_anchor = np.mean(raw_signal[..., -anchor_n:], axis=-1, keepdims=True)

        # Broadcast pads to the leading shape of raw_signal.
        head_shape = list(raw_signal.shape); head_shape[-1] = pad_n
        tail_shape = list(raw_signal.shape); tail_shape[-1] = pad_n
        head_pad = np.broadcast_to(head_anchor, head_shape).copy()
        tail_pad = np.broadcast_to(tail_anchor, tail_shape).copy()

        padded = np.concatenate([head_pad, raw_signal, tail_pad], axis=-1)
        # padtype='constant' is reasonable here; sosfiltfilt still adds a small
        # internal pad but with the DC value, which is harmless.
        filtered_padded = sosfiltfilt(sos, padded, padtype='constant')
        # Crop the explicit pads back off; output shape matches input shape.
        return filtered_padded[..., pad_n:-pad_n]

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


def _evaluate_active_params(args: Tuple[np.ndarray, "WholeCellRecording"]) -> Tuple[np.ndarray, float, float]:
    """
    Worker: takes (theta, recording_copy), mutates the copy, runs analysis,
    returns the SOS-negative-Δgsyn loss summed over all clamps and stimuli.

    Returns three values per evaluation:
        theta:          the candidate point (echoed back for the caller)
        loss_filtered:  SSE on negative-part of filtered Δgsyn -- this is the
                        loss the optimizer minimizes
        loss_raw:       SSE on negative-part of raw (unfiltered) Δgsyn -- a
                        diagnostic comparator; not used for selection

    The two are mathematically distinct quantities. Filtering smears values
    across time, so the filtered Δgsyn has different magnitudes (and possibly
    different signs at any given timepoint) than the raw Δgsyn. We track both
    so we can verify they're consistent in practice; if they diverge
    substantially, the optimizer might be tuning α/β against filter artifacts
    rather than the underlying physical violations. See TODO.md item 2 for
    discussion of the broader negative-Δg interpretation question.
    """
    theta, rec = args
    Iinj_sorted: np.ndarray = rec._Iinj_sorted_for_opt
    Et_by_Iinj, x_beta = _theta_to_Et_dict(theta, Iinj_sorted)

    rec.Et = Et_by_Iinj
    rec.x_beta = x_beta
    rec.run_analysis(verbose=False, complete=False)

    loss_filtered: float = 0.0
    loss_raw:      float = 0.0
    for stim in rec.stimuli.values():
        dgsyn_filt: np.ndarray = stim.timeseries["dgsyn filtered"].to_numpy()
        dgsyn_raw:  np.ndarray = stim.timeseries["dgsyn"].to_numpy()
        loss_filtered += float(np.sum(np.square(np.minimum(dgsyn_filt, 0.0))))
        loss_raw      += float(np.sum(np.square(np.minimum(dgsyn_raw,  0.0))))

    return theta, loss_filtered, loss_raw


class WholeCellRecording:
    def __init__(
        self,
        parameters: pd.DataFrame,
        stimuli: Dict[str, pd.DataFrame],
        filters: Dict[str, LowPassFilter],
        n_workers: int,
        Er_Rin_cfg: ErRinCfg,
        Ee_Ei_cfg: EeEiCfg,
        Et_opt_cfg: EtOptCfg,
        numerics_cfg: NumericsCfg,
        manual_cluster_data: Optional[ManualClusterDataCfg] = None,
        cached_params: Optional[Dict[str, Any]] = None,
    ) -> None:
        assert len(stimuli) > 0
        assert set(filters.keys()) == {"Vm", "Im", "dg", "Eeff", "dgsyn"}

        self.filters: Dict[str, LowPassFilter] = filters

        # Cfg objects -- stored so workers and stimuli can read them after deepcopy.
        self.Er_Rin_cfg: ErRinCfg = Er_Rin_cfg
        self.Ee_Ei_cfg: EeEiCfg = Ee_Ei_cfg
        self.Et_opt_cfg: EtOptCfg = Et_opt_cfg
        self.numerics_cfg: NumericsCfg = numerics_cfg

        # Optional manual cluster data (cluster_assignments.json + steady_states/
        # in the case folder). If present, the resting-state pipeline uses this
        # directly instead of the unsupervised weighted-mode + GMM/BIC. None
        # means "fall back to unsupervised mode". See readers.ManualClusterDataCfg.
        self.manual_cluster_data: Optional[ManualClusterDataCfg] = manual_cluster_data

        # Validate manual data paradigm names against actual stimulus sheets.
        if manual_cluster_data is not None:
            for p in manual_cluster_data.cluster_for_paradigm:
                assert p in stimuli, (
                    f"cluster_assignments.json references paradigm '{p}' which is "
                    f"not present in this case's xlsx. "
                    f"Available paradigms: {sorted(stimuli.keys())}"
                )

        # Tracks whether the resting-state result came from manual data (True)
        # or unsupervised GMM/BIC (False). Used for pretty-print labeling. In the
        # current design, all clusters from a single run share an origin (either
        # all manual or all fit), so this is a per-cluster flag rather than mixed.
        self.cluster_is_manual: Dict[int, bool] = {}

        self.Cm: float = parameters["Cm"][0]

        example_t: pd.Series = list(stimuli.values())[0]["times"]
        self.dt: float = example_t[1] - example_t[0]

        self.Ee: float = np.nan
        self.Ei: float = np.nan

        # Cm-quality diagnostic results, populated by evaluate_Cm_quality()
        # at the end of run_analysis. See docstring on that method for details.
        # rho_negatives is the primary diagnostic; rho_all is supportive.
        self.Cm_rho_all:       float = np.nan
        self.Cm_rho_negatives: float = np.nan

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
            # Step 2: resting-state pipeline.
            # If manual cluster data is present, use it directly. Otherwise
            # run the unsupervised weighted-mode + GMM/BIC + cluster refit.
            if self.manual_cluster_data is not None:
                self._estimate_resting_state_manual()
            else:
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
            origin = "manual" if self.cluster_is_manual.get(k, False) else "fit"
            cluster_lines.append(
                f"\n\t  cluster {k} ({origin}): Er={self.Er_by_cluster[k]*1e3:.1f} mV, "
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
            "Cm_rho_all":           _nan_to_none(self.Cm_rho_all),
            "Cm_rho_negatives":     _nan_to_none(self.Cm_rho_negatives),
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
        self.Cm_rho_all       = _none_to_nan(d.get("Cm_rho_all"))
        self.Cm_rho_negatives = _none_to_nan(d.get("Cm_rho_negatives"))
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

            1. Per stimulus: Vss = median of first Vss_duration_seconds of Vm
               for each (paradigm, Iinj). Assumes the trace is in equilibrium
               for at least that long before any stimulus arrives. STOPGAP --
               replaces the older weighted-mode-on-full-trace estimator, which
               was unreliable when traces had no clean steady-state period.
               -> self.Vss_per_stimulus[paradigm][Iinj]

            2. Per stimulus: 2-parameter OLS regression of (Iinj, Vss) within
               that stimulus.  Stimuli with <2 clamps fall back to NaN here.
               -> self.Er_Rin_per_stimulus[paradigm]

            3. Cluster paradigms by Vss-pooling residual: greedy agglomerative
               merging until any further merge would induce a |pooled - observed|
               residual exceeding max_residual_volts at some (paradigm, Iinj).
               -> self.cluster_assignment[paradigm]

            4. Per cluster: re-fit (Er, Rin) by pooling first-Vss_duration_seconds
               Vm samples across all paradigms in the cluster at each Iinj,
               taking the median, then OLS.
               -> self.Er_by_cluster, self.Rin_by_cluster, self.Vss[...]
        """
        self._estimate_Vss_per_stimulus()
        self._fit_Er_Rin_per_stimulus()
        self._cluster_drift_states()
        self._refit_clusters_pooled()

    def _estimate_resting_state_manual(self) -> None:
        """
        Populate Vss, cluster_assignment, Er_by_cluster, Rin_by_cluster from
        user-provided manual cluster data, skipping the unsupervised pipeline.

        Steps:
          1. Map each cluster_name to an integer cluster_id (dense, 0..K-1).
          2. cluster_assignment[paradigm] = id of its named cluster.
          3. Vss[paradigm][Iinj] = manual Vss for (cluster, Iinj).
          4. (Er, Rin) per cluster by OLS on the manual (Iinj, Vss) pairs.

        The manual cluster names are preserved in self.cluster_labels for
        printout. cluster_is_manual is set True for every cluster.
        """
        mc = self.manual_cluster_data
        assert mc is not None, "_estimate_resting_state_manual called without manual data"

        # Assign dense integer IDs in sorted-by-name order for stable identity.
        cluster_names = sorted(set(mc.cluster_for_paradigm.values()))
        name_to_id: Dict[str, int] = {n: i for i, n in enumerate(cluster_names)}
        self.cluster_labels: Dict[int, str] = {
            name_to_id[n]: mc.cluster_labels.get(n, n) for n in cluster_names
        }

        # Cluster assignment + Vss per (paradigm, Iinj).
        for paradigm, stim in self.stimuli.items():
            cluster_name = mc.cluster_for_paradigm[paradigm]
            cid = name_to_id[cluster_name]
            self.cluster_assignment[paradigm] = cid
            self.cluster_is_manual[cid] = True

            self.Vss[paradigm] = {}
            for Iinj in np.squeeze(stim.Iinj, axis=-1):
                Iinj_f = float(Iinj)
                self.Vss[paradigm][Iinj_f] = mc.lookup_vss(paradigm, Iinj_f)

        # OLS for (Er, Rin) per cluster using manual (Iinj, Vss) pairs.
        for cname, cid in name_to_id.items():
            Iinjs_arr = np.asarray(list(mc.Vss_by_cluster[cname].keys()),   dtype=np.float64)
            Vsss_arr  = np.asarray(list(mc.Vss_by_cluster[cname].values()), dtype=np.float64)
            if Iinjs_arr.size >= 2 and np.std(Iinjs_arr) > 0:
                Er_hat, Rin_hat = self._ols_Er_Rin(Iinjs_arr, Vsss_arr)
            else:
                # Single Iinj level in the manual data: can't get Rin. This is a
                # user-data issue; warn and fall back to NaN so downstream code
                # surfaces the problem.
                print(
                    f"  Warning: cluster '{cname}' has only one Iinj level in "
                    f"manual data; cannot estimate Rin. Provide steady-state "
                    f"data at >=2 Iinj levels per cluster."
                )
                Er_hat = float(Vsss_arr[0]) if Vsss_arr.size > 0 else float("nan")
                Rin_hat = float("nan")
            self.Er_by_cluster[cid]  = Er_hat
            self.Rin_by_cluster[cid] = Rin_hat

        # Per-stimulus Er/Rin and Vss are not estimated in manual mode (we
        # didn't run weighted-mode on the xlsx data). Leave them empty.
        # _refit_clusters_pooled's bookkeeping is similarly unused.

    def _estimate_Vss_per_stimulus(self) -> None:
        """
        Per-stimulus Vss = median of the first Vss_duration_seconds of raw Vm
        samples at each clamp. Assumes the trace is in equilibrium for at
        least that long before any stimulus arrives. STOPGAP -- replaces the
        previous weighted-mode-on-full-trace estimator, which was unreliable
        when traces lacked a clean steady-state period.

        Uses RAW Vm rather than filtered Vm: the start-of-trace filter
        transient is exactly what we're avoiding. The median is robust to
        noise and single-sample outliers without needing filtering.
        """
        n_vss_samples: int = max(1, int(round(self.Er_Rin_cfg.Vss_duration_seconds / self.dt)))

        for paradigm, stim in self.stimuli.items():
            self.Vss_per_stimulus[paradigm] = {}
            for k, Iinj in enumerate(np.squeeze(stim.Iinj, axis=-1)):
                Vss_hat: float = float(np.median(stim.Vm[k, :n_vss_samples]))
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
        Greedy agglomerative clustering by Vss-pooling-residual constraint.

        Per-stimulus Vss is the 5-ms-median estimate, Vss_obs[paradigm, Iinj].
        When a cluster pools across paradigms, the cluster's pooled Vss at each
        Iinj is the median across paradigms-at-that-Iinj, Vss_pooled[cluster,
        Iinj]. The residual we care about is the difference these induce for
        each (paradigm, Iinj) in the cluster:

            r[paradigm, Iinj] = | Vss_pooled[cluster, Iinj] - Vss_obs[paradigm, Iinj] |

        A merge is acceptable only if every such residual is within
        max_residual_volts. This isolates the question of whether two paradigms
        share a leak state from the question of whether an OLS line fits them
        well (slope mismatch alone doesn't disqualify, but per-Iinj disagreement
        does).

        Algorithm:
          - Start with K = N (every paradigm its own cluster).
          - Repeatedly find the pair of clusters whose merge has the smallest
            "merged max pooling residual" across (paradigm, Iinj) in the union.
          - If the best-merge residual is <= max_residual_volts, merge and continue.
          - Otherwise, stop.

        K = N is always feasible: a singleton cluster's pooled-at-Iinj IS the
        paradigm's own observation, so the residual is 0 everywhere.

        Edge case: a (paradigm, Iinj) that is the only observation at that Iinj
        in the merged cluster contributes a trivially-zero residual. This is
        correct -- it has no peer to disagree with at that clamp.

        Sets self.cluster_assignment and self.cluster_is_manual (all False
        in this code path). Non-fittable paradigms (NaN Er/Rin) are placed
        in the largest cluster after the agglomerative step.
        """
        thresh: float = self.Er_Rin_cfg.max_residual_volts

        # Only paradigms with a complete Vss-per-stimulus dict can participate
        # in clustering.
        fittable: List[str] = [
            p for p, (Er, Rin) in self.Er_Rin_per_stimulus.items()
            if np.isfinite(Er) and np.isfinite(Rin)
        ]

        def pooling_max_residual(paradigms: List[str]) -> float:
            """
            For a candidate merged cluster:
              - At each Iinj used in this cluster, pool the first-5ms Vm
                samples across all paradigms that have that clamp, take the
                median -- this is the exact Vss the cluster will receive after
                _refit_clusters_pooled. Compare against each paradigm's own
                per-stimulus Vss (also the median of its own first-5ms samples).
              - Residual = |pooled - observed| for every (paradigm, Iinj) in
                the cluster.
              - Return the max residual across all pairs.

            A singleton cluster gets residual 0 (paradigm's pooled-at-Iinj IS
            its own observation). An Iinj observed by only one paradigm in the
            cluster also contributes residual 0 for that paradigm.

            Uses the same first-5ms raw-Vm concatenate-then-median as
            _refit_clusters_pooled, so the cluster acceptance criterion is
            exactly the residual that pooling will induce -- no approximation.
            """
            n_vss_samples: int = max(1, int(round(self.Er_Rin_cfg.Vss_duration_seconds / self.dt)))

            # Group raw first-5ms Vm samples by Iinj: Iinj -> [(paradigm, samples), ...]
            samples_by_iinj: Dict[float, List[Tuple[str, np.ndarray]]] = {}
            for p in paradigms:
                stim = self.stimuli[p]
                for clamp_idx, Iinj in enumerate(np.squeeze(stim.Iinj, axis=-1)):
                    samples_by_iinj.setdefault(float(Iinj), []).append(
                        (p, stim.Vm[clamp_idx, :n_vss_samples])
                    )

            max_resid: float = 0.0
            for Iinj, entries in samples_by_iinj.items():
                pooled_samples = np.concatenate([s for _, s in entries])
                pooled_vss = float(np.median(pooled_samples))
                for paradigm, samples in entries:
                    obs_vss = float(np.median(samples))  # same value as Vss_per_stimulus
                    r = abs(pooled_vss - obs_vss)
                    if r > max_resid:
                        max_resid = r
            return max_resid

        # Agglomerative loop.
        partition: List[List[str]] = [[p] for p in fittable]

        while len(partition) >= 2:
            best_pair: Optional[Tuple[int, int]] = None
            best_res:  float = float("inf")
            for i in range(len(partition)):
                for j in range(i + 1, len(partition)):
                    merged = partition[i] + partition[j]
                    res = pooling_max_residual(merged)
                    if res < best_res:
                        best_res = res
                        best_pair = (i, j)
            if best_pair is None or best_res > thresh:
                break
            i, j = best_pair
            merged = partition[i] + partition[j]
            partition = [c for idx, c in enumerate(partition) if idx not in (i, j)]
            partition.append(merged)

        # Assign cluster ids in dense order.
        self.cluster_assignment = {}
        self.cluster_is_manual = {}
        for k, members in enumerate(partition):
            self.cluster_is_manual[k] = False
            for p in members:
                self.cluster_assignment[p] = k

        # Place non-fittable stimuli in the largest cluster (or 0 if none exist).
        unassigned = [p for p in self.stimuli if p not in self.cluster_assignment]
        if unassigned:
            if partition:
                idx_largest = int(np.argmax([len(c) for c in partition]))
                majority = idx_largest
            else:
                majority = 0
                self.cluster_is_manual[0] = False
            for p in unassigned:
                self.cluster_assignment[p] = majority

    def _refit_clusters_pooled(self) -> None:
        """
        For each cluster, pool the first Vss_duration_seconds of raw Vm samples
        across all cluster members at each Iinj level, then take the median to
        get the cluster's authoritative Vss at that Iinj. OLS over (Iinj, Vss)
        pairs gives the cluster's (Er, Rin). Each paradigm then inherits its
        cluster's Vss values at its own Iinj levels.

        This runs only in unsupervised mode. When manual cluster data is
        provided, Vss / cluster fit comes from there directly and this
        method is skipped.

        STOPGAP version: previously this method pooled FULL filtered Vm traces
        and ran the weighted-mode estimator. That estimator was unreliable for
        traces without long clean baselines, so we now pool the first
        Vss_duration_seconds of raw Vm (same equilibrium assumption as the
        per-stimulus step) and take the median.
        """
        n_vss_samples: int = max(1, int(round(self.Er_Rin_cfg.Vss_duration_seconds / self.dt)))
        cluster_ids = sorted(set(self.cluster_assignment.values()))
        # Per cluster: Iinj -> pooled Vss_hat
        Vss_by_cluster: Dict[int, Dict[float, float]] = {}

        for k in cluster_ids:
            paradigms_in_k = [p for p, c in self.cluster_assignment.items() if c == k]

            # Pool the first-5ms RAW Vm samples by Iinj across this cluster's paradigms.
            Vm_head_by_Iinj: Dict[float, List[np.ndarray]] = {}
            for p in paradigms_in_k:
                stim = self.stimuli[p]
                for clamp_idx, Iinj in enumerate(np.squeeze(stim.Iinj, axis=-1)):
                    Vm_head_by_Iinj.setdefault(float(Iinj), []).append(
                        stim.Vm[clamp_idx, :n_vss_samples]
                    )

            Vss_by_cluster[k] = {}
            for Iinj, head_segments in Vm_head_by_Iinj.items():
                pooled_head = np.concatenate(head_segments)
                Vss_by_cluster[k][Iinj] = float(np.median(pooled_head))

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
        Estimate (Ee, Ei) from quantiles of the pooled Eeff(t) distribution,
        with a hard constraint on the physiological difference Ee - Ei.

        Operational assumptions (see project notes / EeEiCfg docstring):
          - Δ notation is retained: Δge, Δgi are signed deviations from a
            prestimulus baseline whose existence we acknowledge but do not
            attempt to estimate from this dataset.
          - We assume decreases from baseline are artifactual (disinhibition
            and disexcitation contribute negligibly). This is the operational
            stance for this paper, supported by prior duration-coder findings
            in the same brain region (Rose et al. 2016). Negative Δg after
            Iact correction is analyzed for correlation with active membrane
            events as supporting evidence.
          - Absolute reversal potentials are not recoverable per cell because
            the recordings have multiple uncontrolled voltage offsets
            (uncertain LJP protocol, electrode half-cell potentials, etc.).
            We constrain only the DIFFERENCE Ee - Ei, which is offset-
            invariant and is the only quantity that affects the downstream
            Δg inversion.

        Algorithm:
          1. Pool Eeff(t) across stimuli.
          2. Ee_hat = high quantile, Ei_hat = low quantile (rationale: when
             Δgi ≈ 0 and Δge dominates, Eeff ≈ Ee; converse for low quantile).
          3. Floor Ee_hat at max(Er_by_cluster). (Ee must be above resting V.)
          4. Enforce dE_min ≤ Ee_hat - Ei_hat ≤ dE_max by symmetric adjustment.
             - If Ee - Ei < dE_min: expand to dE_min, splitting the deficit.
             - If Ee - Ei > dE_max: contract to dE_max, splitting the excess.
             - This preserves the midpoint and modifies only the spread.

        Returns (Ee, Ei) in volts.
        """
        cfg = self.Ee_Ei_cfg

        # Step 1+2: quantile-based initial estimate.
        Eeff_pool: List[float] = []
        for stimulus in self.stimuli.values():
            Eeff_pool.extend(stimulus.timeseries["Eeff filtered"])

        Ei_hat: float = float(np.nanquantile(Eeff_pool, cfg.low_quantile))
        Ee_hat: float = float(np.nanquantile(Eeff_pool, cfg.high_quantile))

        # Step 3: Ee must be ≥ Er for every paradigm; use the most depolarized
        # cluster Er as a floor.
        Er_floor: float = float(np.nanmax(list(self.Er_by_cluster.values())))
        Ee_hat = max(Ee_hat, Er_floor)

        # Step 4: enforce physiological difference constraint Ee - Ei ∈ [dE_min, dE_max].
        # Symmetric adjustment: shift each reversal by half the deficit/excess,
        # in opposite directions, preserving the midpoint.
        dE: float = Ee_hat - Ei_hat
        if dE < cfg.dE_min:
            half_deficit: float = (cfg.dE_min - dE) / 2.0
            Ee_hat += half_deficit
            Ei_hat -= half_deficit
        elif dE > cfg.dE_max:
            half_excess: float = (dE - cfg.dE_max) / 2.0
            Ee_hat -= half_excess
            Ei_hat += half_excess

        return Ee_hat, Ei_hat

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
                stimulus.compute_dg_stats()

            if verbose:
                # Per-stimulus stats table. All values in nS.
                print("Per-stimulus statistics (nS):")
                hdr = f"  {'paradigm':<12} {'mean Δge':>10} {'mean Δgi':>10} {'net Δge':>10} {'net Δgi':>10}"
                print(hdr)
                print("  " + "-" * (len(hdr) - 2))
                for paradigm, stim in self.stimuli.items():
                    s = stim.stats
                    print(
                        f"  {paradigm:<12} "
                        f"{s['mean_dge']*1e9:>10.3f} "
                        f"{s['mean_dgi']*1e9:>10.3f} "
                        f"{s['net_dge']*1e9:>10.3f} "
                        f"{s['net_dgi']*1e9:>10.3f}"
                    )

            # Cm quality diagnostic. Runs last so all conductance estimates are
            # available. Cheap (one Spearman per stimulus on Δgsyn vs cross-clamp
            # |dV/dt| sensitivity); always run, never gated, since the result
            # is informative even when benign.
            cm_result = self.evaluate_Cm_quality()
            self.Cm_rho_all       = cm_result["rho_all"]
            self.Cm_rho_negatives = cm_result["rho_negatives"]

            if verbose:
                print(
                    f"Cm diagnostic: "
                    f"ρ_negatives = {self.Cm_rho_negatives:+.3f}, "
                    f"ρ_all = {self.Cm_rho_all:+.3f}"
                )
                # Soft warning thresholds. Calibration from synthetic ground-truth
                # tests in the docstring: |ρ_negatives| ~ 0.04 when Cm is correct,
                # ~0.8 when Cm is off 2x. Pick 0.3 as a midpoint that suggests
                # something is off without crying wolf on small biases.
                if np.isfinite(self.Cm_rho_negatives) and abs(self.Cm_rho_negatives) > 0.3:
                    print(
                        f"  ⚠ |ρ_negatives| > 0.3 suggests Cm may be misspecified "
                        f"(or other PSP-timescale model error). Consider checking Cm "
                        f"via membrane time constant from a hyperpolarizing pulse."
                    )

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

        # Diagnostic: accumulate filtered and raw losses across all evaluations.
        # See _evaluate_active_params docstring. Reported at the end of the
        # optimization to compare the two loss formulations.
        all_loss_filtered: List[float] = []
        all_loss_raw:      List[float] = []

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

                for theta, val_filt, val_raw in results:
                    # Optimizer uses the filtered loss (current behavior).
                    if val_filt < best_val:
                        best_val = val_filt
                        best_theta = theta
                    # Diagnostic: track both losses at every evaluation so we
                    # can compare them at the end of the optimization.
                    all_loss_filtered.append(val_filt)
                    all_loss_raw.append(val_raw)

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

        # ---- Diagnostic: compare filtered vs raw SSE-negative-Δgsyn loss ----
        # This is a temporary diagnostic to determine whether the choice of
        # filtered Δgsyn (current behavior) materially affects the optimization
        # compared to using raw Δgsyn. If the two losses are highly correlated
        # across evaluations and their values at the chosen optimum are close,
        # the filtering doesn't bias the optimization. If they diverge, the
        # current loss may be tracking filter artifacts rather than physical
        # violations.
        if len(all_loss_filtered) >= 2:
            lf = np.asarray(all_loss_filtered, dtype=np.float64)
            lr = np.asarray(all_loss_raw,      dtype=np.float64)
            finite = np.isfinite(lf) & np.isfinite(lr)
            if finite.sum() >= 2 and lf[finite].std() > 0 and lr[finite].std() > 0:
                pearson_r: float = float(np.corrcoef(lf[finite], lr[finite])[0, 1])
            else:
                pearson_r = float("nan")
            # Loss values at the chosen optimum
            best_filt = float(best_val)
            best_raw_at_best = float(lr[int(np.argmin(lf))]) if finite.sum() > 0 else float("nan")
            ratio = best_raw_at_best / best_filt if best_filt > 0 else float("nan")

            print(
                f"  loss diagnostic (n={len(all_loss_filtered)} evals): "
                f"filtered vs raw Pearson r = {pearson_r:.3f}; "
                f"at chosen optimum filtered={best_filt:.3e}, raw={best_raw_at_best:.3e} "
                f"(raw/filtered = {ratio:.2f})"
            )

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

        # (* 1e-3 is to scale Vm from millivolts to volts)
        self.Vm: np.ndarray = data[Iinj_colnames].to_numpy(dtype=np.float64).T * 1e-3       # units: Volts, shape: [Nclamps, Nsamples]
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

        # Per-stimulus scalar statistics, populated by compute_dg_stats() at
        # the end of run_analysis. Keys: "mean_dge", "mean_dgi", "net_dge",
        # "net_dgi" (all in siemens). See compute_dg_stats docstring for
        # definitions.
        self.stats: Dict[str, float] = {}

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

        self.timeseries["Eeff"] = Eeff
        self.timeseries["Eeff filtered"] = Eeff_filtered
        self.timeseries["dgsyn"] = dgsyn_safe
        self.timeseries["dgsyn filtered"] = dgsyn_filtered
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

    def compute_dg_stats(self) -> None:
        """
        Compute the four scalar summary statistics for this stimulus and store
        them in self.stats:

          mean_dge = ⟨[Δge]₊⟩_t
          mean_dgi = ⟨[Δgi]₊⟩_t
          net_dge  = ⟨[[Δge]₊ − [Δgi]₊]₊⟩_t
          net_dgi  = ⟨[[Δge]₊ − [Δgi]₊]₋⟩_t  (i.e. absolute value of the negative part)

        Conventions:
          - All four use the filtered Δge and Δgi.
          - Operates over all timepoints in the trace (no separate "response
            window" concept; the methods phrase "mean of their time courses
            for the duration of response" is interpreted as the whole trace).
          - Three rectifications are stacked in the net quantities: each
            input conductance individually (the inner [·]₊), then the
            positive/negative parts of their difference (the outer [·]₊
            and [·]₋). This matches the existing lab implementation. The
            current paper assumes decreases from baseline are artifactual,
            so rectifying inputs at zero is consistent with that assumption.
          - Units: siemens.
        """
        dge: np.ndarray = np.asarray(self.timeseries["dge filtered"], dtype=np.float64)
        dgi: np.ndarray = np.asarray(self.timeseries["dgi filtered"], dtype=np.float64)

        # Inner rectification: [·]₊ on each input.
        dge_pos: np.ndarray = np.maximum(dge, 0.0)
        dgi_pos: np.ndarray = np.maximum(dgi, 0.0)

        # Outer rectification: positive and negative parts of (dge_pos - dgi_pos).
        diff: np.ndarray = dge_pos - dgi_pos
        net_e: np.ndarray = np.maximum(diff, 0.0)         # positive part
        net_i: np.ndarray = np.maximum(-diff, 0.0)        # |negative part|

        self.stats["mean_dge"] = float(np.nanmean(dge_pos))
        self.stats["mean_dgi"] = float(np.nanmean(dgi_pos))
        self.stats["net_dge"]  = float(np.nanmean(net_e))
        self.stats["net_dgi"]  = float(np.nanmean(net_i))


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

        # Build a recording-level color map indexed by Iinj. This way, the same
        # current clamp value gets the same color across all paradigm subplots,
        # making it easier to track a single clamp's behavior across stimuli.
        # Default matplotlib qualitative cycle (tab10) gives 10 distinct colors;
        # we cycle if there are more unique Iinj values than that.
        all_Iinj = recording._Iinj_sorted_for_opt
        default_cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
        if not default_cycle:
            default_cycle = [f"C{i}" for i in range(10)]
        Iinj_color_map: Dict[float, str] = {
            float(I): default_cycle[i % len(default_cycle)]
            for i, I in enumerate(all_Iinj)
        }

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

            # Per-clamp colors keyed by this paradigm's Iinj values.
            colors = [Iinj_color_map[float(I)] for I in stimulus.Iinj.flatten()]

            for j in range(stimulus.Nclamps):
                axs[0, idx].plot(times, filtered_Vm[j, :], color=colors[j])

            Vss = [recording.Vss[paradigm][float(i)] for i in stimulus.Iinj.flatten()]
            for j in range(stimulus.Nclamps):
                axs[0, idx].plot([times[0], times[-1]], [Vss[j]] * 2, color="grey", ls="--")

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
            axs[2, idx].grid(True)

            axs[3, idx].plot(times, dge, c="r", label="Δge")
            axs[3, idx].plot(times, dgi, c="b", label="Δgi")
            # Δgsyn: fill between 0 and the trace, semitransparent black.
            # Reads as a "shadow" of the total synaptic drive against the
            # red/blue components -- shaded area = total Δgsyn magnitude.
            axs[3, idx].fill_between(times, 0, dgsyn, color="k", alpha=0.2, label="Δgsyn")
            axs[3, idx].plot(times, times * 0, "--k", linewidth=1)
            axs[3, idx].grid(True)

            Iact_pred: np.ndarray = np.stack(stimulus.timeseries["Iact"]).astype(np.float64).T  # type: ignore
            for j in range(stimulus.Nclamps):
                axs[4, idx].plot(times, Iact_pred[j, :], color=colors[j])
            axs[4, idx].plot(times, times * 0, "k--")
            axs[4, idx].grid(True)

            axs[4, idx].set_xlabel("time (s)")

        # Y-axis labels only on the leftmost column (rows share the y axis via
        # sharey="row", so this is purely cosmetic -- avoids redundant labels).
        axs[0, 0].set_ylabel("Vm (V)")
        axs[1, 0].set_ylabel("Eeff (V)")
        axs[2, 0].set_ylabel("Im (A)")
        axs[3, 0].set_ylabel("ΔG (S)")
        axs[4, 0].set_ylabel("Iact (A)")

        # Figure-level legend showing which color corresponds to which Iinj
        # (in nA, the most readable unit for typical patch-clamp injections).
        # Drawn as a horizontal strip below the suptitle so it doesn't crowd
        # any subplot. Each entry is a short colored line + Iinj value.
        legend_handles = [
            Line2D([0], [0], color=Iinj_color_map[float(I)], lw=2, label=f"{float(I)*1e9:.3f} nA")
            for I in all_Iinj
        ]
        fig.legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.97),
            ncol=min(len(legend_handles), 8),
            frameon=False,
            fontsize=9,
        )

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
            case_dir: Path = path_to_spreadsheet.parent

            # Read the xlsx first so we have the paradigm/Iinj structure needed
            # to validate optional manual cluster data.
            rdr: XLReader = XLReader(path_to_spreadsheet)
            stimuli: Dict[str, pd.DataFrame] = {}
            for paradigm in rdr.get_paradigms():
                df = rdr.get_paradigm_data(paradigm)
                aux = {"times", "stimulus", "representative"}
                n_iinj_cols = sum(1 for c in df.columns if c not in aux)
                if n_iinj_cols == 0:
                    print(
                        f"  Warning: paradigm '{paradigm}' in "
                        f"{path_to_spreadsheet.name} has no enabled Iinj columns "
                        f"(all 'use data' flags False); skipping."
                    )
                    continue
                stimuli[paradigm] = df
            assert len(stimuli) > 0, (
                f"No usable paradigms in {path_to_spreadsheet.name} -- every "
                f"paradigm has all 'use data' flags set to False."
            )

            # Build paradigm -> Iinj-list for validation.
            paradigm_iinjs: Dict[str, List[float]] = {}
            aux = {"times", "stimulus", "representative"}
            for paradigm, df in stimuli.items():
                paradigm_iinjs[paradigm] = [
                    float(c) for c in df.columns if c not in aux
                ]

            # Load optional per-case manual cluster data (cluster_assignments.json
            # + steady_states/ in the case folder). If present, the pipeline
            # uses it directly. If absent, falls back to unsupervised. If only
            # one of the two is present, raises (error message in from_case_dir).
            manual_cluster_data: Optional[ManualClusterDataCfg] = (
                ManualClusterDataCfg.from_case_dir(case_dir, paradigm_iinjs)
            )
            if manual_cluster_data is not None:
                n_clusters = len(set(manual_cluster_data.cluster_for_paradigm.values()))
                print(
                    f"  Loaded manual cluster data: {n_clusters} cluster(s) "
                    f"covering {len(manual_cluster_data.cluster_for_paradigm)} paradigm(s)."
                )
            manual_hash: str = (
                manual_cluster_data.hash_for_cache() if manual_cluster_data is not None
                else "no-manual-data"
            )

            # Cache lookup -- AFTER manual data is loaded, since the manual
            # data hash is part of cache invalidation.
            cached_params: Optional[Dict[str, Any]] = None
            cache_path: Optional[Path] = None
            if cache_dir is not None:
                cache_path = cache_dir / f"{path_to_spreadsheet.stem}.json"
                cached_params = self._try_load_cache(
                    cache_path, path_to_spreadsheet, cfg_hash, manual_hash
                )

            recording: WholeCellRecording = WholeCellRecording(
                parameters=rdr.get_paradigm_parameters(),
                stimuli=stimuli,
                filters=filters,
                n_workers=self.cfg.compute.n_workers,
                Er_Rin_cfg=self.cfg.Er_Rin_estimation,
                Ee_Ei_cfg=self.cfg.Ee_Ei_estimation,
                Et_opt_cfg=self.cfg.Et_optimization,
                numerics_cfg=self.cfg.numerics,
                manual_cluster_data=manual_cluster_data,
                cached_params=cached_params,
            )
            del rdr

            recording.run_analysis()

            # Save cache after a successful run (only when we just computed it).
            if cache_dir is not None and cached_params is None and cache_path is not None:
                self._save_cache(cache_path, recording, path_to_spreadsheet, cfg_hash, manual_hash)

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
        manual_hash: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Returns the cached params dict if the cache is valid (file exists, cfg
        hash matches, spreadsheet mtime matches, manual_cluster_data hash
        matches, and cache_invalidate is False). Otherwise returns None.
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

        cached_cfg_hash = payload.get("cfg_hash")
        cached_mtime    = payload.get("spreadsheet_mtime")
        # Old cache files (predating the manual cluster data system) may not
        # have this key. Treat as mismatch -- they should be regenerated.
        cached_manual_hash = payload.get("manual_hash")
        actual_mtime = spreadsheet_path.stat().st_mtime

        if cached_cfg_hash != cfg_hash:
            print(f"  cache: cfg_hash mismatch (cached={cached_cfg_hash}, current={cfg_hash}); will recompute.")
            return None
        if cached_mtime is None or abs(float(cached_mtime) - actual_mtime) > 1e-6:
            print(f"  cache: spreadsheet mtime changed; will recompute.")
            return None
        if cached_manual_hash != manual_hash:
            print(
                f"  cache: manual cluster data changed "
                f"(cached={cached_manual_hash}, current={manual_hash}); will recompute."
            )
            return None

        print(f"  cache: HIT -> {cache_path}")
        return payload["params"]

    def _save_cache(
        self,
        cache_path: Path,
        recording: WholeCellRecording,
        spreadsheet_path: Path,
        cfg_hash: str,
        manual_hash: str,
    ) -> None:
        """Write recording params + cache metadata to disk."""
        payload = {
            "cfg_hash":          cfg_hash,
            "manual_hash":       manual_hash,
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