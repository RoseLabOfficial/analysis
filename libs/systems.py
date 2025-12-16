import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import os

from scipy.optimize import nnls
from pyhelpers.store import save_fig
from pathlib import Path

from libs.readers import XLReader, AnalyzerCfg

from typing import Dict, Optional, List, Tuple

class WholeCellRecording:
    def __init__(self, data: pd.DataFrame, parameters: pd.DataFrame, current_clamps: Optional[List[float]]=None) -> None:
        self.data: pd.DataFrame = data
        self.parameters: pd.DataFrame = parameters

        if current_clamps is not None:
            assert isinstance(current_clamps, list)
            assert all(isinstance(x, float) for x in current_clamps)

            assert all(x in parameters["Iinj"] for x in self.parameters) 
            
            self.parameters = parameters[parameters["Iinj"].isin(current_clamps)]
    
    def compute_passive_conductances(self, bin_s: float = 5e-3):
        # --- timebase ---
        t = self.data["times"].to_numpy()
        dt = float(t[1] - t[0])
        if bin_s < dt:
            raise ValueError("bin_s must be >= dt")

        # --- clamp levels and voltage matrix ---
        Iinj: np.ndarray = np.asarray(self.parameters["Iinj"].to_numpy(), dtype=float)  # (Nclamps,)
        Iinj_colnames: List[str] = [f"{x:.3e}" for x in Iinj]

        v_mV = self.data[Iinj_colnames].to_numpy()          # (Nsamples, Nclamps) in mV
        v = v_mV * 1e-3                                     # convert to V
        Nsamples, Nclamps = v.shape

        # --- constants (assuming stored per clamp colname) ---
        Cm = float(self.parameters["Cm"][0])             # F
        gl = 1.0 / float(self.parameters["Rin"][0])      # S
        Er = float(self.parameters["Er"][0])             # V
        Ee = float(self.parameters["Ee"][0])             # V
        Ei = float(self.parameters["Ei"][0])             # V

        # --- binning indices ---
        bin_len = int(round(bin_s / dt))
        bin_len = max(1, bin_len)
        t0 = np.arange(0, Nsamples, bin_len)                 # (Nbins,)
        tf = np.minimum(t0 + bin_len, Nsamples)              # (Nbins,)
        Nbins = t0.size

        # --- helper: cumulative trapezoid integral along time for each clamp column ---
        # prefix[k] = integral from 0 to time index k (exclusive-ish; see below)
        def cumtrapz_prefix(arr: np.ndarray) -> np.ndarray:
            # arr: (Nsamples, Nclamps)
            # returns pref: (Nsamples, Nclamps) where pref[j] = ∫_{0}^{j} arr(t) dt with trapezoid
            # implement via trapezoid areas per interval
            # area[j] corresponds to interval (j-1 -> j)
            area = 0.5 * (arr[1:, :] + arr[:-1, :]) * dt              # (Nsamples-1, Nclamps)
            pref = np.zeros_like(arr)
            pref[1:, :] = np.cumsum(area, axis=0)
            return pref

        # Better: define pref such that pref[j] = ∫_0^{j} arr dt using samples 0..j with trapezoid
        # With the construction above, pref[idx] is ∫_0^{idx} arr dt over intervals up to idx.
        # Then ∫_{t0}^{tf-1} arr dt = pref[tf-1] - pref[t0]
        # For our bin integrals over [t0, tf) we want ∫_{t0}^{tf} (continuous) approx via trapezoid,
        # which corresponds to pref[tf-1] - pref[t0] plus the last half-interval isn't included.
        # Easiest/robust: compute pref over intervals and use tf-1 indexing; bins of length >=2 behave well.

        # --- build needed integrands ---
        leak = gl * (Er - v)                                  # (Nsamples, Nclamps)
        inj = np.tile(Iinj.reshape(1, -1), (Nsamples, 1))     # (Nsamples, Nclamps)

        Ae_arr = (Ee - v)                                     # (Nsamples, Nclamps)
        Ai_arr = (Ei - v)                                     # (Nsamples, Nclamps)

        # --- prefix integrals ---
        pref_leak_inj = cumtrapz_prefix(leak + inj)
        pref_Ae = cumtrapz_prefix(Ae_arr)
        pref_Ai = cumtrapz_prefix(Ai_arr)

        # --- bin integrals using prefix differences ---
        # ∫_{t0}^{tf-1} f dt approximated
        int_leak_inj = pref_leak_inj[tf - 1, :] - pref_leak_inj[t0, :]   # (Nbins, Nclamps)
        int_Ae = pref_Ae[tf - 1, :] - pref_Ae[t0, :]                     # (Nbins, Nclamps)
        int_Ai = pref_Ai[tf - 1, :] - pref_Ai[t0, :]                     # (Nbins, Nclamps)

        # --- Δv per bin per clamp ---
        dv = v[tf - 1, :] - v[t0, :]                                     # (Nbins, Nclamps)

        # --- y per bin per clamp ---
        y = Cm * dv - int_leak_inj                                       # (Nbins, Nclamps)

        # --- solve per bin with NNLS ---
        ge_bins = np.zeros(Nbins)
        gi_bins = np.zeros(Nbins)
        resnorm_bins = np.full(Nbins, np.nan)
        cond_bins = np.full(Nbins, np.nan)

        for k in range(Nbins):
            Xk = np.column_stack([int_Ae[k, :], int_Ai[k, :]])           # (Nclamps, 2)
            yk = y[k, :]                                                 # (Nclamps,)

            # conditioning diagnostic
            XtX = Xk.T @ Xk
            cond_bins[k] = np.linalg.cond(XtX) if np.all(np.isfinite(XtX)) else np.nan

            gk, rnorm = nnls(Xk, yk)
            ge_bins[k], gi_bins[k] = gk
            resnorm_bins[k] = rnorm

        # --- expand ge/gi to sample grid ---
        ge = np.repeat(ge_bins, bin_len)[:Nsamples]  # (Nsamples,)
        gi = np.repeat(gi_bins, bin_len)[:Nsamples]  # (Nsamples,)

        self.data["excitation"] = ge
        self.data["inhibition"] = gi

        self.data["bin_resnorm"] = np.repeat(resnorm_bins, bin_len)[:Nsamples]
        self.data["bin_cond_XtX"] = np.repeat(cond_bins, bin_len)[:Nsamples]

        # --- per-clamp leakage current ---
        Il = gl * (Er - v)  # (Nsamples, Nclamps), in A
        for j, col in enumerate(Iinj_colnames):
            self.data[f"Il_{col}"] = Il[:, j]

        # --- per-clamp "Im" as bin-consistent capacitive current (no dv/dt noise)
        # Im_bin[i,k] = C * Δv / Δt_bin, then expanded
        bin_durations = (tf - t0) * dt                 # (Nbins,)
        Im_bins = Cm * (dv / bin_durations[:, None])   # (Nbins, Nclamps), A
        Im = np.repeat(Im_bins, bin_len, axis=0)[:Nsamples, :]  # (Nsamples, Nclamps)
        for j, col in enumerate(Iinj_colnames):
            self.data[f"Im_{col}"] = Im[:, j]

        # --- forward-simulated Vpred per clamp (Euler; vectorized across clamps) ---
        Vpred = np.empty_like(v) # (Nsamples, Nclamps), Volts
        Vpred[0, :] = v[0, :]

        for ti in range(Nsamples - 1):
            ge_t = ge[ti]
            gi_t = gi[ti]

            def f(vstate):
                return (
                    ge_t * (Ee - vstate) +
                    gi_t * (Ei - vstate) +
                    gl   * (Er - vstate) +
                    Iinj
                ) / Cm

            k1 = f(Vpred[ti, :])
            k2 = f(Vpred[ti, :] + 0.5 * dt * k1)
            k3 = f(Vpred[ti, :] + 0.5 * dt * k2)
            k4 = f(Vpred[ti, :] + dt * k3)

            Vpred[ti + 1, :] = Vpred[ti, :] + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

        # save in mV to match original v columns (optional; just be consistent)
        Vpred_mV = Vpred * 1e3
        for j, col in enumerate(Iinj_colnames):
            self.data[f"Vpred_{col}"] = Vpred_mV[:, j]
        
    def get_clamp_near_0(self) -> Tuple[int, float]:
        index_of_minimum_injected_current: int = int(np.argmin(np.abs(self.parameters["Iinj"])))
        minimum_injected_current: float = self.parameters["Iinj"][index_of_minimum_injected_current]
        return index_of_minimum_injected_current, minimum_injected_current

    def estimate_conductances(self):
        self.compute_passive_conductances()
        return self.data

class Analyzer:
    def __init__(self, cfg: AnalyzerCfg):
        self.cfg: AnalyzerCfg = cfg
    
    def estimation_without_optim_activation_potential(self, recordings: Dict[str, WholeCellRecording]) -> Dict[str, WholeCellRecording]:
        assert len(recordings) > 0

        for paradigm in recordings:
            recordings[paradigm].estimate_conductances()
            print(f"{paradigm} done")

        return recordings

    def analyze(self, path_to_spreadsheet: Path, current_clamps: Optional[List[float]]=None):
        reader = XLReader(path_to_spreadsheet)
        recordings: Dict[str, WholeCellRecording] = {}
        for _, paradigm in enumerate(reader.get_paradigms()):
            recordings[paradigm] = WholeCellRecording(
                reader.get_paradigm_data(paradigm), 
                reader.get_paradigm_parameters(paradigm),
                current_clamps=current_clamps
            )

        recordings = self.estimation_without_optim_activation_potential(recordings)
       
        return recordings

    def plot_dev(
        self,
        recordings,
        filename: Path,
        filetype: str = "png",
        current_clamps: Optional[List[float]] = None,
        cond_warn: float = 1e8,
        cond_bad: float = 1e12,
        resnorm_factor_warn: float = 10.0,
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

        ncols = len(recordings)
        nrows = 6
        fig, axs = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            sharex="all",
            sharey="row",
            figsize=(15, 10),
            constrained_layout=True,
        )

        if ncols == 1:
            axs = np.expand_dims(axs, axis=1)

        for idx, paradigm in enumerate(recordings):
            rec = recordings[paradigm]
            paradigm_iinj: List[float] = list(rec.parameters["Iinj"])

            if current_clamps is not None and set(current_clamps) <= set(paradigm_iinj):
                paradigm_iinj = list(set(paradigm_iinj).intersection(current_clamps))
                assert len(paradigm_iinj) > 0, (
                    "Cannot plot. Specified current clamps have no intersection with "
                    "current clamps listed in parameters."
                )

            Iinj_cols = [f"{x:.3e}" for x in paradigm_iinj]
            times = rec.data["times"].to_numpy()

            # measured Vm
            v = rec.data[Iinj_cols].to_numpy() * 1e-3

            # predicted Vm
            vpred_cols = [f"Vpred_{c}" for c in Iinj_cols]
            vpred = rec.data[vpred_cols].to_numpy() * 1e-3

            # Im and Il
            im_cols = [f"Im_{c}" for c in Iinj_cols]
            il_cols = [f"Il_{c}" for c in Iinj_cols]
            Im = rec.data[im_cols].to_numpy()
            Il = rec.data[il_cols].to_numpy()

            # conductances
            ge = rec.data["excitation"].to_numpy()
            gi = rec.data["inhibition"].to_numpy()

            # diagnostics
            resnorm = rec.data["bin_resnorm"].to_numpy() if "bin_resnorm" in rec.data else None
            cond = rec.data["bin_cond_XtX"].to_numpy() if "bin_cond_XtX" in rec.data else None

            axs[0, idx].set_title(paradigm)

            # ---- Row 0: Vm + Vpred (same colors) ----
            curves = axs[0, idx].plot(times, v)  # solid Vm traces
            colors = [line.get_color() for line in curves]

            for j in range(v.shape[1]):
                axs[0, idx].plot(times, vpred[:, j], linestyle=":", color=colors[j])  # dotted Vpred

            # Pull reference voltages from parameters
            Er = float(rec.parameters["Er"][0])
            Ee = float(rec.parameters["Ee"][0])
            Ei = float(rec.parameters["Ei"][0])

            Ess = rec.parameters["Ess"]
            Ess = np.asarray(Ess.to_numpy(), dtype=float)
            # if Vss includes more clamps than we're plotting, subset by the paradigm_iinj indices
            # (assumes same ordering as rec.parameters["Iinj"])
            all_cols = [f"{x:.3e}" for x in list(rec.parameters["Iinj"])]
            idxs = [all_cols.index(c) for c in Iinj_cols]
            Ess = Ess[idxs]

            # Horizontal reference lines: Er, Ee, Ei (fixed colors)
            axs[0, idx].axhline(Er, linestyle="--", color="k", linewidth=1, label="Er" if idx == 0 else None)
            axs[0, idx].axhline(Ee, linestyle="--", color="r", linewidth=1, label="Ee" if idx == 0 else None)
            axs[0, idx].axhline(Ei, linestyle="--", color="b", linewidth=1, label="Ei" if idx == 0 else None)

            # Horizontal Vss lines: per clamp, same color as that clamp trace
            for j in range(v.shape[1]):
                axs[0, idx].axhline(Ess[j], linestyle="--", color=colors[j], linewidth=0.8)

            axs[0, idx].set_ylabel("Vm (V)")
            axs[0, idx].grid(True)
            if idx == 0:
                axs[0, idx].legend(loc="upper right")

            # ---- Row 1: Im ----
            for j in range(Im.shape[1]):
                axs[1, idx].plot(times, Im[:, j], color=colors[j])
            axs[1, idx].set_ylabel("Im (A)")
            axs[1, idx].grid(True)

            # ---- Row 2: Il ----
            for j in range(Il.shape[1]):
                axs[2, idx].plot(times, Il[:, j], color=colors[j])
            axs[2, idx].set_ylabel("Il (A)")
            axs[2, idx].grid(True)

            # ---- Row 3: conductances ----
            axs[3, idx].plot(times, ge, c="r", label="ge")
            axs[3, idx].plot(times, gi, c="b", label="gi")
            axs[3, idx].plot(times, times * 0, "--k", linewidth=1)
            axs[3, idx].set_ylabel("G (S)")
            axs[3, idx].grid(True)
            if idx == 0:
                axs[3, idx].legend(loc="upper right")

            # ---- Row 4: residual norm + warning threshold ----
            ax_r = axs[4, idx]
            ax_r.grid(True)
            ax_r.set_ylabel("resnorm")

            if resnorm is not None:
                ax_r.plot(times, resnorm, linewidth=1)

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

            if cond is not None:
                cond_plot = np.where(np.isfinite(cond) & (cond > 0), cond, np.nan)
                ax_c.plot(times, cond_plot, linewidth=1)

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

        plt.show()

    def run(self) -> None:
        for i in range(len(self.cfg.paths_to_spreadsheets)):
            iinj_clamps_to_use: Optional[List[float]] = self.cfg.iinj_clamps_to_use if self.cfg.iinj_clamps_to_use is None else self.cfg.iinj_clamps_to_use[i]
            recordings = self.analyze(
                self.cfg.paths_to_spreadsheets[i], 
                current_clamps=iinj_clamps_to_use
            )
            self.plot_dev(recordings, self.cfg.image_save_dir / self.cfg.paths_to_spreadsheets[i].name, filetype=self.cfg.image_save_type, current_clamps=iinj_clamps_to_use)