import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 

from scipy.optimize import nnls
from pyhelpers.store import save_fig
from pathlib import Path

from libs.readers import XLReader, AnalyzerCfg

from typing import Dict, Optional, List

class WholeCellRecording:
    def __init__(self, data: pd.DataFrame, parameters: pd.DataFrame, current_clamps: Optional[List[float]]) -> None:
        self.data: pd.DataFrame = data
        self.parameters: pd.DataFrame = parameters

        # Mask parameters to only include those corresponding to selected usable Iinj levels
        if current_clamps is not None:
            assert isinstance(current_clamps, list)
            assert all(isinstance(x, float) for x in current_clamps)
            assert all(x in parameters["Iinj"] for x in self.parameters) 
            self.parameters = parameters[parameters["Iinj"].isin(current_clamps)]

        # Scale voltages from millivolts to volts
        for iinj in self.parameters["Iinj"]:
            self.data[f"{iinj:.3d}"] *= 1e3
            
    def estimate_conductances(self, bin_s: float = 5e-3):
        """
        Estimate excitatory and inhibitory synaptic conductances from whole-cell
        current-clamp recordings using a bin-integrated, voltage-domain formulation.

        --------------------------------------------------------------------------
        1. Biophysical model
        --------------------------------------------------------------------------

        For each injected current level i, the membrane potential v_i(t) is assumed
        to obey a single-compartment current balance equation:

            C_m * dv_i(t)/dt =
                g_e(t) * (E_e - v_i(t))
            + g_i(t) * (E_i - v_i(t))
            + g_l * (E_r - v_i(t))
            + I_inj,i

        where:
            C_m     membrane capacitance (F)
            g_e(t)  excitatory synaptic conductance (S)
            g_i(t)  inhibitory synaptic conductance (S)
            g_l     leak conductance (S)
            E_e     excitatory reversal potential (V)
            E_i     inhibitory reversal potential (V)
            E_r     resting (leak) reversal potential (V)
            I_inj,i injected current for clamp i (A)

        The synaptic conductances g_e(t) and g_i(t) are assumed to be identical across
        current-clamp levels for a given stimulus (i.e., the neuron receives the same
        synaptic input regardless of injected current), and to be non-negative.

        No assumptions are made about spike-generating currents; the analysis is
        intended for subthreshold membrane dynamics.

        --------------------------------------------------------------------------
        2. Numerical formulation and solution
        --------------------------------------------------------------------------

        Time is partitioned into contiguous bins of duration bin_s. Within each bin
        k, synaptic conductances are assumed to be constant:

            g_e(t) = g_e^(k),   g_i(t) = g_i^(k)    for t in bin k.

        The membrane equation is integrated over each bin [t_k, t_{k+1}] for each
        current clamp i, yielding:

            C_m * (v_i(t_{k+1}) - v_i(t_k)) =
                g_e^(k) * ∫_{t_k}^{t_{k+1}} (E_e - v_i(t)) dt
            + g_i^(k) * ∫_{t_k}^{t_{k+1}} (E_i - v_i(t)) dt
            + ∫_{t_k}^{t_{k+1}} [ g_l (E_r - v_i(t)) + I_inj,i ] dt

        Rearranging gives a linear system for each bin k:

            y_{i,k} = g_e^(k) * A_{e,i,k} + g_i^(k) * A_{i,i,k}

        where:
            y_{i,k}     = C_m * Δv_{i,k}
                        - ∫_{t_k}^{t_{k+1}} [ g_l (E_r - v_i(t)) + I_inj,i ] dt
            A_{e,i,k}   = ∫_{t_k}^{t_{k+1}} (E_e - v_i(t)) dt
            A_{i,i,k}   = ∫_{t_k}^{t_{k+1}} (E_i - v_i(t)) dt

        For each bin, the system is solved across all current clamps i using
        non-negative least squares (NNLS):

            minimize || X_k g_k - y_k ||_2
            subject to g_k >= 0

        where g_k = [g_e^(k), g_i^(k)]^T and X_k has columns A_e and A_i.

        All integrals are computed using prefix (cumulative) trapezoidal integration
        on the uniformly sampled voltage traces. This avoids numerical differentiation
        of v(t) and yields stable estimates of bin-wise conductances.

        The resulting bin-level conductances are expanded back to the original time
        grid as piecewise-constant time series.

        In addition to conductances, the function computes:
            - the leak current I_l(t) = g_l (E_r - v(t))
            - a bin-consistent estimate of capacitive membrane current
            I_m(t) = C_m * Δv / Δt_bin
            - a forward-simulated membrane potential V_pred(t), obtained by
            numerically integrating the membrane equation using the inferred
            conductances

        Per-bin diagnostics are also stored, including:
            - the NNLS residual norm
            - the condition number of X_k^T X_k (a measure of identifiability of
            excitatory vs inhibitory contributions)

        These diagnostics can be used for model validation and for unsupervised
        selection of reversal potentials in higher-level optimization loops.

        --------------------------------------------------------------------------
        Assumptions and limitations
        --------------------------------------------------------------------------

        - Conductances are assumed constant within bins.
        - Synaptic inputs are assumed identical across current-clamp levels.
        - Conductances are constrained to be non-negative.
        - Reversal potentials are treated as fixed inputs to this function.
        - The method does not enforce sparsity or smoothness priors on conductances;
        any such criteria should be applied only at the hyperparameter selection
        stage.

        Parameters
        ----------
        bin_s : float
            Duration of each time bin in seconds.

        Returns
        -------
        self.data : pandas.DataFrame
            Updated data table containing estimated conductances, reconstructed
            currents, predicted membrane potential, and diagnostic time series.
        """

        """ === DATA FETCHING & SETUP === """
        # --- timebase ---
        dt: float = float(self.data["times"][1] - self.data["times"][0])
        assert bin_s >= dt, f"Bin size ({bin_s}) cannot be smaller than sampling rate ({dt})."

        # --- Experimental Controlled & Measured Variables: Current Injection & Membrane Potential ---
        Iinj: np.ndarray = np.asarray(self.parameters["Iinj"].to_numpy(), dtype=float) # shape: (Nclamps,), units: Amperes   
        Iinj_colnames: List[str] = [f"{x:.3e}" for x in Iinj]

        Vm: np.ndarray = self.data[Iinj_colnames].to_numpy() # shape: (Nsamples, Nclamps), units: Volts
        Nsamples, Nclamps = Vm.shape

        # --- PARAMTERS ---
        # --- Measurables ---
        Cm: float = float(self.parameters["Cm"][0])         # units: Farads
        gl: float = 1.0 / float(self.parameters["Rin"][0])  # units: Siemens
        Er: float = float(self.parameters["Er"][0])         # units: Volts
        
        # --- Hyperparameters ---
        Ee: float = float(self.parameters["Ee"][0])         # units: Volts
        Ei: float = float(self.parameters["Ei"][0])         # units: Volts

        # --- BINNING CONVENTIONS ---
        bin_len: int = int(round(bin_s / dt))
        left_edges: np.ndarray = np.arange(0, Nsamples, bin_len)                # shape: (Nbins,)
        right_edges: np.ndarray = np.minimum(left_edges + bin_len, Nsamples)    # shape: (Nbins,)
        Nbins: int = left_edges.size

        # --- Integration Helper ---
        def trapez_prefix_integral(arr: np.ndarray) -> np.ndarray:
            """
            In main loop, instead of computing each bin integral individually 
            (∫_{t_i}^{t_{i+1}} for each i) we compute the full integral function
            (prefix[k] = ∫_{t_k}^{t_{k+1}}) and compute bin integrals by subtraction
            ∫_{t_i}^{t_{i+1}} = prefix[i+1] - prefix[i]. Faster & simpler. 
            """
            area: np.ndarray = 0.5 * (arr[1:, :] + arr[:-1, :]) * dt    # shape: (Nsamples-1, Nclamps)
            pref: np.ndarray = np.zeros_like(arr)                       # shape: (Nsamples, Nclamps)
            pref[1:, :] = np.cumsum(area, axis=0)
            return pref
        

        """" === MAIN COMPUTATIONS === """
        # --- build needed integrands ---
        Il: np.ndarray = gl * (Er - Vm)                                 # shape: (Nsamples, Nclamps), units: Amperes
        Iinj: np.ndarray = np.tile(Iinj.reshape(1, -1), (Nsamples, 1))  # shape: (Nsamples, Nclamps), units: Amperes

        # --- prefix integrals ---
        pref_leak_inj: np.ndarray = trapez_prefix_integral(Il + Iinj)   # shape: (Nsamples, Nclamps), units: Coulombs
        pref_epotential: np.ndarray = trapez_prefix_integral(Ee - Vm)   # shape: (Nsamples, Nclamps), units: Webers (Volt * Second)
        pref_ipotential: np.ndarray = trapez_prefix_integral(Ei - Vm)   # shape: (Nsamples, Nclamps), units: Webers (Volt * Second)

        # --- bin integrals using prefix differences ---
        int_leak_inj: np.ndarray = pref_leak_inj[right_edges - 1, :] - pref_leak_inj[left_edges, :]         # shape: (Nbins, Nclamps), units: Coulombs
        int_epotential: np.ndarray = pref_epotential[right_edges - 1, :] - pref_epotential[left_edges, :]   # shape: (Nbins, Nclamps), units: Webers (Volt * Second)
        int_ipotential: np.ndarray = pref_ipotential[right_edges - 1, :] - pref_ipotential[left_edges, :]   # shape: (Nbins, Nclamps), units: Webers (Volt * Second)

        # --- Δv per bin per clamp ---
        dv: np.ndarray = Vm[right_edges - 1, :] - Vm[left_edges, :] # shape: (Nbins, Nclamps), units: Volts

        # --- y per bin per clamp ---
        y: np.ndarray = Cm * dv - int_leak_inj # (Nbins, Nclamps), units: Coulombs

        # --- solve per bin with NNLS ---
        ge_bins: np.ndarray = np.empty(Nbins)       # shape: (Nbins,), units: Siemens
        gi_bins: np.ndarray = np.empty(Nbins)       # shape: (Nbins,), units: Siemens
        resnorm_bins: np.ndarray = np.empty(Nbins)  # shape: (Nbins,), units: Coulombs
        cond_bins: np.ndarray = np.empty(Nbins)     # shape: (Nbins,)

        for k in range(Nbins):
            Xk: np.ndarray = np.column_stack([int_epotential[k, :], int_ipotential[k, :]])  # shape: (Nclamps, 2), units: Webers (Volt * Second)
            yk: np.ndarray = y[k, :]                                                        # shape: (Nclamps,), units: Coulombs

            # conditioning diagnostic
            XtX: np.ndarray = Xk.T @ Xk # shape: (2, 2), units: Webers^2 
            cond_bins[k] = np.linalg.cond(XtX) if np.all(np.isfinite(XtX)) else np.nan

            gk, rnorm = nnls(Xk, yk)
            ge_bins[k], gi_bins[k] = gk
            resnorm_bins[k] = rnorm
        

        """ SAVE RESULTS """
        # --- expand ge/gi to sample grid ---
        ge: np.ndarray = np.repeat(ge_bins, bin_len)[:Nsamples]  # shape: (Nsamples,), units: Siemens
        gi: np.ndarray = np.repeat(gi_bins, bin_len)[:Nsamples]  # shape: (Nsamples,), units: Siemens

        self.data["excitation"] = ge
        self.data["inhibition"] = gi

        # --- expand & save diagnostic trances ---
        self.data["bin_resnorm"] = np.repeat(resnorm_bins, bin_len)[:Nsamples]
        self.data["bin_cond_XtX"] = np.repeat(cond_bins, bin_len)[:Nsamples]

        # --- per-clamp leakage current ---
        for j, col in enumerate(Iinj_colnames):
            self.data[f"Il_{col}"] = Il[:, j]

        # --- per-clamp capacitive current ---
        Im_bins = Cm * dv / bin_s                               # shape: (Nbins, Nclamps), units: Amperes
        Im = np.repeat(Im_bins, bin_len, axis=0)[:Nsamples, :]  # shape: (Nsamples, Nclamps), units: Amperes
        for j, col in enumerate(Iinj_colnames):
            self.data[f"Im_{col}"] = Im[:, j]

        # --- forward-simulated Vpred per clamp (RK4; vectorized across clamps) ---
        Vpred: np.ndarray = np.empty_like(Vm) # shape: (Nsamples, Nclamps), units: Volts
        Vpred[0, :] = Vm[0, :] # Initial Conditions

        # RK4 algo
        for ti in range(Nsamples - 1):
            ge_t: float = ge[ti]
            gi_t: float = gi[ti]

            def f(vstate: np.ndarray) -> np.ndarray:
                return (
                    ge_t * (Ee - vstate) +
                    gi_t * (Ei - vstate) +
                    gl   * (Er - vstate) +
                    Iinj
                ) / Cm

            k1: np.ndarray = f(Vpred[ti, :])
            k2: np.ndarray = f(Vpred[ti, :] + 0.5 * dt * k1)
            k3: np.ndarray = f(Vpred[ti, :] + 0.5 * dt * k2)
            k4: np.ndarray = f(Vpred[ti, :] + dt * k3)

            Vpred[ti + 1, :] = Vpred[ti, :] + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

        for j, col in enumerate(Iinj_colnames):
            self.data[f"Vpred_{col}"] = Vpred[:, j]
        
class Analyzer:
    def __init__(self, cfg: AnalyzerCfg):
        self.cfg: AnalyzerCfg = cfg

    def plot_timeseries(
        self,
        recordings,
        filename: Path,
        filetype: str = "png",
        current_clamps: Optional[List[float]] = None,
        cond_warn: float = 1e6,
        cond_bad: float = 1e10,
        resnorm_factor_warn: float = 5.0,
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
        n_files: int = len(self.cfg.paths_to_spreadsheets)
        iinj_clamps_to_use: list = [None] * n_files if self.cfg.iinj_clamps_to_use is None else self.cfg.iinj_clamps_to_use
        for i in range(n_files):
            rdr: XLReader = XLReader(self.cfg.paths_to_spreadsheets[i])
            recordings: Dict[str, WholeCellRecording] = {}
            for paradigm in rdr.get_paradigms():
                recording: WholeCellRecording = WholeCellRecording(rdr.get_paradigm_data(paradigm), rdr.get_paradigm_parameters(paradigm), iinj_clamps_to_use[i])
                recording.estimate_conductances()
                recordings[paradigm] = recording
                print(f"{paradigm} done")
            
            self.plot_timeseries(recordings, self.cfg.image_save_dir / self.cfg.paths_to_spreadsheets[i].name, filetype=self.cfg.image_save_type, current_clamps=iinj_clamps_to_use)