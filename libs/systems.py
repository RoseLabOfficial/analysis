import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import nnls
from pyhelpers.store import save_fig
from pathlib import Path

from libs.readers import XLReader, AnalyzerCfg
from libs.test import posterior_C_batch_grid, summarize_posterior_batch

from typing import Dict, Optional, List, Tuple

import numpy as np

def weighted_median(values: np.ndarray, weights: np.ndarray, quantile: float):    
    sort_indices = np.argsort(values)
    
    values_sorted = values[sort_indices]
    weights_sorted = weights[sort_indices]  

    cumsum = weights_sorted.cumsum()
    cutoff = weights_sorted.sum() * quantile
    
    return values_sorted[cumsum >= cutoff][0]

class WholeCellRecording:
    def __init__(self, parameters: pd.DataFrame, bin_s: float, stimuli: Dict[str, pd.DataFrame]) -> None:
        assert len(stimuli) > 0
        
        self.Cm: float = parameters["Cm"][0]        # units: Farads
        self.gl: float = 1 / parameters["Rin"][0]   # units: Siemens
        self.Er: float = parameters["Er"][0]        # units: Volts

        example_t: pd.Series = list(stimuli.values())[0]["times"] # units: Seconds, shape: (Nsamples,)
        self.dt: float = example_t[1] - example_t[0] # time assumed sampled at constant interval; units: Seconds

        self.bin_nsamples: int = round(bin_s / self.dt)
        self.bin_s: float = self.bin_nsamples * self.dt # units: Seconds

        self.Ee: float # units: Volts
        self.Ei: float # units: Volts

        self.stimuli: Dict[str, WholeCellStimulus] = {name: WholeCellStimulus(self, data) for name, data in stimuli.items()}

    def estimate_Ee_Ei(self) -> Tuple[float, float]:
        Eeff_pool: List[float] = []
        gsyn_pool: List[float] = []
        for stimulus in self.stimuli.values():
            Eeff_pool.extend(stimulus.binned_timeseries["Eeff unbiased"])
            gsyn_pool.extend(stimulus.binned_timeseries["gsyn"])

        Ei_hat: float = float(weighted_median(np.array(Eeff_pool), np.array(gsyn_pool)**2, 0.1)) # units: Volts
        Ee_hat: float = float(weighted_median(np.array(Eeff_pool), np.array(gsyn_pool)**2, 0.9)) # units: Volts

        return Ee_hat, Ei_hat
        
    def run_analysis(self):
        print("Estimating reversal potentials... ")
        for paradigm, stimulus in self.stimuli.items():
            print(f"\t...{paradigm}")
            stimulus.calculate_target_Qsyn()
            stimulus.estimate_Eeff()
        print("Done")
        self.Ee, self.Ei = self.estimate_Ee_Ei()

        print("Calculating synaptic conductances... ")
        for paradigm, stimulus in self.stimuli.items():
            print(f"\t...{paradigm}")
            stimulus.estimate_ge_gi()
            stimulus.calculate_predicted_Vm()
        print("Done")
        

class WholeCellStimulus:
    def __init__(self, recording: WholeCellRecording, data: pd.DataFrame) -> None:
        self.recording: WholeCellRecording = recording

        self.times: np.ndarray = data["times"].to_numpy(dtype=np.float64)

        Iinj_colnames: List[str] = list(data.columns)
        Iinj_colnames.remove("times")        
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

    def calculate_target_Qsyn(self) -> None:
        
        def cumtrapz_prefix_integral(arr: np.ndarray) -> np.ndarray:
            pref: np.ndarray = np.zeros((self.Nclamps, self.Nsamples + 1), dtype=arr.dtype) # shape: (Nclamps, Nsamples + 1)
            area: np.ndarray = 0.5 * (arr[:, 1:] + arr[:, :-1]) * self.recording.dt         # shape: (Nclamps, Nsamples)
            pref[:, 2:] = np.cumsum(area, axis=1)
            return pref
        
        l_bin_idxs: np.ndarray = self.binned_timeseries["left index"].to_numpy(np.int32)           # shape: (Nbins,)
        r_bin_idxs: np.ndarray = self.binned_timeseries["right index"].to_numpy(np.int32)           # shape: (Nbins,)
        bin_duration: np.ndarray = self.binned_timeseries["bin duration"].to_numpy(np.float64)  # units: Seconds, shape: (Nbins,)

        Vm_prefix_integral: np.ndarray = cumtrapz_prefix_integral(self.Vm)  # units: Webers, shape: (Nclamps, Nsamples + 1)
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
        integral_Vm = np.stack(self.binned_timeseries["integral Vm"].to_numpy()).astype(np.float64).T  # (Nclamps, Nbins)
        target_Qsyn = np.stack(self.binned_timeseries["target Qsyn"].to_numpy()).astype(np.float64).T  # (Nclamps, Nbins)

        Phi = integral_Vm
        Q   = target_Qsyn

        # Means per bin
        Phi_m = Phi.mean(axis=0)   # (Nbins,)
        Q_m   = Q.mean(axis=0)

        # Centered
        Phi_c = Phi - Phi_m[None, :]
        Q_c   = Q - Q_m[None, :]

        # Regression slope b and intercept a for each bin
        denom = np.sum(Phi_c * Phi_c, axis=0)         # var * (Nclamps-1) up to scale
        numer = np.sum(Phi_c * Q_c, axis=0)

        # Handle degenerate bins where Phi has no variation across clamps
        eps = np.finfo(np.float64).tiny
        b = np.where(np.abs(denom) > eps, numer / denom, np.nan)   # slope (Nbins,)
        a = Q_m - b * Phi_m                                        # intercept

        # Your derived params
        gsyn = -b                                                  # Siemens
        # Avoid divide-by-zero when gsyn ~ 0
        gsyn_safe = np.where(np.abs(gsyn) > eps, gsyn, np.nan)

        Eeff = a / (gsyn_safe * float(self.recording.bin_s))       # Volts

        # Diagnostics per bin
        # Your SSE formula: sum (Q - gsyn*(Eeff - Phi))^2
        # We can compute predicted Q directly from a + b*Phi (same fit)
        Q_hat = a[None, :] + b[None, :] * Phi
        resid = Q - Q_hat
        sse = np.sum(resid * resid, axis=0)

        # Proper per-bin R^2: 1 - SSE / SST, SST = sum (Q - mean(Q))^2 within the bin
        sst = np.sum((Q - Q_m[None, :]) ** 2, axis=0)
        r2 = np.where(sst > 0, 1.0 - (sse / sst), np.nan)

        C_grid = np.linspace(-0.130, 0.060, 4001)

        pdf, logpost = posterior_C_batch_grid(
            a[:, np.newaxis], gsyn_safe[:, np.newaxis], C_grid,
            sigma_A=5e-11, sigma_B=5e-10,
            mu_C=self.recording.Er, tau_C=0.05,
            mu_B=np.array([1e-9]*np.size(gsyn_safe)),  # per-problem B prior mean
            tau_B=1e-9
        )
        summ = summarize_posterior_batch(C_grid, pdf, cred_mass=0.95)

        # Store
        self.binned_timeseries["Eeff unbiased"] = Eeff
        self.binned_timeseries["Eeff bayesian"] =  summ["mean"]
        self.binned_timeseries["gsyn"] = gsyn
        self.binned_timeseries["SSE least-squares Qsyn"] = sse
        self.binned_timeseries["r2 least-squares Qsyn"] = r2

    def estimate_ge_gi(self) -> None:
        bin_duration: np.ndarray = self.binned_timeseries["bin duration"].to_numpy(dtype=np.float64)                # units: Seconds, shape: (Nbins,)
        integral_Vm: np.ndarray = np.stack(self.binned_timeseries["integral Vm"].to_numpy()).astype(np.float64).T   #type: ignore , units: Webers, shape: (Nclamps, Nbins)
        target_Qsyn: np.ndarray  = np.stack(self.binned_timeseries["target Qsyn"].to_numpy()).astype(np.float64).T  #type: ignore , units: Coulombs, shape: (Nclamps, Nbins)

        Ee: float = self.recording.Ee # units: Volts
        Ei: float = self.recording.Ei # units: Volts

        Phie: np.ndarray = Ee * bin_duration[np.newaxis, :] - integral_Vm # units: Webers, shape: (Nclamps, Nbins)
        Phii: np.ndarray = Ei * bin_duration[np.newaxis, :] - integral_Vm # units: Webers, shape: (Nclamps, Nbins)

        # dot products per bin (sum over clamps axis=0)
        eTe: np.ndarray = np.sum(Phie * Phie, axis=0)                   # units: Weber^2, shape: (Nbins,)
        iTi: np.ndarray = np.sum(Phii * Phii, axis=0)                   # units: Weber^2, shape: (Nbins,)
        eTi: np.ndarray = np.sum(Phie * Phii, axis=0)                   # units: Weber^2, shape: (Nbins,)
        eTq: np.ndarray = np.sum(Phie * target_Qsyn, axis=0)            # units: Joule * Second, shape: (Nbins,)
        iTq: np.ndarray = np.sum(Phii * target_Qsyn, axis=0)            # units: Joule * Second, shape: (Nbins,)
        qTq: np.ndarray = np.sum(target_Qsyn * target_Qsyn, axis=0)     # units: Coulomb^2, shape: (Nbins,)

        det: np.ndarray = eTe * iTi - eTi * eTi # units: Weber^4
        safe_det: np.ndarray = np.where(np.abs(det) > np.finfo(np.float64).tiny, det, np.nan)

        # unconstrained LS candidate
        ge_u: np.ndarray = ( iTi * eTq - eTi * iTq) / safe_det # units: Siemens, shape: (Nbins,)
        gi_u: np.ndarray = (-eTi * eTq + eTe * iTq) / safe_det # units: Siemens, shape: (Nbins,)

        r2_u = (
            qTq
            - 2.0 * ge_u * eTq
            - 2.0 * gi_u * iTq
            + (ge_u * ge_u) * eTe
            + 2.0 * ge_u * gi_u * eTi
            + (gi_u * gi_u) * iTi
        )

        # boundary candidates
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

        resnorm = np.sqrt(np.maximum(r2, 0.0))

        # cond(XtX) from eigenvalues of 2x2 gram matrix
        tr = eTe + iTi
        disc = np.sqrt((eTe - iTi) ** 2 + 4.0 * (eTi ** 2))
        lam1 = 0.5 * (tr + disc)
        lam2 = 0.5 * (tr - disc)
        cond = np.where(lam2 > 0, lam1 / lam2, np.inf)

        self.binned_timeseries["ge"] = ge
        self.binned_timeseries["gi"] = gi
        self.binned_timeseries["ge/gi residual norm"] = resnorm
        self.binned_timeseries["ge/gi matrix conditioning"] = cond

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

        ge_bins: np.ndarray = self.binned_timeseries["ge"].to_numpy(dtype=np.float64) # units: Siemens, shape: (Nbins,)
        gi_bins: np.ndarray = self.binned_timeseries["gi"].to_numpy(dtype=np.float64) # units: Siemens, shape: (Nbins,)

        Vpred: np.ndarray = np.empty_like(Vm, dtype=np.float64) # units: Volts, shape: (Nclamps, Nsamples)
        Vpred[:, 0] = Vm[:, 0]  # initial condition from measured Vm

        for r, l, ge, gi in zip(r_idx + 1, l_idx + 1, ge_bins, gi_bins):
            v0: np.ndarray = Vpred[:, l - 1].reshape(-1, 1)
            G: np.ndarray = ge + gi + gl
            vss: np.ndarray = (ge * Ee + gi * Ei + gl * Er + Iinj) / G
            t: np.ndarray = np.arange(1, r - l + 1) * dt
            Vpred[:, l:r] = (v0 - vss) * np.exp(-G * t / Cm) + vss

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
        nrows = 6
        fig, axs = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            sharex="all",
            sharey="row",
            figsize=(15, 10),
            constrained_layout=True,
        )
        fig.suptitle(filename.name)

        if ncols == 1:
            axs = np.expand_dims(axs, axis=1)

        for idx, paradigm in enumerate(recording.stimuli):
            stimulus = recording.stimuli[paradigm]

            times: np.ndarray = stimulus.times
            bin_times: np.ndarray = stimulus.binned_timeseries["left time"].to_numpy(np.float64)

            Vpred: np.ndarray = np.stack(stimulus.timeseries["predicted Vm"].to_numpy()).astype(np.float64).T

            Eeff: np.ndarray = stimulus.binned_timeseries["Eeff bayesian"].to_numpy(np.float64) 
            gsyn: np.ndarray = stimulus.binned_timeseries["gsyn"].to_numpy(np.float64)

            ge = stimulus.binned_timeseries["ge"].to_numpy(np.float64)
            gi = stimulus.binned_timeseries["gi"].to_numpy(np.float64)
            
            # diagnostics
            resnorm = stimulus.binned_timeseries["ge/gi residual norm"].to_numpy()
            cond = stimulus.binned_timeseries["ge/gi matrix conditioning"].to_numpy()

            axs[0, idx].set_title(paradigm)

            # ---- Row 0: Vm + Vpred (same colors) ----
            curves = axs[0, idx].plot(times, stimulus.Vm.T)  # solid Vm traces
            colors = [line.get_color() for line in curves]

            for j in range(stimulus.Nclamps):
                axs[0, idx].plot(times, Vpred[j, :], linestyle=":", color=colors[j])  # dotted Vpred

            axs[0, idx].set_ylabel("Vm (V)")
            axs[0, idx].grid(True)
            if idx == 0:
                axs[0, idx].legend(loc="upper right")

            # ---- Row 1: Eeff ----
            # axs[1, idx].plot(bin_times, Eeff, c="black")
            axs[1, idx].plot(bin_times, (1e9 * Eeff + stimulus.recording.Er / gsyn) / (1e9 + 1 / gsyn), c="black")
            axs[1, idx].axhline(recording.Er, linestyle="--", color="k", linewidth=1, label="Er" if idx == 0 else None)
            axs[1, idx].axhline(recording.Ee, linestyle="--", color="r", linewidth=1, label="Ee" if idx == 0 else None)
            axs[1, idx].axhline(recording.Ei, linestyle="--", color="b", linewidth=1, label="Ei" if idx == 0 else None)

            # ---- Row 3: conductances ----
            axs[3, idx].plot(bin_times, ge, c="r", label="ge")
            axs[3, idx].plot(bin_times, gi, c="b", label="gi")
            
            axs[3, idx].plot(bin_times, bin_times * 0, "--k", linewidth=1)
            axs[3, idx].set_ylabel("G (S)")
            axs[3, idx].grid(True)
            if idx == 0:
                axs[3, idx].legend(loc="upper right")

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
