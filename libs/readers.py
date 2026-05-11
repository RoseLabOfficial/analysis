import pandas as pd
import os
import json
import hashlib
from pathlib import Path
from dataclasses import dataclass, field, fields
from typing import List, Optional, Dict, Any


# ---------------------------------------------------------------------------
# Sub-configs
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PathsCfg:
    spreadsheets_input_dir: Path
    image_save_dir: Path
    image_save_type: str
    files_to_analyze: Optional[List[str]] = None
    iinj_clamps_to_use: Optional[List[List[float]]] = None
    # Parameter cache: per-spreadsheet JSON of computed (Er, Rin, Vss, Et,
    # x_beta, ...) so re-runs can skip the expensive optimization. Set
    # cache_dir to null to disable. cache_invalidate=true forces recompute
    # and overwrites the cache. The cache lives at <cache_dir>/<spreadsheet
    # stem>.json by default.
    cache_dir: Optional[Path] = None
    cache_invalidate: bool = False

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "PathsCfg":
        _expect_keys(cls, d, "paths")

        spreadsheets_input_dir: Path = Path(d["spreadsheets_input_dir"])
        image_save_dir: Path = Path(d["image_save_dir"])

        if not os.path.isdir(spreadsheets_input_dir):
            spreadsheets_input_dir = spreadsheets_input_dir.expanduser()
        if not os.path.isdir(image_save_dir):
            image_save_dir = image_save_dir.expanduser()

        assert os.path.isdir(spreadsheets_input_dir), spreadsheets_input_dir
        assert os.path.isdir(image_save_dir), image_save_dir
        assert d["image_save_type"] in {"png", "emf"}

        files_to_analyze = d["files_to_analyze"]
        if files_to_analyze is not None:
            assert isinstance(files_to_analyze, list)
            assert all(isinstance(x, str) for x in files_to_analyze)
            assert len(files_to_analyze) > 0
            assert all(os.path.splitext(x)[1] == ".xlsx" for x in files_to_analyze)
            assert all((spreadsheets_input_dir / x).is_file() for x in files_to_analyze)

        iinj_clamps_to_use = d["iinj_clamps_to_use"]
        if iinj_clamps_to_use is not None:
            assert isinstance(iinj_clamps_to_use, list)
            assert all(isinstance(x, list) for x in iinj_clamps_to_use)
            assert all(len(x) >= 2 for x in iinj_clamps_to_use)
            assert all(isinstance(y, float) for x in iinj_clamps_to_use for y in x)

        cache_dir_raw = d["cache_dir"]
        if cache_dir_raw is None:
            cache_dir: Optional[Path] = None
        else:
            assert isinstance(cache_dir_raw, str)
            cache_dir = Path(cache_dir_raw)
            if not cache_dir.exists():
                cache_dir = cache_dir.expanduser()
            cache_dir.mkdir(parents=True, exist_ok=True)
            assert cache_dir.is_dir()

        cache_invalidate = d["cache_invalidate"]
        assert isinstance(cache_invalidate, bool)

        return cls(
            spreadsheets_input_dir=spreadsheets_input_dir,
            image_save_dir=image_save_dir,
            image_save_type=d["image_save_type"],
            files_to_analyze=files_to_analyze,
            iinj_clamps_to_use=iinj_clamps_to_use,
            cache_dir=cache_dir,
            cache_invalidate=cache_invalidate,
        )

    @property
    def paths_to_spreadsheets(self) -> List[Path]:
        if self.files_to_analyze is None:
            paths_to_spreadsheets = [x for x in self.spreadsheets_input_dir.iterdir() if x.suffix == ".xlsx"]
        else:
            paths_to_spreadsheets = [self.spreadsheets_input_dir / x for x in self.files_to_analyze]

        assert all(os.path.exists(i) for i in paths_to_spreadsheets)
        assert len(paths_to_spreadsheets) > 0

        return paths_to_spreadsheets


@dataclass(frozen=True)
class ComputeCfg:
    n_workers: int

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ComputeCfg":
        _expect_keys(cls, d, "compute")
        assert isinstance(d["n_workers"], int) and d["n_workers"] >= 1
        return cls(n_workers=d["n_workers"])


@dataclass(frozen=True)
class FilterCfg:
    passband_hz: float
    stopband_hz: float
    attenuation_db: float
    ripple_db: float

    @classmethod
    def from_dict(cls, d: Dict[str, Any], name: str) -> "FilterCfg":
        _expect_keys(cls, d, f"filters.{name}")
        assert d["passband_hz"] < d["stopband_hz"], (
            f"filters.{name}: passband_hz must be < stopband_hz "
            f"(got {d['passband_hz']} >= {d['stopband_hz']})"
        )
        for key in ("passband_hz", "stopband_hz", "attenuation_db", "ripple_db"):
            assert d[key] > 0, f"filters.{name}.{key} must be positive"
        return cls(
            passband_hz=float(d["passband_hz"]),
            stopband_hz=float(d["stopband_hz"]),
            attenuation_db=float(d["attenuation_db"]),
            ripple_db=float(d["ripple_db"]),
        )


@dataclass(frozen=True)
class FiltersCfg:
    # Note: dg, dgsyn are filters for the *signed* deviation-from-baseline
    # quantities (Δge, Δgi, Δgsyn). The estimator returns Δg, not absolute g;
    # negative values within [-g0, 0] are biophysically valid (disinhibition,
    # withdrawal of tonic excitation) and should not be assumed to be artifacts.
    Vm: FilterCfg
    Im: FilterCfg
    dg: FilterCfg
    Eeff: FilterCfg
    dgsyn: FilterCfg

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "FiltersCfg":
        _expect_keys(cls, d, "filters")
        return cls(
            Vm=FilterCfg.from_dict(d["Vm"], "Vm"),
            Im=FilterCfg.from_dict(d["Im"], "Im"),
            dg=FilterCfg.from_dict(d["dg"], "dg"),
            Eeff=FilterCfg.from_dict(d["Eeff"], "Eeff"),
            dgsyn=FilterCfg.from_dict(d["dgsyn"], "dgsyn"),
        )


@dataclass(frozen=True)
class ErRinCfg:
    # Vss histogram-mode estimation parameters.
    # Bin width and smoothing sigma are in volts (rather than counts of bins).
    # Volts-based makes the estimator behaviour independent of each clamp's
    # actual Vm range -- a depolarizing clamp with a narrow distribution and a
    # hyperpolarizing clamp with a wide one will get the same effective resolution.
    Vss_bin_width:    float        # volts; e.g. 5e-4 = 0.5 mV
    Vss_smooth_sigma: float        # volts; e.g. 1e-3 = 1.0 mV
    # Weighting of the histogram by exp(-|Im_filtered| / Vss_Im_scale) so quiet
    # baseline samples (where |Im| ~ 0) contribute fully and PSP-transient
    # samples (where |Im| is large) contribute less. At |Im| = Vss_Im_scale,
    # the weight is exp(-1) ~ 0.37; at |Im| = 3*Vss_Im_scale, weight ~ 0.05.
    # Choose Vss_Im_scale around the |Im| magnitude of a typical PSP transient
    # in your recordings (e.g. 50-100 pA). Set to +inf to disable weighting.
    Vss_Im_scale: float            # amperes

    # BIC penalty multiplier for drift-cluster selection. Standard BIC has alpha=1;
    # alpha>1 makes adding clusters harder, biasing toward fewer drift states.
    cluster_penalty_alpha: float
    # Physical scales for Er (volts) and Rin (ohms). Per-stimulus (Er, Rin)
    # estimates are normalized by these before clustering. Pick to roughly match
    # the smallest drift you'd want to detect (e.g. 5 mV in Er, 50 MOhm in Rin).
    cluster_scale_Er:  float        # volts
    cluster_scale_Rin: float        # ohms
    # Cluster variance floor in normalized units (square of per-stim-noise /
    # cluster_scale). Without this, GMM drives variance to zero and log-
    # likelihood to +inf, which lets it create spurious clusters indistinguishable
    # from noise. E.g., if your typical Er-estimate noise is ~0.5 mV and
    # cluster_scale_Er is 5 mV, set this to (0.5/5)^2 = 0.01.
    cluster_noise_floor: float

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ErRinCfg":
        _expect_keys(cls, d, "Er_Rin_estimation")
        assert d["Vss_bin_width"] > 0
        assert d["Vss_smooth_sigma"] > 0
        assert d["Vss_Im_scale"] > 0
        assert d["cluster_penalty_alpha"] > 0
        assert d["cluster_scale_Er"]  > 0
        assert d["cluster_scale_Rin"] > 0
        assert d["cluster_noise_floor"] > 0
        return cls(
            Vss_bin_width=float(d["Vss_bin_width"]),
            Vss_smooth_sigma=float(d["Vss_smooth_sigma"]),
            Vss_Im_scale=float(d["Vss_Im_scale"]),
            cluster_penalty_alpha=float(d["cluster_penalty_alpha"]),
            cluster_scale_Er=float(d["cluster_scale_Er"]),
            cluster_scale_Rin=float(d["cluster_scale_Rin"]),
            cluster_noise_floor=float(d["cluster_noise_floor"]),
        )


@dataclass(frozen=True)
class EeEiCfg:
    """
    Hard physiological priors for the (Ee, Ei) estimator.

    The estimator solves for (Ee, Ei) in [Ee_min, Ee_max] x [Ei_min, Ei_max]
    that minimizes the integrated bound-violation loss
        L = ∫ [max(0, lower_bound - Δge(t))² + max(0, lower_bound - Δgi(t))²] dt
    where lower_bound = -(g_l'_j - g_l_min) for cluster j. This penalizes
    inferred Δge or Δgi that fall below the physically allowed floor (set by
    the cluster's effective leak gl' minus a global pure-leak floor gl_min).

    When the loss has a flat minimum region (data underdetermines the choice
    -- e.g., the box's interior contains many (Ee, Ei) that produce no
    bound violations), tiebreak by Euclidean distance to (Ee_prior_center,
    Ei_prior_center). The data fundamentally cannot identify (Ee, Ei) per
    cell from a single recording; the prior center is an honest fallback.

    Physiological priors (post-LJP):
      Ee:  AMPA/NMDA-driven, true reversal in [-10, +10] mV.
      Ei:  Lumped Cl- (around -75 mV) and K+ (around -100 mV) reversals;
           "effective" Ei drifts toward whichever inhibitory conductance
           dominates. Box [-110, -65] mV captures both extremes.
      Center (0, -75) mV is a sensible mid-physiological default.
      gl_min: a global lower bound on the cell's pure (non-synaptic) leak
           conductance. Used to set the per-cluster Δg lower bound:
           Δg ≥ -(gl'_j - gl_min). Smaller gl_min → looser Δg bound (more
           generous to disinhibition); larger gl_min → tighter bound. For
           IC neurons in Rana pipiens with typical Rin ~400 MΩ (gl' ~ 2.5 nS),
           gl_min = 0.5 nS is a defensible floor based on K+ leak
           contributions in central neurons.
    """
    Ee_min:            float
    Ee_max:            float
    Ei_min:            float
    Ei_max:            float
    Ee_prior_center:   float
    Ei_prior_center:   float
    gl_min:            float

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "EeEiCfg":
        _expect_keys(cls, d, "Ee_Ei_estimation")
        assert d["Ee_min"] < d["Ee_max"], "Ee_min must be < Ee_max"
        assert d["Ei_min"] < d["Ei_max"], "Ei_min must be < Ei_max"
        assert d["Ee_min"] <= d["Ee_prior_center"] <= d["Ee_max"], (
            "Ee_prior_center must lie within [Ee_min, Ee_max]"
        )
        assert d["Ei_min"] <= d["Ei_prior_center"] <= d["Ei_max"], (
            "Ei_prior_center must lie within [Ei_min, Ei_max]"
        )
        assert d["gl_min"] > 0, "gl_min must be positive (units: Siemens)"
        return cls(
            Ee_min=float(d["Ee_min"]),
            Ee_max=float(d["Ee_max"]),
            Ei_min=float(d["Ei_min"]),
            Ei_max=float(d["Ei_max"]),
            Ee_prior_center=float(d["Ee_prior_center"]),
            Ei_prior_center=float(d["Ei_prior_center"]),
            gl_min=float(d["gl_min"]),
        )


@dataclass(frozen=True)
class EtOptCfg:
    x_beta_max: float
    eps_v: float
    n_epochs: int
    n_samples_per_epoch: int
    shrink: float
    rng_seed: int

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "EtOptCfg":
        _expect_keys(cls, d, "Et_optimization")
        assert d["x_beta_max"] > 0
        assert d["eps_v"] > 0
        assert isinstance(d["n_epochs"], int) and d["n_epochs"] >= 1
        assert isinstance(d["n_samples_per_epoch"], int) and d["n_samples_per_epoch"] >= 1
        assert 0.0 < d["shrink"] < 1.0
        assert isinstance(d["rng_seed"], int)
        return cls(
            x_beta_max=float(d["x_beta_max"]),
            eps_v=float(d["eps_v"]),
            n_epochs=int(d["n_epochs"]),
            n_samples_per_epoch=int(d["n_samples_per_epoch"]),
            shrink=float(d["shrink"]),
            rng_seed=int(d["rng_seed"]),
        )


@dataclass(frozen=True)
class NumericsCfg:
    # dgsyn_zero_threshold: absolute-value floor for division by Δgsyn when
    # computing Eeff = -a / Δgsyn. Purely numerical (prevents blow-up for
    # near-zero |Δgsyn|); has no semantic constraint on the sign of Δgsyn.
    dgsyn_zero_threshold: float
    box_volume_eps: float
    box_volume_floor: float

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "NumericsCfg":
        _expect_keys(cls, d, "numerics")
        for key in ("dgsyn_zero_threshold", "box_volume_eps", "box_volume_floor"):
            assert d[key] > 0, f"numerics.{key} must be positive"
        return cls(
            dgsyn_zero_threshold=float(d["dgsyn_zero_threshold"]),
            box_volume_eps=float(d["box_volume_eps"]),
            box_volume_floor=float(d["box_volume_floor"]),
        )


@dataclass(frozen=True)
class RecordingCfg:
    """
    Recording-prep parameters that adjust raw measurements.

    Currently a single global field for the liquid junction potential. This is
    a temporary simplification: in practice LJP differs per pipette solution
    (e.g., K-gluconate vs. KF), and ideally would be specified per-cell in
    spreadsheet metadata. For now we apply the same value to all cells; the
    user is responsible for using a cfg with the matching LJP for the pipette
    solution used in the recordings being analyzed.

    LJP convention (from LJPcalc / Marino et al. 2014, Barry & Lynch):
    LJP is reported as the bath potential relative to the pipette. The
    correction subtracts this value from amplifier readings:
        V_true = V_measured - V_LJP
    A cell measured at -70 mV with V_LJP = +17 mV is actually at -87 mV.
    """
    # Liquid junction potential in volts (e.g., +17.21e-3 = +17.21 mV).
    # Applied to all measured Vm and to the user-supplied Et_measured at ingestion.
    # Set to 0.0 to disable correction.
    liquid_junction_potential_volts: float

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RecordingCfg":
        _expect_keys(cls, d, "recording")
        return cls(
            liquid_junction_potential_volts=float(d["liquid_junction_potential_volts"]),
        )


# ---------------------------------------------------------------------------
# Top-level config
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AnalyzerCfg:
    paths: PathsCfg
    compute: ComputeCfg
    recording: RecordingCfg
    filters: FiltersCfg
    Er_Rin_estimation: ErRinCfg
    Ee_Ei_estimation: EeEiCfg
    Et_optimization: EtOptCfg
    numerics: NumericsCfg

    @classmethod
    def from_json(cls, path_to_json: Path) -> "AnalyzerCfg":
        with open(path_to_json, "r") as f:
            kwargs: Dict[str, Any] = json.load(f)

        _expect_keys(cls, kwargs, "<root>")

        return cls(
            paths=PathsCfg.from_dict(kwargs["paths"]),
            compute=ComputeCfg.from_dict(kwargs["compute"]),
            recording=RecordingCfg.from_dict(kwargs["recording"]),
            filters=FiltersCfg.from_dict(kwargs["filters"]),
            Er_Rin_estimation=ErRinCfg.from_dict(kwargs["Er_Rin_estimation"]),
            Ee_Ei_estimation=EeEiCfg.from_dict(kwargs["Ee_Ei_estimation"]),
            Et_optimization=EtOptCfg.from_dict(kwargs["Et_optimization"]),
            numerics=NumericsCfg.from_dict(kwargs["numerics"]),
        )

    # Convenience pass-through (kept for back-compat with code that read this directly off the cfg).
    @property
    def paths_to_spreadsheets(self) -> List[Path]:
        return self.paths.paths_to_spreadsheets

    def hash_for_cache(self) -> str:
        """
        Hash of cfg fields that affect computed parameters. Excludes:
        - paths (output-only or runtime configuration)
        - compute (n_workers; performance only)
        Used to invalidate cached parameters when analysis settings change.
        """
        from dataclasses import asdict
        relevant = {
            "recording":         asdict(self.recording),
            "filters":           asdict(self.filters),
            "Er_Rin_estimation": asdict(self.Er_Rin_estimation),
            "Ee_Ei_estimation":  asdict(self.Ee_Ei_estimation),
            "Et_optimization":   asdict(self.Et_optimization),
            "numerics":          asdict(self.numerics),
        }
        # JSON dump with sorted keys gives a deterministic representation.
        payload = json.dumps(relevant, sort_keys=True).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()[:16]   # 16 hex chars is plenty for cache keying


def _expect_keys(cls, d: Dict[str, Any], section_name: str) -> None:
    """Verify dict has exactly the dataclass field names; helpful error otherwise.

    Keys starting with '_' are ignored, allowing user comments in JSON
    (e.g., "_comment": "explanation of the values below").
    """
    expected = set(f.name for f in fields(cls))
    actual = set(k for k in d.keys() if not k.startswith("_"))
    missing = expected - actual
    extra = actual - expected
    assert not missing and not extra, (
        f"settings section [{section_name}]: "
        f"{'missing keys: ' + str(sorted(missing)) + '. ' if missing else ''}"
        f"{'unexpected keys: ' + str(sorted(extra)) + '. ' if extra else ''}"
    )


# ---------------------------------------------------------------------------
# Spreadsheet reader (unchanged)
# ---------------------------------------------------------------------------

class XLReader:
    def __init__(self, filepath: Path):
        self.filepath: Path = filepath
        self.data_pointer = pd.ExcelFile(self.filepath, engine="openpyxl")

    @property
    def sheet_names(self) -> List[str]:
        assert isinstance(self.data_pointer.sheet_names, list)
        assert all(isinstance(x, str) for x in self.data_pointer.sheet_names)
        return self.data_pointer.sheet_names  # type: ignore

    def get_paradigms(self) -> List[str]:
        paradigms: List[str] = list(filter(lambda x: not any(y in x.lower() for y in ["parameters", "stats", "results"]), self.sheet_names))
        assert len(paradigms) > 0
        return paradigms

    def get_paradigm_data(self, paradigm: str) -> pd.DataFrame:
        data: pd.DataFrame = pd.read_excel(self.filepath, sheet_name=paradigm, header=0)
        assert "times" in data, f"No 'times' column present in {paradigm} sheet."
        for key in data.keys():
            if key not in {"times", "stimulus", "representative"}:
                data.rename(columns={key: f"{float(key):.3e}"}, inplace=True)
        return data

    def get_paradigm_parameters(self, paradigm: str):
        paradigm = f"parameters_{paradigm}"
        df = pd.read_excel(self.filepath, sheet_name=paradigm, header=0)
        df = df.dropna(how='all').dropna(axis=1, how='all')
        df = df.dropna(how='any')
        assert all(df["Eact"] > df["Ess"]), f"Fix spreadsheet {self.filepath} {paradigm}: Eact must be greater than Ess"
        return df