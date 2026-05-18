import pandas as pd
import numpy as np
import os
import json
import hashlib
import re
from pathlib import Path
from dataclasses import dataclass, field, fields
from typing import List, Optional, Dict, Any


# ---------------------------------------------------------------------------
# Sub-configs
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PathsCfg:
    """
    Input directory structure (one .xlsx per "case", one case per subfolder):

        spreadsheets_input_dir/
            2005-11-01_1_baseline_prr_baseline_pulse_number/
                2005-11-01_1_baseline_prr_baseline_pulse_number.xlsx
                <raw .smr/.smrx files, .json spike times, .txt notes, etc.>
            2005-11-02_3_baseline_prr_baseline_pulse_number/
                2005-11-02_3_baseline_prr_baseline_pulse_number.xlsx
                ...
            ...

    Each case lives in its own subfolder. Each subfolder must contain exactly
    one .xlsx file (plus whatever auxiliary files are convenient to keep
    alongside it). The xlsx filename does not need to match the subfolder
    name -- discovery is by .xlsx extension within each subfolder.

    files_to_analyze, if set, restricts analysis to specific cases. Entries
    are SUBFOLDER NAMES (not filenames), so a case is identified by its
    directory rather than its xlsx file. Example:
        "files_to_analyze": ["2005-11-01_1_baseline_prr_baseline_pulse_number"]

    Subfolders with zero .xlsx files are skipped with a warning. Subfolders
    with more than one .xlsx file are an error -- the convention is one
    case per subfolder, and ambiguity should be resolved upstream rather
    than guessed at here.
    """
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
            # Now subfolder names rather than xlsx filenames.
            for case_name in files_to_analyze:
                case_dir = spreadsheets_input_dir / case_name
                assert case_dir.is_dir(), (
                    f"files_to_analyze entry '{case_name}' is not a "
                    f"subdirectory of {spreadsheets_input_dir}."
                )

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
        """
        Discover the .xlsx file for each case subfolder.

        For each subdirectory of spreadsheets_input_dir:
          - 0 xlsx files: skip with warning (subfolder has no analyzable
            spreadsheet -- maybe raw data only, or work in progress).
          - 1 xlsx file: include it.
          - 2+ xlsx files: raise (ambiguous; one case per subfolder).

        If files_to_analyze is set, only the listed subfolders are visited.
        """
        # Pick which subfolders to scan.
        if self.files_to_analyze is None:
            case_dirs: List[Path] = sorted(
                p for p in self.spreadsheets_input_dir.iterdir() if p.is_dir()
            )
        else:
            case_dirs = [self.spreadsheets_input_dir / name for name in self.files_to_analyze]

        paths_to_spreadsheets: List[Path] = []
        for case_dir in case_dirs:
            xlsx_files: List[Path] = [
                p for p in case_dir.iterdir() if p.suffix.lower() == ".xlsx"
            ]
            if len(xlsx_files) == 0:
                print(
                    f"  Warning: case subfolder '{case_dir.name}' contains no .xlsx file; "
                    f"skipping."
                )
                continue
            if len(xlsx_files) > 1:
                names = sorted(p.name for p in xlsx_files)
                raise AssertionError(
                    f"Case subfolder '{case_dir.name}' contains {len(xlsx_files)} "
                    f"xlsx files: {names}. Each case subfolder must contain "
                    f"exactly one xlsx."
                )
            paths_to_spreadsheets.append(xlsx_files[0])

        assert len(paths_to_spreadsheets) > 0, (
            f"No analyzable .xlsx files found under {self.spreadsheets_input_dir} "
            f"(either no subfolders, or every subfolder is empty of xlsx files)."
        )
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
    # Vss is currently estimated as the median of the first 5 ms of raw Vm
    # samples per (paradigm, Iinj). The fields below are vestigial from the
    # earlier weighted-mode-on-full-trace estimator; they are unused by the
    # current pipeline but kept around so old config files still load (and in
    # case we ever revert to weighted-mode for some traces). Pass any positive
    # values.
    Vss_bin_width:    float        # volts
    Vss_smooth_sigma: float        # volts
    Vss_Im_scale:     float        # amperes

    # Drift-cluster selection: pick the smallest K such that, for every cluster
    # at K, every (paradigm, clamp)'s Vss observation lies within
    # max_residual_volts of the cluster's OLS-fit line Vss_hat(Iinj) = Er + Rin
    # * Iinj. K=N (every paradigm its own cluster) is always feasible since a
    # 2-clamp paradigm's OLS line passes exactly through its two points, so the
    # search always terminates. This replaces the prior GMM/BIC-with-Er/Rin-
    # scaling clustering -- residuals are the directly meaningful quantity for
    # whether two paradigms share a leak state.
    max_residual_volts: float       # e.g. 2e-3 = 2 mV

    # Window at the start of each trace used to estimate Vss. Currently Vss =
    # median of raw Vm over the first Vss_duration_seconds. Set this to as long
    # as you're confident the cell is at equilibrium before the first stimulus
    # arrives. Longer windows give more samples to median over (lower noise);
    # too long and you start including stimulus-evoked deflections.
    Vss_duration_seconds: float     # e.g. 5e-3 = 5 ms

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ErRinCfg":
        _expect_keys(cls, d, "Er_Rin_estimation")
        assert d["Vss_bin_width"] > 0
        assert d["Vss_smooth_sigma"] > 0
        assert d["Vss_Im_scale"] > 0
        assert d["max_residual_volts"] > 0
        assert d["Vss_duration_seconds"] > 0
        return cls(
            Vss_bin_width=float(d["Vss_bin_width"]),
            Vss_smooth_sigma=float(d["Vss_smooth_sigma"]),
            Vss_Im_scale=float(d["Vss_Im_scale"]),
            max_residual_volts=float(d["max_residual_volts"]),
            Vss_duration_seconds=float(d["Vss_duration_seconds"]),
        )


@dataclass(frozen=True)
class EeEiCfg:
    """
    Configuration for the (Ee, Ei) estimator (quantile-of-Eeff method with
    physiological difference constraint).

    The estimator takes a low and high quantile of the pooled Eeff(t)
    distribution as initial estimates of Ei and Ee respectively. It then
    enforces a hard physiological constraint on the difference (Ee - Ei):

        dE_min ≤ Ee - Ei ≤ dE_max

    Rather than constraining Ee and Ei in absolute terms. Our recordings have
    multiple uncontrolled sources of voltage offset (most importantly an
    uncertain liquid junction potential; see project notes), so absolute
    reversal potentials are not recoverable per cell. The *difference*
    Ee - Ei is invariant under offset shifts and is the only well-determined
    quantity for the downstream conductance inversion.

    If the quantile estimate gives Ee - Ei outside [dE_min, dE_max], project
    to the nearest in-bound pair by adjusting both reversals symmetrically
    (preserving the midpoint, modifying only the spread).

    Defaults:
      low/high quantile: 0.001 / 0.999 (small percentage clipping for noise robustness)
      dE_min = 55 mV  (Ee = -10, Ei = -65 -- the tightest physiological separation)
      dE_max = 120 mV (Ee = +10, Ei = -110 -- the widest physiological separation)

    The Δ notation is retained throughout the pipeline. This is consistent
    with the assumption that decreases from baseline are artifactual (the
    operational stance for this paper, supported by prior duration-coder
    findings in the same brain region; see Rose et al. 2016) -- i.e., we
    believe baselines exist but assume disinhibition / disexcitation does not
    contribute meaningfully to the observed signals. Residual negative Δg
    after Iact correction is treated as model misspecification and analyzed
    for correlation with active membrane events.
    """
    low_quantile:  float
    high_quantile: float
    dE_min:        float
    dE_max:        float

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "EeEiCfg":
        _expect_keys(cls, d, "Ee_Ei_estimation")
        assert 0.0 <= d["low_quantile"] < d["high_quantile"] <= 1.0, (
            "must have 0 ≤ low_quantile < high_quantile ≤ 1"
        )
        assert 0.0 < d["dE_min"] < d["dE_max"], (
            "must have 0 < dE_min < dE_max"
        )
        return cls(
            low_quantile=float(d["low_quantile"]),
            high_quantile=float(d["high_quantile"]),
            dE_min=float(d["dE_min"]),
            dE_max=float(d["dE_max"]),
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


# ---------------------------------------------------------------------------
# Top-level config
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AnalyzerCfg:
    paths: PathsCfg
    compute: ComputeCfg
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
    """
    Reads .xlsx files in the lab's per-case format.

    File structure (one .xlsx per recording / "case"):

      - One sheet per stimulus paradigm (e.g. "40pps", "60pps", "1pulse").
        Layout within each paradigm sheet:
          Row 0:           "use data" label in column A; per-Iinj-column flags
                           in columns B onward. A flag is True (include this
                           clamp), False (exclude), or blank/NaN (include --
                           treated as missing = default True, which also covers
                           the case where the cell holds a =TRUE() formula whose
                           cached value isn't present yet).
          Row 1:           column headers. Column A is "times" (sample times in
                           seconds); columns B onward are Iinj values in
                           amperes, one per current-clamp condition.
          Row 2 onward:    data rows. Vm values in millivolts.

      - One "parameters" sheet with cell-level metadata. Currently:
          Cm (column header) and Cm value (single row in farads).
        Other parameters (Et, Eact, Ess) are no longer required -- they are
        either estimated by the pipeline or unused.
    """
    PARAMETERS_SHEET_NAME: str = "parameters"

    def __init__(self, filepath: Path):
        self.filepath: Path = filepath
        self.data_pointer = pd.ExcelFile(self.filepath, engine="openpyxl")

    @property
    def sheet_names(self) -> List[str]:
        assert isinstance(self.data_pointer.sheet_names, list)
        assert all(isinstance(x, str) for x in self.data_pointer.sheet_names)
        return self.data_pointer.sheet_names  # type: ignore

    def get_paradigms(self) -> List[str]:
        # Exclude metadata sheets. The "parameters" sheet name is reserved;
        # legacy files may also have "stats" / "results" sheets which we ignore.
        excluded = {"parameters", "stats", "results"}
        paradigms: List[str] = [s for s in self.sheet_names if s.lower() not in excluded]
        assert len(paradigms) > 0, f"No stimulus paradigm sheets found in {self.filepath}"
        return paradigms

    def get_paradigm_data(self, paradigm: str) -> pd.DataFrame:
        """
        Read a paradigm sheet and return a DataFrame with columns: "times" plus
        one column per included Iinj (named as the Iinj value as a float-string
        like "-3.000e-11"). Iinj columns flagged False in the "use data" row
        are dropped.

        Layout:
          - Row 0: "use data" label in col A; per-Iinj flags in cols B onward.
            Blank/NaN flag is treated as True (default include). The "use data"
            row only governs the Iinj columns; the "times" column is always kept.
          - Row 1: headers. Col A = "times", col B+ = Iinj values.
          - Row 2+: data.
        """
        # Read row 0 alone to get the flags. The label in A1 ("use data") is
        # ignored; flags start at B1 and correspond 1:1 to Iinj columns.
        flags_df = pd.read_excel(
            self.filepath, sheet_name=paradigm, header=None, nrows=1
        )
        # Interpret each flag cell: blank/NaN -> True (include by default);
        # otherwise cast to bool. This handles =TRUE()/=FALSE() formula cells
        # whose cached values aren't present (openpyxl returns NaN) by
        # defaulting to include, which is the safer of the two.
        flags_for_iinjs: List[bool] = []
        for v in flags_df.iloc[0, 1:].tolist():
            if v is None or (isinstance(v, float) and np.isnan(v)):
                flags_for_iinjs.append(True)
            else:
                flags_for_iinjs.append(bool(v))

        # Read the data with header=1 (row 1 as headers). Row 0 leaks in as
        # the first data row, but it's a row of NaNs under the now-numeric Iinj
        # columns, so dropna(how='all') clears it.
        data: pd.DataFrame = pd.read_excel(
            self.filepath, sheet_name=paradigm, header=1
        )
        data = data.dropna(how='all').reset_index(drop=True)

        assert "times" in data.columns, (
            f"No 'times' column found in sheet '{paradigm}' of {self.filepath.name}. "
            f"Columns: {list(data.columns)}"
        )

        # The Iinj columns are everything after "times".
        iinj_columns: List[str] = [c for c in data.columns if c != "times"]
        assert len(flags_for_iinjs) == len(iinj_columns), (
            f"'use data' row in sheet '{paradigm}' of {self.filepath.name} has "
            f"{len(flags_for_iinjs)} flag cells but the data has {len(iinj_columns)} "
            f"Iinj column(s). Check the spreadsheet layout."
        )
        cols_to_drop = [c for c, flag in zip(iinj_columns, flags_for_iinjs) if not flag]
        if cols_to_drop:
            data = data.drop(columns=cols_to_drop)

        # Rename Iinj columns to canonical float-string form for downstream code.
        for key in list(data.columns):
            if key not in {"times", "stimulus", "representative"}:
                data.rename(columns={key: f"{float(key):.3e}"}, inplace=True)
        return data

    def get_paradigm_parameters(self) -> pd.DataFrame:
        """
        Read the per-case parameters sheet. Currently expected to contain just
        one column ("Cm") and one row of values. Returns a DataFrame so the
        downstream code can keep its `parameters["Cm"][0]` access pattern.
        """
        df = pd.read_excel(
            self.filepath, sheet_name=self.PARAMETERS_SHEET_NAME, header=0
        )
        df = df.dropna(how='all').dropna(axis=1, how='all').reset_index(drop=True)
        assert "Cm" in df.columns, (
            f"No 'Cm' column in parameters sheet of {self.filepath.name}; "
            f"got columns: {list(df.columns)}"
        )
        assert len(df) >= 1, f"parameters sheet of {self.filepath.name} is empty"
        return df



# ---------------------------------------------------------------------------
# Per-case manual cluster + steady-state Vss system
# ---------------------------------------------------------------------------

# Iinj filename grammar:
#   "(<value><unit>)_<index>.txt"
#   <value>: signed decimal float (e.g. "-0.01", "0", "+1.5")
#   <unit>:  SI prefix from {f, p, n, u, m, ''} immediately followed by "A"
#            (empty prefix is bare amperes)
#   <index>: any nonnegative integer; not used semantically (just for uniqueness)
#
# Examples that match:
#   (-0.01nA)_1.txt    -> -0.01e-9 A = -10 pA
#   (0nA)_3.txt        ->  0 A
#   (-30pA)_2.txt      -> -30e-12 A
#   (1.5uA)_0.txt      ->  1.5e-6 A
_IINJ_FILENAME_RE = re.compile(
    r"^\(\s*"
    r"(?P<value>[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)"
    r"\s*(?P<unit>[fpnum]?)A\s*\)"
    r"_(?P<index>\d+)"
    r"\.txt$"
)

_SI_PREFIX_FACTORS: Dict[str, float] = {
    "":  1.0,
    "m": 1e-3,
    "u": 1e-6,
    "n": 1e-9,
    "p": 1e-12,
    "f": 1e-15,
}


def _parse_iinj_from_filename(filename: str) -> Optional[float]:
    """
    Parse an Iinj value (in amperes) from a steady-state .txt filename.
    Returns None if the name doesn't match the expected grammar.
    """
    m = _IINJ_FILENAME_RE.match(filename)
    if m is None:
        return None
    value: float = float(m.group("value"))
    unit:  str   = m.group("unit")
    return value * _SI_PREFIX_FACTORS[unit]


def _read_spike2_txt(path: Path) -> np.ndarray:
    """
    Read a Spike2 default-export .txt file. Format:
        "Time"   "Vm"
        0.00000000  -422.66846
        0.00010000  -420.37964
        ...
    Returns a 1D float array of Vm samples in VOLTS (input is mV, scaled by 1e-3).
    Time column is discarded; only the Vm values are used for steady-state mean.
    """
    # Use pandas: tolerant of whitespace, comments, header line.
    df = pd.read_csv(
        path,
        sep=r"\s+",
        engine="python",
        comment="#",
        skiprows=1,           # skip the "Time" "Vm" header
        header=None,
        names=["time", "vm_mV"],
    )
    vm_mV = df["vm_mV"].to_numpy(dtype=np.float64)
    return vm_mV * 1e-3       # mV -> V


@dataclass(frozen=True)
class ManualClusterDataCfg:
    """
    Per-case manual cluster definitions plus their steady-state Vss values,
    computed from user-provided Spike2 excerpts.

    Two pieces of input data, both in the case folder:

    1. `cluster_assignments.json`:
       {
         "stimulus_clusters": {
           "5pps":  "baseline",
           "10pps": "baseline",
           "40pps": "post-NBQX"
         },
         "cluster_labels": {                    // optional, decorative
           "baseline":  "Baseline (pre-drug)",
           "post-NBQX": "Post-NBQX wash"
         }
       }

    2. `steady_states/<cluster_name>/(<Iinj><unit>)_<idx>.txt`:
       Each .txt is a Spike2 default-export file (header line + columns of
       time, Vm in mV) for a stretch of confirmed steady-state recording at
       a known Iinj level. Multiple files per (cluster, Iinj) are pooled
       (samples concatenated) and reduced by mean to give Vss.

    The pipeline uses these directly to obtain (cluster_assignment, Vss),
    skipping the unsupervised weighted-mode + GMM/BIC. (Er, Rin) is then
    computed by OLS on (Iinj, Vss) pairs per cluster as usual.

    Both files (the JSON and the steady_states/ folder) must be present
    together; if only one exists, treat it as a user error and raise.
    If neither exists, the pipeline falls back to fully unsupervised mode.

    Validation (all enforced at load time):
      - cluster_assignments.json is well-formed; stimulus_clusters maps every
        xlsx paradigm to a cluster name; cluster names referenced in
        stimulus_clusters all exist as subfolders of steady_states/.
      - Each (cluster, Iinj) pair used by stimuli is covered by at least one
        .txt file.
      - Each .txt file's name parses as an Iinj value.
      - Each .txt file's content parses as a Spike2 export.
      - Warn if any .txt file's Iinj value isn't used by any stimulus in
        its cluster (extra data; may indicate mislabeling).
    """
    # cluster_name -> Iinj_amperes -> Vss_volts (mean across all .txt files)
    Vss_by_cluster:    Dict[str, Dict[float, float]]
    # paradigm_name -> cluster_name
    cluster_for_paradigm: Dict[str, str]
    # cluster_name -> human-readable label (or the cluster_name itself if no label)
    cluster_labels:    Dict[str, str]
    # Content hash for cache invalidation
    _content_hash:     str

    @classmethod
    def from_case_dir(
        cls,
        case_dir: Path,
        paradigm_iinjs: Dict[str, List[float]],
    ) -> Optional["ManualClusterDataCfg"]:
        """
        Attempt to load manual cluster data from `case_dir`.

        `paradigm_iinjs` is {paradigm_name: [Iinj amperes per clamp]} from the
        xlsx; used to validate that every (cluster, Iinj) pair the analysis
        actually needs is covered by the steady-state files.

        Returns None if neither the JSON nor the steady_states/ folder exists
        (i.e., the case is configured for fully-unsupervised analysis).
        Raises AssertionError on validation failures or if exactly one of the
        two inputs is present.
        """
        json_path:    Path = case_dir / "cluster_assignments.json"
        ss_dir:       Path = case_dir / "steady_states"

        json_exists = json_path.is_file()
        ss_exists   = ss_dir.is_dir()

        if not json_exists and not ss_exists:
            return None

        assert json_exists and ss_exists, (
            f"Inconsistent manual cluster data in {case_dir.name}: "
            f"found {'cluster_assignments.json' if json_exists else 'steady_states/'} "
            f"but not the other. Both must be present together, or neither."
        )

        # ---- Load JSON ----
        with open(json_path, "r") as f:
            jd = json.load(f)
        assert isinstance(jd, dict), f"{json_path}: top level must be an object"
        keys = {k for k in jd.keys() if not k.startswith("_")}
        assert keys.issubset({"stimulus_clusters", "cluster_labels"}), (
            f"{json_path}: unexpected keys {keys - {'stimulus_clusters', 'cluster_labels'}}. "
            f"Allowed: 'stimulus_clusters', 'cluster_labels'."
        )
        assert "stimulus_clusters" in keys, f"{json_path}: missing required 'stimulus_clusters'"
        stim_to_cluster_raw = jd["stimulus_clusters"]
        assert isinstance(stim_to_cluster_raw, dict), (
            f"{json_path}: 'stimulus_clusters' must be an object mapping paradigm -> cluster name"
        )
        for k, v in stim_to_cluster_raw.items():
            assert isinstance(k, str) and isinstance(v, str), (
                f"{json_path}: stimulus_clusters keys/values must be strings; got ({k!r}: {v!r})"
            )
        cluster_for_paradigm: Dict[str, str] = dict(stim_to_cluster_raw)

        labels_raw = jd.get("cluster_labels", {})
        assert isinstance(labels_raw, dict), f"{json_path}: 'cluster_labels' must be an object"
        for k, v in labels_raw.items():
            assert isinstance(k, str) and isinstance(v, str), (
                f"{json_path}: cluster_labels keys/values must be strings"
            )

        # Every paradigm from the xlsx must be assigned.
        missing = [p for p in paradigm_iinjs if p not in cluster_for_paradigm]
        assert not missing, (
            f"{json_path}: paradigms present in xlsx but missing from stimulus_clusters: {missing}"
        )
        # No extra paradigms in the JSON that aren't in the xlsx.
        extra = [p for p in cluster_for_paradigm if p not in paradigm_iinjs]
        assert not extra, (
            f"{json_path}: stimulus_clusters references paradigms not present in xlsx: {extra}"
        )

        cluster_names: set = set(cluster_for_paradigm.values())

        # cluster_labels keys must be a subset of actually-used cluster names.
        for cname in labels_raw:
            assert cname in cluster_names, (
                f"{json_path}: cluster_labels has '{cname}' which isn't used in stimulus_clusters"
            )
        cluster_labels: Dict[str, str] = {c: labels_raw.get(c, c) for c in cluster_names}

        # ---- Scan steady_states/ ----
        # Each cluster name must be a subfolder; each subfolder contains .txt files.
        Vss_by_cluster: Dict[str, Dict[float, float]] = {}
        # Track files actually used vs files seen, for the "extra data" warning.
        files_seen: List[Path] = []
        files_used: set = set()

        for cluster_name in sorted(cluster_names):
            cluster_dir = ss_dir / cluster_name
            assert cluster_dir.is_dir(), (
                f"steady_states/{cluster_name} not found in {case_dir.name}. "
                f"Every cluster named in stimulus_clusters must have a subfolder."
            )

            # Group all .txt files by parsed Iinj.
            samples_by_iinj: Dict[float, List[np.ndarray]] = {}
            for p in sorted(cluster_dir.iterdir()):
                if p.suffix.lower() != ".txt":
                    continue
                files_seen.append(p)
                iinj = _parse_iinj_from_filename(p.name)
                assert iinj is not None, (
                    f"{p}: filename does not parse as '(<value><unit>A)_<index>.txt'. "
                    f"Expected forms like '(-0.01nA)_1.txt', '(0pA)_3.txt'."
                )
                # Read and store; we'll average later.
                vm = _read_spike2_txt(p)
                samples_by_iinj.setdefault(iinj, []).append(vm)

            # Compute Vss = mean across all pooled samples per Iinj.
            Vss_by_cluster[cluster_name] = {}
            for iinj, arrays in samples_by_iinj.items():
                pooled = np.concatenate(arrays)
                Vss_by_cluster[cluster_name][iinj] = float(np.mean(pooled))

        # ---- Cross-validate: every (cluster, Iinj) used by stimuli is covered ----
        # And warn if any .txt file's Iinj isn't used by any stimulus in its cluster.
        IINJ_MATCH_TOL: float = 1e-15
        def _find_iinj_match(target: float, available: List[float]) -> Optional[float]:
            for a in available:
                if abs(a - target) <= IINJ_MATCH_TOL:
                    return a
            return None

        for paradigm, iinjs in paradigm_iinjs.items():
            cluster_name = cluster_for_paradigm[paradigm]
            available = list(Vss_by_cluster[cluster_name].keys())
            for iinj in iinjs:
                match = _find_iinj_match(iinj, available)
                assert match is not None, (
                    f"steady_states/{cluster_name}: no .txt file covers Iinj = "
                    f"{iinj*1e12:.1f} pA (needed by paradigm '{paradigm}'). "
                    f"Available Iinj in this cluster: "
                    f"{[f'{a*1e12:.1f} pA' for a in sorted(available)]}."
                )

        # Identify Iinj values present in .txt files but not used by any stimulus.
        for cluster_name, vss_dict in Vss_by_cluster.items():
            stims_in_cluster = [p for p, c in cluster_for_paradigm.items() if c == cluster_name]
            needed_iinjs: set = set()
            for p in stims_in_cluster:
                for iinj in paradigm_iinjs[p]:
                    needed_iinjs.add(round(iinj, 18))   # canonical key for set membership
            for iinj_avail in vss_dict:
                if round(iinj_avail, 18) not in needed_iinjs:
                    # Approximate match check, since floating point.
                    found = any(abs(iinj_avail - n) <= IINJ_MATCH_TOL for n in needed_iinjs)
                    if not found:
                        print(
                            f"  Warning: steady_states/{cluster_name}/ has data for "
                            f"Iinj = {iinj_avail*1e12:.1f} pA which isn't used by any "
                            f"stimulus in cluster '{cluster_name}'. Check for mislabeling."
                        )

        # ---- Content hash for cache invalidation ----
        # Hash the canonical JSON content AND each .txt file's mtime + first few bytes
        # (full file mtime + size is enough to detect changes without reading everything).
        hasher = hashlib.sha256()
        canonical = json.dumps({
            "stimulus_clusters": cluster_for_paradigm,
            "cluster_labels":    {c: cluster_labels[c] for c in sorted(cluster_labels)},
        }, sort_keys=True).encode("utf-8")
        hasher.update(canonical)
        for f in sorted(files_seen):
            stat = f.stat()
            hasher.update(f.relative_to(case_dir).as_posix().encode("utf-8"))
            hasher.update(str(stat.st_size).encode("utf-8"))
            hasher.update(str(stat.st_mtime).encode("utf-8"))
        content_hash = hasher.hexdigest()[:16]

        return cls(
            Vss_by_cluster=Vss_by_cluster,
            cluster_for_paradigm=cluster_for_paradigm,
            cluster_labels=cluster_labels,
            _content_hash=content_hash,
        )

    def hash_for_cache(self) -> str:
        return self._content_hash

    # Convenience for the pipeline:
    def lookup_vss(self, paradigm: str, iinj: float) -> float:
        """
        Return the manually-defined Vss for (paradigm, iinj) in volts.
        Tolerates small floating-point mismatch in the Iinj value.
        """
        cluster_name = self.cluster_for_paradigm[paradigm]
        vss_dict = self.Vss_by_cluster[cluster_name]
        for avail_iinj, vss in vss_dict.items():
            if abs(avail_iinj - iinj) <= 1e-15:
                return vss
        raise KeyError(
            f"Manual Vss not available for paradigm '{paradigm}' (cluster "
            f"'{cluster_name}') at Iinj = {iinj*1e12:.1f} pA. "
            f"Available: {sorted(vss_dict.keys())}"
        )