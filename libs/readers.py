import pandas as pd
import os
import json
from pathlib import Path
from dataclasses import dataclass, fields
from typing import List, Optional, Dict, Any

@dataclass(frozen=True)
class AnalyzerCfg:
    spreadsheets_input_dir: Path 
    image_save_dir: Path
    image_save_type: str
    files_to_analyze: Optional[List[str]]=None
    iinj_clamps_to_use: Optional[List[List[float]]]=None

    @classmethod
    def from_json(cls, path_to_json: Path) -> 'AnalyzerCfg':
        with open(path_to_json, "r") as f:
            kwargs: Dict[str, Any] = json.load(f)

        assert set(x.name for x in fields(cls)) == set(kwargs.keys()), "Mismatch between analyzercfg kwargs and those provided in json settings."
        
        spreadsheets_input_dir: Path = Path(kwargs["spreadsheets_input_dir"])
        image_save_dir: Path = Path(kwargs["image_save_dir"])

        if not os.path.isdir(spreadsheets_input_dir):
            spreadsheets_input_dir = spreadsheets_input_dir.expanduser()
        if not os.path.isdir(image_save_dir):
            image_save_dir = image_save_dir.expanduser()

        assert os.path.isdir(spreadsheets_input_dir), spreadsheets_input_dir
        assert os.path.isdir(image_save_dir), image_save_dir
        assert kwargs["image_save_type"] in {'png', 'emf'}
        if kwargs["files_to_analyze"] is not None:
            assert isinstance(kwargs["files_to_analyze"], list)
            assert all(isinstance(x, str) for x in kwargs["files_to_analyze"])
            assert len(kwargs["files_to_analyze"]) > 0
            assert all(os.path.splitext(x)[1] == ".xlsx" for x in kwargs["files_to_analyze"])
            assert all((spreadsheets_input_dir / x).is_file() for x in kwargs["files_to_analyze"])
        if kwargs["iinj_clamps_to_use"] is not None:
            assert isinstance(kwargs["iinj_clamps_to_use"], list)
            assert all(isinstance(x, list) for x in kwargs["iinj_clamps_to_use"])
            assert all(len(x) >= 2 for x in kwargs["iinj_clamps_to_use"])
            assert all(isinstance(y, float) for x in kwargs["iinj_clamps_to_use"] for y in x)

        return AnalyzerCfg(spreadsheets_input_dir, image_save_dir, kwargs["image_save_type"], kwargs["files_to_analyze"], kwargs["iinj_clamps_to_use"])

    @property
    def paths_to_spreadsheets(self) -> List[Path]:
        paths_to_spreadsheets: List[Path]
        if self.files_to_analyze is None:
            paths_to_spreadsheets = [x for x in self.spreadsheets_input_dir.iterdir() if x.suffix == ".xlsx"]
        else:
            paths_to_spreadsheets = [self.spreadsheets_input_dir / x for x in self.files_to_analyze]
        
        assert all([os.path.exists(i) for i in paths_to_spreadsheets])
        assert len(paths_to_spreadsheets) > 0
        
        return paths_to_spreadsheets


class XLReader:
    def __init__(self, filepath: Path):
        self.filepath: Path = filepath
        self.data_pointer = pd.ExcelFile(self.filepath, engine="openpyxl")

    @property
    def sheet_names(self) -> List[str]:
        assert isinstance(self.data_pointer.sheet_names, list)
        assert all(isinstance(x, str) for x in self.data_pointer.sheet_names)
        return self.data_pointer.sheet_names #type: ignore (type checker can't infer that this is List[str], but above assertions do guarantee this)

    def get_paradigms(self) -> List[str]:

        paradigms: List[str] = list(filter(lambda x: not any(y in x.lower() for y in ["parameters", "stats", "results"]), self.sheet_names))
        assert len(paradigms) > 0
        return paradigms

    def get_paradigm_data(self, paradigm: str) -> pd.DataFrame:
        data: pd.DataFrame = pd.read_excel(self.filepath, sheet_name=paradigm, header=0)
        assert "times" in data, f"No 'times' column present in {paradigm} sheet." 
        for key in data.keys():
            if key not in {"times", "stimulus", "representative"}: # stimulus is often (optionally) stored alongside membrane potential averages
                data.rename(columns={key: f"{float(key):.3e}"}, inplace=True)                
        return data
    
    def get_paradigm_parameters(self, paradigm: str):
        paradigm = f"parameters_{paradigm}"
        df = pd.read_excel(self.filepath, sheet_name=paradigm, header=0)
        df = df.dropna(how='all').dropna(axis=1, how='all')
        df = df.dropna(how='any')
        assert all(df["Eact"] > df["Ess"]), f"Fix spreadsheet {self.filepath} {paradigm}: Eact must be greater than Ess"
        return df
    