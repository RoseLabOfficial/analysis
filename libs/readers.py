import pandas as pd
import os
from dataclasses import dataclass
from typing import List, Optional, Dict

@dataclass(frozen=True)
class FilterCfg:
    passband: float 
    ripple: float 
    stopband: float 
    attenuation: float

@dataclass(frozen=True)
class AnalyzerCfg:
    spreadsheets_input_dir: str 
    image_save_dir: str 
    image_save_type: Optional[str]
    timeseries_filter_cfgs: Dict[str, FilterCfg]
    files_to_analyze: str | List[str]='ALL' # ALL runs all files in input directory. Single str runs just that file. List[str] runs all specified files.
    optimization_level: int=0 # 0: No optimization 1: Calculate Eact based on max depolarization paradigm 2: Calculate Eact for each paradigm.
    iinj_clamps_to_use: Optional[List[float] | List[List[float]]]=None # If None, uses all current clamps. If List[float], uses specified clamps

    @classmethod
    def new(
        cls, 
        spreadsheets_input_dir: Optional[str]=None, 
        image_save_dir: Optional[str]=None, 
        iamge_save_type: Optional[str]=None,
        timeseries_filter_cfgs: Optional[Dict[str, Dict[str, float]]]=None, 
        files_to_analyze: Optional[str | List[str]]=None, 
        optimization_level: int=0, 
        iinj_clamps_to_use: Optional[List[float] | List[List[float]]]=None
    ) -> 'AnalyzerCfg':
        # Default values: (can't implement in usual default kwarg value because args are set equal to None in initialzation rather than being omitted)
        spreadsheets_input_dir = './inputs' if spreadsheets_input_dir is None else spreadsheets_input_dir
        image_save_dir = './outputs' if image_save_dir is None else image_save_dir
        iamge_save_type = 'png' if iamge_save_type is None else iamge_save_type
        timeseries_filter_cfgs = {
            "membrane_potentials": {"passband": 200, "ripple": 0.01, "stopband": 400, "attenuation": 80},
            "membrane_currents": {"passband": 40, "ripple": 0.01, "stopband": 60, "attenuation": 80},
            "activation_currents": {"passband": 40, "ripple": 0.01, "stopband": 60, "attenuation": 80}
        } if timeseries_filter_cfgs is None else timeseries_filter_cfgs
        files_to_analyze = 'ALL' if files_to_analyze is None else files_to_analyze
        optimization_level = 0 if optimization_level is None else optimization_level

        assert os.path.isdir(spreadsheets_input_dir)
        assert os.path.isdir(image_save_dir)
        assert iamge_save_type in {'png', 'emf'}
        if isinstance(files_to_analyze, list):
            assert all(os.path.splitext(x)[1] == ".xlsx" for x in files_to_analyze)
        else:
            assert isinstance(files_to_analyze, str)
            assert files_to_analyze == "ALL" or os.path.splitext(files_to_analyze)[1] == ".xlsx"
        assert optimization_level in {0, 1, 2}

        filtercfgs: Dict[str, FilterCfg] = {key:FilterCfg(**val) for key, val in timeseries_filter_cfgs.items()}

        return AnalyzerCfg(spreadsheets_input_dir, image_save_dir, iamge_save_type, filtercfgs, files_to_analyze, optimization_level, iinj_clamps_to_use)

    @property
    def paths_to_all_analysis_spreadsheets(self) -> List[str]:
        filenames: List[str]
        if isinstance(self.files_to_analyze, str):
            if self.files_to_analyze.upper() == 'ALL':
                filenames = list(filter(lambda x: os.path.splitext(x)[1] == ".xlsx", os.listdir(self.spreadsheets_input_dir)))
            else:
                filenames = [self.files_to_analyze]
        elif isinstance(self.files_to_analyze, list):
            filenames = self.files_to_analyze
        else:
            assert False, f"self.run is of incorrect type: {type(self.files_to_analyze)}"

        paths_to_spreadsheets: List[str] = [os.path.join(self.spreadsheets_input_dir, i) for i in filenames]
        
        assert all([os.path.exists(i) for i in paths_to_spreadsheets])
        assert len(paths_to_spreadsheets) > 0
        
        return paths_to_spreadsheets


class XLReader:
    def __init__(self, filepath: str):
        self.filepath = filepath
        _, filename_with_extension = os.path.split(self.filepath)
        self.filename, _ = os.path.splitext(filename_with_extension)
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
            if key != "times":
                data.rename(columns={key: f"{float(key):.3e}"}, inplace=True)                
        return data
    
    def get_paradigm_parameters(self, paradigm: str):
        paradigm = f"parameters_{paradigm}"
        df = pd.read_excel(self.filepath, sheet_name=paradigm, header=0)
        df = df.dropna(how='all').dropna(axis=1, how='all')
        df = df.dropna(how='any')
        assert all(df["Eact"] > df["Ess"]), f"Fix spreadsheet '{paradigm}': Eact must be greater than Ess"
        return df
    