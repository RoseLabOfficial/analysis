from . import *
from libs.utils import SystemOperations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

@dataclass(frozen=True)
class FilterCfg:
    passband: float 
    ripple: float 
    stopband: float 
    attenuation: float

@dataclass(frozen=True)
class Configs:
    input_directory: str 
    output_directory: str 

    filters: Dict[str, FilterCfg]

    run: str | List[str]='ALL' # ALL runs all files in input directory. Single str runs just that file. List[str] runs all specified files.
    optimization_level: int=0 # 0: No optimization 1: Calculate Eact based on max depolarization paradigm 2: Calculate Eact for each paradigm.
    use_clamps: Optional[List[float]]=None # If None, uses all current clamps. If List[float], uses specified clamps

    @classmethod
    def new(
        cls, input_directory: Optional[str]=None, output_directory: Optional[str]=None,
        filters: Optional[Dict[str, Dict[str, float]]]=None, run: Optional[str | List[str]]=None, 
        optimization_level: int=0, use_clamps: Optional[List[float]]=None
    ) -> 'Configs':
        input_directory = './inputs' if input_directory is None else input_directory
        output_directory = './outputs' if output_directory is None else output_directory

        filters = {
            "membrane_potentials": {"passband": 200, "ripple": 0.01, "stopband": 400, "attenuation": 80},
            "membrane_currents": {"passband": 40, "ripple": 0.01, "stopband": 60, "attenuation": 80},
            "activation_currents": {"passband": 40, "ripple": 0.01, "stopband": 60, "attenuation": 80}
        } if filters is None else filters
        filtercfgs: Dict[str, FilterCfg] = {key:FilterCfg(**val) for key, val in filters.items()}

        run = 'ALL' if run is None else run
        optimization_level = 0.0 if optimization_level is None else optimization_level

        return Configs(input_directory, output_directory, filtercfgs, run, optimization_level, use_clamps)
    
    @property
    def analyzer_run_kwargs(self) -> Dict[str, int | Optional[List[float]] | Dict[str, FilterCfg]]:
        return {
            "optimization_level": self.optimization_level,
            "current_clamps": self.use_clamps, 
            "filter_configurations": self.filters
        }

    @property
    def full_run_paths(self) -> List[str]:
        if self.run.upper() == 'ALL':
            filenames: List[str] = os.listdir(self.input_directory)
        elif isinstance(self.run, list):
            filenames = self.run
        else:
            filenames = [self.run]
        return [os.path.join(self.input_directory, i) for i in filenames]

class XLReader:
    def __init__(self, filepath: str):
        self.filepath = filepath
        _, filename_with_extension = os.path.split(self.filepath)
        self.filename, _ = os.path.splitext(filename_with_extension)
        self.data_pointer = pd.ExcelFile(self.filepath)
        pass

    def get_sheet_names(self):
        return self.data_pointer.sheet_names

    def get_paradigms(self):
        sheet_names = self.get_sheet_names()
        return [sheet_name for sheet_name in sheet_names if "parameters" not in sheet_name if "stats" not in sheet_name]

    def get_paradigm_data(self, paradigm: str):
        return pd.read_excel(self.filepath, sheet_name=paradigm, header=0)
    
    def get_paradigm_parameters(self, paradigm: str):
        paradigm = "parameters_"+paradigm
        df = pd.read_excel(self.filepath, sheet_name=paradigm, header=0)
        df = df.dropna(how='all').dropna(axis=1, how='all')
        return df

class AnalysisConfigurations:
    def __init__(self, configurations_json_address: Path):
        self.settings = jh.read(configurations_json_address)
        self.sysops_handler = SystemOperations()
    
    def get_input_directory(self):
        return Path(self.settings["inputdir"])
    
    def set_input_directory(self, newdir: str) -> None:
        self.settings["inputdir"] = str(newdir)
    
    def get_output_directory(self):
        return Path(self.settings["outputdir"])
    
    def set_output_directory(self, newdir: str) -> None:
        self.settings["outputdir"] = str(newdir)
    
    def get_input_files(self):
        return self.sysops_handler.list_files(self.get_input_directory())
    
    def get_user_files(self):
        return self.settings["userfiles"]
    
    def get_ignore_files(self):
        return self.settings["ignorefiles"]
    
    def get_optimiizer_level(self):
        return self.settings["optimizer"]["level"]
    
    def get_filters(self):
        return self.settings["filters"]