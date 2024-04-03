from . import *
from libs.utils import SystemOperations

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
        return pd.read_excel(self.filepath, sheet_name=paradigm, header=0)


class AnalysisConfigurations:
    def __init__(self, configurations_json_address: Path):
        self.settings = jh.read(configurations_json_address)
        self.sysops_handler = SystemOperations()
    
    def get_input_directory(self):
        return Path(self.settings["inputdir"])
    
    def get_output_directory(self):
        return Path(self.settings["outputdir"])
    
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