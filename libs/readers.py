import pandas as pd
import os
from pathlib import Path
import json
from utils import *

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
    
class JsonReader:
    def __init__(self, json_file_with_path: Path):
        self.file_with_path = json_file_with_path
        self.file_type_checker = FileTypeChecker()
        pass

    def read(self):
        if self.file_type_checker.is_json(self.file_with_path):
            with open(self.file_with_path, "r") as file_pointer:
                content = json.load(file_pointer)
            return content
        else:
            return TypeError(f"{self.file_with_path} is not a json file!")

class ConfigurationsReader(JsonReader):
    def __init__(self, configurations_file_with_path: Path):
        super.__init__(configurations_file_with_path)
        self.settings = self.read()
        self.system_operations = SystemOperations()
        pass

    def get_input_directory(self) -> Path:
        return Path(self.settings["inputdir"])
    
    def get_output_directory(self) -> Path:
        return Path(self.settings["outputdir"])
    
    def list_input_files(self):
        return self.system_operations.list_files(self.settings["inputdir"])
    
    def list_output_files(self):
        return self.system_operations.list_files(self.settings["outputdir"])