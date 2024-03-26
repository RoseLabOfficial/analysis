import json
import os
from pathlib import Path
from datetime import datetime

class SystemOperations:
    def make_directory(self, directory):
        if not os.path.exists(directory):
            os.makedirs(directory)

    def get_time_as_string(self):
        now = datetime.now()
        return now.strftime("_%Y_%m_%d_%H_%M_%S")
    
    def make_timed_directory(self, dir: Path, prefix: str) -> Path:
        now = datetime.now()
        dir_name = now.strftime("_%Y_%m_%d_%H_%M_%S")
        new_directory_path = dir / (prefix + dir_name)
        if not new_directory_path.exists():
            new_directory_path.mkdir(parents=True)
            print(f"Directory created: {new_directory_path}")
        else:
            print(f"Directory already exists: {new_directory_path}")
        return new_directory_path
    
    def read_json(self, jsonfile: str):
        with open(jsonfile, "r") as settings:
            configs = json.load(settings)
        return configs

    def write_json(self, jsonfile: str, json_data: dict):
        with open(jsonfile, "w") as settings:
            configs = json.dump(json_data, settings, indent=4)
        pass
    
    def check_directory(self, directory):
        return os.path.exists(directory)
    
    def list_files(self, directory):
        return os.listdir(directory)
    
    def split_file_address(self, filepaths: Path):
        processed_filepaths = []
        for filepath in filepaths:
            path = filepath.parent
            stem = filepath.stem
            suffix = filepath.suffix
            if suffix == ".xlsx":
                processed_filepaths.append([path, stem, suffix])
        return processed_filepaths

class Configurations:
    def __init__(self, configs_json: str):
        self.configs_json = configs_json
        self.sysops = SystemOperations()
        self.settings = self.sysops.read_json(configs_json)
        pass

    def set_input_directory(self, inputdir: Path):
        self.settings["inputdir"] = str(inputdir)
        self.sysops.write_json(self.configs_json, self.settings)
        pass

    def get_input_directory(self) -> Path:
        return Path(self.settings["inputdir"])

    def set_output_directory(self, outputdir: Path):
        self.settings["outputdir"] = str(outputdir)
        self.sysops.write_json(self.configs_json, self.settings)
        pass
    
    def get_output_directory(self) -> Path:
        return Path(self.settings["outputdir"])

    def get_input_files(self):
        return self.sysops.list_files(self.get_input_directory())
    
    def get_files_to_skip(self):
        return self.settings["skipfiles"]
    
    def get_files_to_analyze(self):
        return [filename for filename in self.get_input_files() if filename not in self.get_files_to_skip()]
    
    def get_filters(self):
        return self.settings["filters"]