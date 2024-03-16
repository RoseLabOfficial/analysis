import json
import os
from pathlib import Path

class SystemOperations:
    def make_directory(self, directory):
        if not os.path.exists(directory):
            os.makedirs(directory)
    
    def read_json(self, jsonfile: str):
        with open(jsonfile, "r") as settings:
            configs = json.load(settings)
        return configs

    def write_json(self, jsonfile: str, json_data: dict):
        with open(jsonfile, "w") as settings:
            configs = json.dump(json_data, indent=4)
        pass
    
    def check_directory(self, directory):
        return os.path.exists(directory)
    
    def list_files(self, directory):
        return os.listdir(directory)

class Configurations:
    def __init__(self, configs_json: str):
        self.configs_json = configs_json
        self.sysops = SystemOperations()
        self.settings = self.sysops.read_json(configs_json)
        pass

    def set_input_directory(self, inputdir: Path):
        self.settings["inputdir"] = inputdir
        self.sysops.write_json(self.configs_json, self.settings)
        pass

    def get_input_directory(self):
        return self.settings["inputdir"]

    def set_output_directory(self, outputdir: Path):
        self.settings["outputdir"] = outputdir
        self.sysops.write_json(self.configs_json, self.settings)
        pass
    
    def get_output_director(self):
        return self.settings["outputdir"]

    def get_input_files(self):
        return self.sysops.list_files(self.get_input_directory())
    
    def get_files_to_skip(self):
        return self.settings["skiplist"]
    
    def get_files_to_analyze(self):
        return [filename for filename in self.get_input_files() if filename not in self.get_files_to_skip()]
