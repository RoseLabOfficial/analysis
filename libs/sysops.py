import json
import os

class SystemOperations:
    def make_directory(self, directory):
        if not os.path.exists(directory):
            os.makedirs(directory)
    
    def read_config_json(self, jsonfile: str):
        with open(jsonfile, "r") as settings:
            configs = json.load(settings)
        return configs
    
    def check_directory(self, directory):
        return os.path.exists(directory)
    
    def list_files(self, directory):
        return os.listdir(directory)

class Configurations:
    def __init__(self, configs_json: str):
        self.sysops = SystemOperations()
        self.configs = self.sysops.read_config_json(configs_json)
        pass

    def get_input_directory(self):
        return self.configs["inputdir"]
    
    def get_output_director(self):
        return self.configs["outputdir"]
    
    def get_files_to_skip(self):
        return self.configs["skiplist"]
    
    def get_files_to_analyze(self):
        files = self.sysops.list_files(self.get_input_directory())
        return [filename for filename in files if filename not in self.get_files_to_skip()]