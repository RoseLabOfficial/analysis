from . import *
from utils import SystemOperations

class Configurations:
    def __init__(self, configurations_json_address: Path):
        self.settings = jh.read(configurations_json_address)
        self.sysops_handler = SystemOperations()
    
    def get_input_directory(self):
        return Path(self.settings["inputdirs"])
    
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