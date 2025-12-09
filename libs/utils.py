from . import *

class SystemOperations:
    def __init__(self):
        pass

    def does_exist(self, dir_or_file_with_path: Path):
        return os.path.exists(dir_or_file_with_path)
    
    def list_files(self, dir: Path):
        if self.does_exist(dir):
            return os.listdir(dir)
        else:
            raise FileExistsError(f"{dir} does not exist!")

class Compliance:
    def __init__(self):
        pass

    def check_parameter_compliance(self, parameters: pd.DataFrame):
        return all(key in parameters for key in ("Iinj", "Cm", "Rin", "Er", "Ee", "Ei", "Et", "Ess", "xalpha", "xbeta", "sps", "Eref", "rate", "npulses", "amplitude"))

    def check_compliance(self, data: pd.DataFrame, parameters: pd.DataFrame):
        if self.check_parameter_compliance():
            for clamp in parameters["Iinj"]:
                if not clamp in data:
                    return False
        return True