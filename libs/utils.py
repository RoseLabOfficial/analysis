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

class FileTypeChecker:
    def __init__(self):
        self.system_operations = SystemOperations()
        pass

    def is_xlsx(self, file_with_path: Path):
        if self.system_operations.does_exist(file_with_path) and file_with_path.endswith('.xlsx'):
            try:
                openpyxl.load_workbook(file_with_path)
                return True
            except openpyxl.utils.exceptions.InvalidFileException:
                return False
        return False

    def is_json(self, file_with_path: Path):
        if self.system_operations.does_exist(file_with_path):
            try:
                with open(file_with_path, 'r') as file:
                    json.load(file)
                return True
            except (json.JSONDecodeError, UnicodeDecodeError):
                return False
        return False

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