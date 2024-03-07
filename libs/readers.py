import pandas as pd

class XLReader:
    def __init__(self, filename: str):
        self.filename = filename
        self.data_pointer = pd.ExcelFile(self.filename)
        pass

    def get_sheet_names(self):
        return self.data_pointer.sheet_names

    def get_paradigms(self):
        sheet_names = self.get_sheet_names()
        return [sheet_name for sheet_name in sheet_names if "parameters" not in sheet_name if "stats" not in sheet_name]

    def get_paradigm_data(self, paradigm: str):
        return pd.read_excel(self.filename, sheet_name=paradigm, header=0)
    
    def get_paradigm_parameters(self, paradigm: str):
        paradigm = "parameters_"+paradigm
        return pd.read_excel(self.filename, sheet_name=paradigm, header=0)
