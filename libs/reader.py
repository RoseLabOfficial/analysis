import pandas as pd
import numpy as np
from clamp import clamp_scm_quadratic as clamp

class reader:
    def __init__(self):
        self.filename:str = None
        self.data:float = None
        pass

    def read_xlsx(self, filename:str, sheet_names:str):
        self.filename = filename
        self.clamps = {}
        for sheet_name in sheet_names:
            properties = pd.read_excel(self.filename, str("parameters_")+sheet_name, engine='openpyxl')
            nClamps = len(properties["clamp"])
            self.clamps[sheet_name] = {}
            for i in range(nClamps):
                self.clamps[sheet_name][properties["clamp"][i]] = clamp(properties["clamp"][i], properties["Cm"][i], properties["Rin"][i], properties["Er"][i], properties["Ee"][i], properties["Ei"][i], properties["Ess"][i], properties["Eact"][i], properties["Et"][i], properties["xα"][i], properties["xβ"][i], properties["sps"][i])
            data = pd.read_excel(self.filename, sheet_name, engine='openpyxl')
            for i in range(nClamps):
                self.clamps[sheet_name][properties["clamp"][i]].insert_time(data["times"])
                self.clamps[sheet_name][properties["clamp"][i]].insert_Vm(data[self.clamps[sheet_name][properties["clamp"][i]].Iinj])
                self.clamps[sheet_name][properties["clamp"][i]].build()
        pass

if __name__ == "__main__":
    filename = "./libs/test.xlsx"
    rdr = reader()
    rdr.read_xlsx(filename, ["10pps"])