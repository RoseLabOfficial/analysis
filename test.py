from libs.readers import XLReader
from libs.systems import WholeCellRecording
import numpy as np
import matplotlib.pyplot as plt
if __name__ == "__main__":
    filename = "./test_data/2012-05-03_2_baseline_prr_v4.xlsx"
    reader = XLReader(filename)
    for paradigm in reader.get_paradigms():
        print(f"Analyzing paradigm: {paradigm}")
        data = WholeCellRecording(reader.get_paradigm_data(paradigm), 
                                reader.get_paradigm_parameters(paradigm))
        data.estimate_conductances(membrane_voltage_filter_cutoff=200, active_currents_filter_cutoff=100, membrane_currents_filter_cuoff=60)
        data.plot("./test_results/"+paradigm)
        print(f"done.")