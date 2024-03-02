from libs.readers import XLReader
from libs.systems import WholeCellRecording
import numpy as np
if __name__ == "__main__":
    filename = "./test_data/2012-05-03_2_baseline_prr_v4.xlsx"
    reader = XLReader(filename)
    data = WholeCellRecording(reader.get_paradigm_data("10pps"), 
                              reader.get_paradigm_parameters("10pps"))
    data.scale_data()
    data.compute_polarizations()
    data.compute_active_conductance_constants()
    data.compute_active_currents()
    data.compute_leakage_currents()
    print(data.compute_membrane_currents())