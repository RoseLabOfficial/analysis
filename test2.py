from libs.readers import XLReader
from libs.systems import WholeCellRecording, NonLinearModel
import numpy as np
import matplotlib.pyplot as plt
if __name__ == "__main__":
    filename = "./test_data/2012-05-03_2_baseline_prr_v4.xlsx"
    reader = XLReader(filename)
    model = NonLinearModel(reader.get_paradigm_parameters("10pps"), 0.0, 1.0, 1e4)
    v = model.make_voltage_grid()
    i = model.compute_current(v)
    fig, ax = plt.subplots(1, 1)
    ax.plot(v, i)
    ax.set_xlabel("voltage")
    ax.set_ylabel("current")
    ax.grid("on")
    plt.show()