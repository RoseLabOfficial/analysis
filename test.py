from libs.readers import XLReader
from libs.systems import WholeCellRecording
from libs.sysops import *
import numpy as np
import matplotlib.pyplot as plt
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Passing named arguments to analyzer.")
    parser.add_argument("--configs", "-c", type=str, default="./settings/batch_processing.json", help="json file specifying inputdir and outputdir")
    parser.add_argument("--list", "-l", type=str, default=None, help="list of all the files in input directory")
    parser.add_argument("--test", "-t", type=str, default=None, help="perform a test run to verify analysis.")
    args = parser.parse_args()
    print(args)
    if args.configs is not None:
        configs = Configurations(args.configs)
        files = configs.get_files_to_analyze()
        for file_id in files:
            reader = XLReader(file_id)
            for paradigm in reader.get_paradigms():
                data = WholeCellRecording(reader.get_paradigm_data(paradigm), 
                                reader.get_paradigm_parameters(paradigm))
                data.estimate_conductances(membrane_voltage_filter_cutoff=200, active_currents_filter_cutoff=100, membrane_currents_filter_cuoff=60)