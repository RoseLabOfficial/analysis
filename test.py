from libs.readers import XLReader
from libs.systems import WholeCellRecording
from libs.sysops import *
import numpy as np
import matplotlib.pyplot as plt
import argparse
import logging
import os
from pathlib import Path

def init_logging():
    current_file = os.path.splitext(os.path.basename(__file__))[0]
    logging.basicConfig(filename=current_file+".log", level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")

def load_configurations():
    configurations_json: Path = "./settings/configurations.json"
    return Configurations(configurations_json)
    

def main():
    configs = load_configurations()
    parser = argparse.ArgumentParser(description="Passing named arguments to analyzer")
    parser.add_argument("--show_configs", "-sc", action="store_true", help="show configurations")
    parser.add_argument("--list_inputs", "-li", action="store_true", help="list files from the input directory")
    parser.add_argument("--list_skips", "-ls", action="store_true", help="list files that are specified to be skipped in configurations")
    parser.add_argument("--list_files", "-la", action="store_true", help="list files that will be analyzed as per current configuraitons")
    parser.add_argument("--set_inputdir", "-si", type=Path, default ="./inputs", help="set input directory in configurations")
    parser.add_argument("--set_outputdir", "-so", type=Path, default ="./outputs", help="set output directory in configurations")
    args = parser.parse_args()
    if args.show_configs:
        print(configs.settings)
    if args.list_inputs:
        print(configs.get_input_files())
    if args.list_skips:
        print(configs.get_files_to_skip())
    if args.list_files:
        print(configs.get_files_to_analyze())
    if args.set_inputdir is not None:
        configs.set_input_directory(args.set_inputdir)
        print(configs.settings)
    if args.set_outputdir is not None:
        configs.set_output_directory(args.set_outputdir)
        print(configs.settings)
    
        
    
    

if __name__ == "__main__":
    main()
    # parser = argparse.ArgumentParser(description="Passing named arguments to analyzer.")
    # parser.add_argument("--show_configs", "-sc", action="store_file", help="show configurations")
    # parser.add_argument("--show_inputs", "-si", action="store_file
    # parser.add_argument("--configs", "-c", type=str, default="./settings/batch_processing.json", help="json file specifying inputdir and outputdir")
    # parser.add_argument("--list", "-l", type=str, default=None, help="list of all the files in input directory")
    # parser.add_argument("--test", "-t", type=str, default=None, help="perform a test run to verify analysis.")
    # args = parser.parse_args()
    # if args.configs is not None:
    #     configs = Configurations(args.configs)
    #     files = configs.get_files_to_analyze()
    #     for file_id in files:
    #         print(configs.get_input_directory()+"/"+file_id)
            # reader = XLReader(file_id)
            # for paradigm in reader.get_paradigms():
            #     data = WholeCellRecording(reader.get_paradigm_data(paradigm), 
            #                     reader.get_paradigm_parameters(paradigm))
            #     data.estimate_conductances(membrane_voltage_filter_cutoff=200, active_currents_filter_cutoff=100, membrane_currents_filter_cuoff=60)
