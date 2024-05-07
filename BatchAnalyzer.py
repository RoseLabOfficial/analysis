from libs.readers import AnalysisConfigurations
from libs.systems import Analyzer
import argparse
from pathlib import Path

def load_configurations():
    configurations_json: Path = "./settings/configurations.json"
    return AnalysisConfigurations(configurations_json)

def main():
    configs = load_configurations()
    parser = argparse.ArgumentParser(description="Passing named arguments to analyzer")
    parser.add_argument("--show_configs", "-sc", action="store_true", help="show configurations")
    parser.add_argument("--list_inputs", "-li", action="store_true", help="list files from the input directory")
    parser.add_argument("--list_skips", "-ls", action="store_true", help="list files that are specified to be skipped in configurations")
    parser.add_argument("--list_files", "-la", action="store_true", help="list files that will be analyzed as per current configuraitons")
    parser.add_argument("--set_inputdir", "-si", type=Path, default = None, help="set input directory in configurations")
    parser.add_argument("--set_outputdir", "-so", type=Path, default = None, help="set output directory in configurations")
    parser.add_argument("--run", "-r", action="store_true", help="Run all the files in a dictory and put the results in outputs")
    parser.add_argument("--run-file", "-rf", type=Path, default = None, help="run a single file")
    parser.add_argument("--run-user-files", "-ru", action="store_true", help="run user files from configuration")
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
    if args.run:
        inputdir = configs.get_input_directory()
        inputfiles = configs.get_input_files()
        files = [inputdir / filename for filename in inputfiles]
        analyzer = Analyzer(files, configs.get_output_directory(), configs.get_filters())
        analyzer.run()
    if args.run_file is not None:
        inputdir = configs.get_input_directory()
        filepath = inputdir / args.run_file
        if filepath.exists():
            analyzer = Analyzer([filepath], configs.get_output_directory(), configs.get_filters())
            analyzer.run(configs.get_optimization()["level"])
        else:
            BufferError(f"{filepath} doesn't exist!")
    if args.run_user_files is not None:
        inputdir = configs.get_input_directory()
        filepaths = [inputdir / x for x in configs.get_user_files()]
        for filepath in filepaths:
            if filepath.exists():
                analyzer = Analyzer([filepath], configs.get_output_directory(), configs.get_filters())
                analyzer.run(configs.get_optimization()["level"])
            else:
                Warning(f"{filepath} does not exist. Skipping!")

if __name__ == "__main__":
    main()