from libs.systems import Analyzer
from libs.readers import Configs

import json, os, sys
from argparse import ArgumentParser

from typing import Any, Dict

if __name__ == "__main__":
    if len(sys.argv) == 2:
        with open(sys.argv[1], 'r') as f:
            cfg_kwargs: Dict[str, Any] = json.load(f)

    elif len(sys.argv) > 2:
        parser: ArgumentParser = ArgumentParser(description="Passing named arguments to analyzer")
        parser.add_argument("--set_inputdir", "-si", type=str, default='./inputs', help="set input directory in configurations")
        parser.add_argument("--set_outputdir", "-so", type=str, default='./outputs', help="set output directory in configurations")
        parser.add_argument("--run", "-r", action="store_true", help="Run all the files in a dictory and put the results in outputs")
        parser.add_argument("--run_file", "-rf", type=str, default='ALL', help="run a single file")
        parser.add_argument("--optimization_level", "-o", type=int, default=0, help="0: User-specified values of Eact, 1: Compute Eact based on PR with max depolarization, 2: Compute Eact for each PR")
        args = parser.parse_args()
        args.run_file = 'ALL' if args.run else args.run_file 

        cfg_kwargs = dict(
            input_directory=args.set_inputdir, 
            output_directory=args.set_outputdir, 
            run=args.run_file, 
            optimization_level=args.optimization_level
        )
    
    else:
        with open('./settings/js_configs.json', 'r') as f:
            cfg_kwargs = json.load(f)

    cfg: Configs = Configs.new(**cfg_kwargs)
    assert all([os.path.exists(i) for i in cfg.full_run_paths])

    analyzer = Analyzer(cfg.full_run_paths, cfg.output_directory)
    analyzer.run(**cfg.analyzer_run_kwargs)