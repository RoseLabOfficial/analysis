import sys
from pathlib import Path

from libs.systems import Analyzer
from libs.readers import AnalyzerCfg

if __name__ == "__main__":
    assert len(sys.argv) == 2
    cfg: AnalyzerCfg = AnalyzerCfg.from_json(Path(sys.argv[1]))
    analyzer = Analyzer(cfg)
    analyzer.run()