import sys
from pathlib import Path

from libs.systems import Analyzer
from libs.readers import AnalyzerCfg

if __name__ == "__main__":
    assert len(sys.argv) in {2, 3}
    cfg: AnalyzerCfg = AnalyzerCfg.from_json(Path(sys.argv[1]))
    analyzer: Analyzer = Analyzer(cfg)
    if len(sys.argv) == 3:
        assert isinstance(sys.argv[2], bool)
        analyzer.run(sys.argv[2])
    else:
        analyzer.run(True)