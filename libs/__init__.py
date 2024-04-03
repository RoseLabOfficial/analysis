import pandas as pd
import numpy as np
import os
import openpyxl
import logging.config
import json
from pathlib import Path

def init_logging(logger_configurations):
    logging.config.dictConfig(logger_configurations)

class JsonHandler:
    def __init__(self):
        pass

    def read(self, fileaddress: Path):
        try:
            with open(fileaddress, 'r') as filepointer:
                return json.load(filepointer)
        except FileNotFoundError:
            raise FileNotFoundError(f'File not found: {fileaddress}')
        except json.JSONDecodeError as e:
            raise ValueError(f'Failed to decode JSON from {fileaddress}: {e.msg}')

jh = JsonHandler()
logger_configurations = jh.read("./settings/logging.json")
init_logging(logger_configurations)
analysis_logger = logging.getLogger('Analyzer')