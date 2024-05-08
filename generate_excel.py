import pandas as pd
import os, re, sys

from typing import List, Optional, Dict
from dataclasses import dataclass

def numeric(s: str) -> float | int:
    match = re.search(r'-?\d*\.?\d+', s)
    if match:
        group = match.group()
        if '.' in group:
            return float(group)
        return int(group)
    return None

@dataclass(frozen=True)
class WavAvg:
    filename: str
    atten: float
    carrier_freq: float
    pulse_dur: float
    pulse_n: int
    pps: float
    iinj: float
    version: str
    is_stim: bool

    def new(filename: str) -> "WavAvg":
        fields: List[str] = os.path.splitext(os.path.basename(filename))[0].split('_')
        args = [numeric(i) for i in fields[:6]]
        return WavAvg(filename, *args, version=fields[6], is_stim=len(fields) >= 8)

    def compat(self, other: "WavAvg") -> bool:
        if any((
            self.atten != other.atten,
            self.carrier_freq != other.carrier_freq,
            self.pulse_dur != other.pulse_dur,
            self.pulse_n != other.pulse_n
        )):
            return False
        return True

def main(dir: str, output_dir: str, version: Optional[str]=None) -> None:
    data: Dict[WavAvg, pd.DataFrame] = {
        WavAvg.new(filename):pd.read_csv(os.path.join(dir, filename), sep='\t', engine='python')
        for filename in os.listdir(dir) if filename.endswith('.txt')
    }

    compat_groups: List[List[WavAvg]] = []
    not_in = lambda x, nest: all([not x in sub for sub in nest])
    for wavavg in data.keys():
        if not_in(wavavg, compat_groups):
            compat_groups.append([])
        for other_wavavg in data.keys():
            if not_in(other_wavavg, compat_groups) and wavavg.compat(other_wavavg):
                compat_groups[-1].append(other_wavavg)
    
    for group in compat_groups:
        output_filepath: str = os.path.join(output_dir, f"{group[0].filename[:group[0].filename.find('pulses')+len('pulses')]}.xlsx")
        writer = pd.ExcelWriter(output_filepath, engine='openpyxl', mode='w')

        print('\n'+'\n'.join([wavavg.filename for wavavg in group])+f"{output_filepath}\n")
        
        write_data: Dict[str, pd.DataFrame] = {}
        for wavavg in group:
            sheet: str = f"{wavavg.pps}pps"

            if not sheet in write_data:
                write_data[sheet] = pd.DataFrame()
                write_data[sheet][('times', None)] = data[wavavg]['Time']

            if wavavg.is_stim and (wavavg.version == version if version else True):
                if sum(data[wavavg]['1 Stimulus']) != 0:
                    write_data[sheet][('stimulus', None)] = data[wavavg]['1 Stimulus']

            elif not wavavg.is_stim:
                write_data[sheet][(wavavg.iinj, wavavg.version)] = data[wavavg]['1 Spikes+hr']
            
        for sheet, df in write_data.items(): 
            if version:
                chosen_cols: List[int] = [i for i, key in enumerate(df.keys()) if key[1] == version]
            else:
                msg: str = f"\tUse from {sheet}:\n" + ''.join(f'\t\t{i}: {v}'+'\n' for i, key in enumerate(df.keys()) if key[1] is not None)
                chosen_cols = [int(i) for i in input(msg).replace(' ', '').split(',') if i.isdigit()]

            if len(chosen_cols) != 0:
                for i, key in enumerate(df.keys()):
                    if i in chosen_cols or key[1] is None:
                        df.rename(columns={key: key[0] * 1e-9 if not isinstance(key[0], str) else key[0]}, inplace=True)
                    else:
                        df.drop(key, axis=1, inplace=True)
            
            df = df[sorted(list(df.keys()), key=lambda x: -x if not isinstance(x, str) else float('-inf') if x == 'times' else float('inf'))]

            df.to_excel(writer, sheet_name=sheet, index=False)

        writer._save()

if __name__ == '__main__':
    main(*sys.argv[1:4])

