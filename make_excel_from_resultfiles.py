import pandas as pd
import os, re


from typing import List, Any, Set, Dict
from dataclasses import dataclass

def numeric(s):
    match = re.search(r'-?\d*\.?\d+', s)
    if match:
        return match.group()
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
    stimulus: bool

    def new(filename: str) -> "WavAvg":
        args: List[str] = os.path.splitext(os.path.basename(filename))[0].split('_')
        return WavAvg(
            filename=filename,
            atten=numeric(args[0]),
            carrier_freq=numeric(args[1]),
            pulse_dur=numeric(args[2]),
            pulse_n=numeric(args[3]),
            pps=numeric(args[4]),
            iinj=numeric(args[5]),
            version=args[6],
            stimulus=len(args) >= 8
        )

    def compat(self, other: "WavAvg") -> bool:
        if any((
            self.atten != other.atten,
            self.carrier_freq != other.carrier_freq,
            self.pulse_dur != other.pulse_dur,
            self.pulse_n != other.pulse_n
        )):
            return False
        return True

def txt_to_excel(directory, output_folder, version):
    data: Dict[WavAvg, Any] = {}

    for filename in os.listdir(directory):
        if filename.endswith('.txt'):
            file_path = os.path.join(directory, filename)
            data[WavAvg.new(filename)] = pd.read_csv(file_path, sep='\t', engine='python')
    
    compat_groups: List[Set[WavAvg]] = []
    for wavavg in data.keys():
        if all([wavavg not in group for group in compat_groups]):
            compat_groups.append(set())
        for other_wavavg in data.keys():
            if all([other_wavavg not in group for group in compat_groups]) and wavavg.compat(other_wavavg):
                compat_groups[-1].add(other_wavavg)
    
    for i, group in enumerate(compat_groups):
        writer = pd.ExcelWriter(os.path.join(output_folder, str(i)+'.xlsx'), engine='openpyxl', mode='w')
        print_group = list(group)
        print('\n'.join([elem.filename for elem in print_group]))
        print_group = print_group[0]
        print(os.path.join(output_folder, str(i)+'.xlsx')+'\n'+f"ATTEN: {print_group.atten} FREQ: {print_group.carrier_freq} PD: {print_group.pulse_dur} N_PULSES: {print_group.pulse_n}")
        write_data: Dict[str, pd.DataFrame] = {}
        for key in group:
            sheet_name = f"{key.pps}pps"
            if sheet_name not in write_data:
                write_data[sheet_name] = pd.DataFrame()
                write_data[sheet_name]['times'] = data[key]['Time']
            if key.stimulus:
                write_data[sheet_name]['stimulus'] = data[key]['1 Stimulus']
            else:
                write_data[sheet_name][f"{key.iinj}~{key.version}"] = data[key]['1 Spikes+hr']
            
        for sheet, df in write_data.items(): 
            if version is None:
                confstr = input(f"\tUse from {sheet}:\n" + ''.join(f'\t\t{j}: {version}'+'\n' for j, version in list(enumerate(df.keys())) if version != "times" or version != "stimulus")).split(',')
            else:
                confstr = [j for j, elem in enumerate(df.keys()) if elem.split('~')[-1] == version or elem == 'times' or elem == 'stimulus']
            if len(confstr)-1:   
                confirmed = [int(i) for i in confstr]
                for j, colname in enumerate(df.keys()):
                    if colname != 'times' and colname != 'stimulus':
                        if not j in confirmed:
                            print(list(df.keys()))
                            df.drop(colname, axis=1, inplace=True)
                        else:
                            df.rename(columns={colname: float(colname.split('~')[0])*1e-9}, inplace=True)
            print('\n')
            df.to_excel(writer, sheet_name=sheet, index=False)
        writer._save()

# Specify the directory containing the .txt files and the desired output Excel file
directory_path = '/home/roselab/Documents/Docs/JS/2016-10-26_1/resultfiles/baseline'
output_excel_file = 'output_data'

# Call the function
txt_to_excel(directory_path, output_excel_file, 'v1')
