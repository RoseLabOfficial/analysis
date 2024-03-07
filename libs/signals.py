import numpy as np

class TimeSeries:
    def __init__(self, tstart: float, tend: float, sampling_rate: float, dtype: np.dtype = np.float64):
        self.dtype = dtype
        self.tstart = self.dtype(tstart)
        self.tend = self.dtype(tend)
        self.sampling_rate = self.dtype(sampling_rate)
        pass

    def get_total_samples(self) -> int:
        return int(self.sampling_rate*np.abs(self.tend - self.tstart))

    def make_time(self):
        nsamples = self.get_total_samples()
        return np.linspace(start = self.tstart, stop = self.tend, num = nsamples, dtype = self.dtype)