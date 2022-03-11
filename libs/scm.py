import numpy as np

class clamp:
    def __init__(self, Ess=-65e-3, Eact=-45e-3, Et=-40e-3, Iinj=0e9, xalpha=1, xbeta=1, precision=np.float32):
        self.precision = precision
        self.Ess = self.precision(Ess)
        self.Eact = self.precision(Eact)
        self.Et = self.precision(Et)
        self.Iinj = self.precision(Iinj)
        self.Vm = np.zeros(shape=(0), dtype=self.precision)
        self.alpha = self.precision(0)
        self.beta = self.precision(0)
        self.xalpha = self.precision(xalpha)
        self.xbeta = self.precision(xbeta)
        pass

    def set_Vm(self, Vm):
        self.Vm = np.asarray(Vm, dtype=self.precision)
        pass

    def correct_to_Ess(self):
        if self.Vm.shape[0] > 0:
            self.Vm = self.Vm - self.Vm[0, 0]
            self.Vm = self.Vm + self.Ess
            return 1
        return 0

    def correct_Ess(self):
        if self.Vm.shape[0] > 0:
            self.Ess = self.Vm[0, 0]
            return 1
        return 0

    def compute_active_conductances(self, Rin):
        self.alpha = (1/self.precision(Rin))/(self.precision(2)*(self.Eact - self.Ess))
        self.beta = self.alpha*(self.Et - self.Ess)
        self.alpha = self.alpha*self.xalpha
        self.beta = self.beta*self.xbeta
        return 1

class neuron:
    def __init__(self, Cm=0.14e-11, Rin=0.5e9, Ee=-30e-3, Ei=-100e-3, precision=np.float32):
        self.precision = precision
        self.Cm = self.precision(Cm)
        self.Rin = self.precision(Rin)
        self.Ee = self.precision(Ee)
        self.Ei = self.precision(Ei)
        self.clamps = []
        pass

    def create_new_clamp(self, Ess=-65e-3, Eact=-45e-3, Et=-40e-3, Iinj=0e9, xalpha=1, xbeta=1):
        c = clamp(Ess, Eact, Et, Iinj, xalpha, xbeta, precision=self.precision)
        return self.insert_clamp(c)
    
    def insert_clamp(self, c: clamp):
        for clamp in self.clamps:
            if c.Iinj == clamp.Iinj:
                clamp = c
                return 1
        self.clamps.append(c)
        return 2
    
    def delete_clamp(self, Iinj):
        index = None
        for index, clamp in enumerate(self.clamps):
            if clamp.Iinj == self.precision(Iinj):
                break
        if index is None:
            return 0
        else:
            del self.clamps[index]
            return 1
    
    def update_Iinj(self, Iinj, Iinj_new):
        for clamp in self.clamps:
            if clamp.Iinj == clamp.precision(Iinj):
                clamp.Iinj = clamp.precision(Iinj_new)
                return 1
        return 0
    
    def update_Ess(self, Iinj, Ess):
        for clamp in self.clamps:
            if clamp.Iinj == clamp.precision(Iinj):
                clamp.Ess = clamp.precision(Ess)
                return 1
        return 0
    
    def update_Eact(self, Iinj, Eact):
        for clamp in self.clamps:
            if clamp.Iinj == clamp.precision(Iinj):
                clamp.Eact = clamp.precision(Eact)
                return 1
        return 0
    
    def update_Et(self, Iinj, Et):
        for clamp in self.clamps:
            if clamp.Iinj == clamp.precision(Iinj):
                clamp.Et = clamp.precision(Et)
                return 1
        return 0

    def update_xalpha(self, Iinj, xalpha):
        for clamp in self.clamps:
            if clamp.Iinj == clamp.precision(Iinj):
                clamp.xalpha = clamp.precision(xalpha)
                return 1
        return 0
    
    def update_xbeta(self, Iinj, xbeta):
        for clamp in self.clamps:
            if clamp.Iinj == clamp.precision(Iinj):
                clamp.xbeta = clamp.precision(xbeta)
                return 1
        return 0
    