import tensorflow as tf
import numpy as np

class clamp_scm_quadratic:
    def __init__(self, Fs:float, T:float, Iinj:float, Cm:float, Rin:float, Er:float, Ee:float, Ei:float, Ess:float, Eact:float, Et:float, xalpha:float, xbeta:float, sps:float):
        self.Fs = Fs
        self.T = T
        self.N = int (self.Fs*self.T)
        self.time = None
        self.Iinj = Iinj
        self.Cm = Cm
        self.Rin = Rin
        self.Er = Er
        self.Ee = Ee
        self.Ei = Ei
        self.Ess = Ess
        self.Eact = Eact
        self.Et = Et
        self.xalpha = xalpha
        self.xbeta = xbeta
        self.sps = sps
        pass

    def insert_time(self, time:float):
        if time is None:
            self.time = tf.reshape(tf.linspace(0, self.T, self.N, dtype=tf.float32), [1, self.N])
        else:
            self.time = tf.constant(time, dtype=tf.float32)
        pass

    def insert_Vm(self, Vm:float):
        self.Vm = tf.reshape(Vm, self.Vm.shape)
        pass

    def build(self):
        self.Vm = tf.zeros(self.time.shape, dtype=tf.float32)+self.Ess
        self.dVm_by_dt = tf.zeros(self.time.shape, dtype=tf.float32)
        self.Im = tf.zeros(self.time.shape, dtype=tf.float32)
        self.Iactive = tf.zeros(self.time.shape, dtype=tf.float32)
        self.Ileak = tf.zeros(self.time.shape, dtype=tf.float32)
        if self.time == None:
            self.time = tf.reshape(tf.linspace(0, self.T, self.N, dtype=tf.float32), [1, self.N])
        pass

    def compute_dVm_by_dt(self):
        Vm_leading = tf.concat([tf.constant(self.Ess, dtype=tf.float32), self.Vm[:, :-2]], axis=-1)
        time_leading = tf.concat([tf.constant(self.time[:, 0] - self.time[:, 1], dtype=tf.float32), self.time[:, :-2]], axis=-1)
        self.dVm_by_dt = (self.Vm - Vm_leading) / (self.time-time_leading)
        pass
        
    def compute_Im(self):
        self.Im = self.Cm*self.compute_dVm_by_dt()
        pass

    def compute_Ileak(self):
        self.Ileak = (1/self.Rin)*(self.Vm - self.Er)
        pass

    def compute_Iactive(self):
        self.Iactive = self.xalpha*(self.Vm - self.Et)*(self.Vm - self.Ess) - self.xbeta*(self.Vm-self.Ess)
        pass
        