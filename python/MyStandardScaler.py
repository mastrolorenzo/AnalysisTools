#!/usr/bin/env python
import numpy as np

class StandardScaler(object):
    def __init__(self):
        self.mean_ = None
        self.scale_ = None

    def get_params(self):
        return [self.mean_, self.scale_]

    def fit(self, data):
        self.mean_ = np.mean(np.nan_to_num(data), axis=0)
        self.scale_ = np.sqrt(np.var(np.nan_to_num(data), axis=0))
        scale = self.scale_.copy()
        scale[scale == 0.0] = 1.0
        self.scale_ = scale

    def transform(self, data):
        dataScaled = (data - self.mean_) / self.scale_
        return dataScaled

    def fit_transform(self, data):
        self.fit(data)
        return self.transform(data)
