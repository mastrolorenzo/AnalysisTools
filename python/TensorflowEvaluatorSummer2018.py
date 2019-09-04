#!/usr/bin/env python
from TensorflowDNNClassifierSummer2018 import TensorflowDNNClassifier
import numpy as np
import pickle

# reads tensorflow checkpoints and applies DNN output to ntuples,
# including variations from systematic uncertainties
# needs: tensorflow >=1.4
class TensorflowEvaluator():

    def __init__(self, cfg, scaler_dump, chpt, n_features, verbose=False):
        self.tensorflowConfig = cfg
        self.scalerDump = scaler_dump
        self.checkpoint = chpt
        self.nFeatures = n_features
        self.verbose = verbose
        # create tensorflow graph
        self.reloadModel()

    # read network architecture/hyper parameters from *.cfg file
    def loadModelConfig(self):
        with open(self.tensorflowConfig, 'r') as inputFile:
            lines = inputFile.read()
        return eval(lines)

    # rebuild tensorflow graph
    def reloadModel(self):
        self.parameters = self.loadModelConfig()

        # for evaluation of weights always set limitResources=True which limits threads (and therefore memory)
        self.clf = TensorflowDNNClassifier(parameters=self.parameters, nFeatures=self.nFeatures, limitResources=True)
        self.clf.buildModel()

        # restore rom checkpoint
        self.clf.restore(self.checkpoint)
        with open(self.scalerDump, "r") as inputFile:
            self.scaler = pickle.load(inputFile)

        # initialize arrays for transfering data to graph
        self.data = {'test':{
                                'X': np.full((1, self.clf.nFeatures), 0.0),
                                'y': np.full((1,), 1.0, dtype=np.float32),
                                'sample_weight': np.full((1,), 1.0, dtype=np.float32),
                            }
                    }

    # called from main loop for every event

    array=[1.123,1.123,1.123,1.123,1.123,1.123,1.123,1.123,1.123,1.123,1.123,1.123,1.123,1.123]

    def EvaluateDNN(self, array):

        for i in range(len(array)):
            if self.verbose:
                print array[i]
            self.data['test']['X'][0,i] = array[i] #FIXME: here you access the event vars

        # standardize
        self.data['test']['X_sc'] = self.scaler.transform(self.data['test']['X'])

        # evaluate
        probabilities = self.clf.session.run(self.clf.predictions, feed_dict={
                            self.clf.x: self.data['test']['X_sc'],
                            self.clf.y_: self.data['test']['y'].astype(int),
                            self.clf.w: self.data['test']['sample_weight'],
                            self.clf.keep_prob: np.array([1.0]*len(self.clf.pKeep), dtype=np.float32),
                            self.clf.learning_rate_adam: self.parameters['learning_rate_adam_start'],
                            self.clf.is_training: False,
                        })

        # fill output branches
        if self.verbose:
            print "DNN for this event is : ",  probabilities[0,0]  #FIXME: you want to write somewhere rather than printing

        return probabilities[0,0]

