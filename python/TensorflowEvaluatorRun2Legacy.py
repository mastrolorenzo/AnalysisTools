#!/usr/bin/env python
from __future__ import print_function
import numpy
import tensorflow as tf
import pickle
import os
import json
from MyStandardScaler import StandardScaler

# this is a standalone evaluator class which loads the graph from the proto-buff file instead of creating it
# with the tfDNNclassifier.py class and reloading the weights only.
# both ways for evaluation should lead to the same result!

class TensorflowDNNEvaluator(object):

    def __init__(self, checkpoint, scaler=None):
        self.checkpoint = checkpoint
        self.scaler = scaler
        if self.scaler:
            self.restoreScaler(self.scaler)
        infoFilePath = checkpoint + '.info'
        if os.path.isfile(infoFilePath):
            with open(infoFilePath) as infoFile:
                self.info = json.load(infoFile)
            #print(self.info)
        self.graph = tf.Graph()
        tf.reset_default_graph()
        self.session = tf.Session(graph=self.graph, config=tf.ConfigProto(device_count={"CPU":1}, inter_op_parallelism_threads=1, intra_op_parallelism_threads=1))
        with self.graph.as_default():
            self.session.run([tf.global_variables_initializer()])
            self.vars_from_ckpt(self.checkpoint)
            self.inputs = self.graph.get_tensor_by_name("input_data:0")
            self.outputs  = self.graph.get_tensor_by_name("Softmax:0")
            print("shapes:", self.inputs.shape, "--->", self.outputs.shape)

    def vars_from_ckpt(self, checkpoint):
        metafile = checkpoint + '.meta'
        self.saver = tf.train.import_meta_graph(metafile)
        self.saver.restore(self.session, checkpoint)
    
    def restoreScaler(self, fileName): 
        with open(fileName, "r") as inputFile:
            self.scaler = pickle.load(inputFile)
        return self.scaler

    def MultiEvaluateDNN(self, inputs):
        inputs_scaled = self.scale(inputs) if self.scaler else inputs
        with self.graph.as_default():
            result = self.session.run(self.outputs, feed_dict={self.inputs: inputs_scaled})
        return result

    def EvaluateDNN(self, oneSetOfInputs):
        oneSetOfInputs=[oneSetOfInputs]
        probabilities=self.MultiEvaluateDNN(oneSetOfInputs)
        return probabilities[0,0]

    def EvaluateMultiDNN(self, oneSetOfInputs):
        oneSetOfInputs=[oneSetOfInputs]
        probabilities=self.MultiEvaluateDNN(oneSetOfInputs)
        return probabilities[0]

    def scale(self, inputs):
        return self.scaler.transform(inputs)


