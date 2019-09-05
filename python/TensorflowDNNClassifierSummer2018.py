#!/usr/bin/env python
from __future__ import print_function
import tensorflow as tf
import random
import numpy as np
import math
import pickle
import os
import sys
from copy import deepcopy
import json
import gzip

# use MyStandardScaler if sklearn is not available
from MyStandardScaler import StandardScaler
#from sklearn.preprocessing import StandardScaler

# root is not needed for training, just for producing the control plots
try:
    import ROOT
except:
    print("warning: ROOT not found! some plotting features will not work")

# not available in some versions of sklearn!
try:
    from sklearn.preprocessing import QuantileTransformer
except:
    print("warning: no sklearn found or sklearn version does not support QuantileTransformer")

tf.logging.set_verbosity(tf.logging.INFO)

# ------------------------------------------
# mirror printout to file ("tee")
# ------------------------------------------
class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        try:
            self.log = open(filename, "a")
        except:
            self.log = open(filename, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

# ------------------------------------------
# DNN classifier
# ------------------------------------------
class TensorflowDNNClassifier(object):

    def __init__(self, project=None, input=None, parameters=None, nFeatures=None, limitResources=False):
        self.limitResources = limitResources
        self.parameters = {
            # epochs, data is reload
            # nStepsPerEpoch: repeat every epoch n times (without reloading data)
            # "real" number of epochs is product of both
            'nEpochs': 100,
            'nStepsPerEpoch': 50,

            # adam_epsilon is added in some divisions for numerical stability
            'adam_epsilon': 1.0e-9,

            # set total weight of signal to background, balanceClasses for multi class
            # (needed e.g. for cross_entropy loss function)
            'balanceSignalBackground': True,
            'balanceClasses': False,

            # for plotting: rescale MVA score [0, 1] to match the distribution of raw score in these quantiles
            'mvaScoreRescalingPercentileLow': 0.007,
            'mvaScoreRescalingPercentileHigh': 0.997,

            # list of feature numbers to remove from training, starting from 0
            'removeFeature': [],

            # batched training (saves a bit of memory)
            'batchSize': -1,

            # network architecture
            'nNodes': [64, 32, 16, 16],

            # dropout
            'pDropout': [0.3, 0.5, 0.5, 0.5],

            # multiply p_dropout by this factor every epoch
            'dropoutDecay': 1.0,

            # additional weights between not neighboring layers
            # example: {5: [1, 2]} add skip connections from layer 5 to 1 and 2
            'skipConnections': {},

            # random weight/bias initialization scale factors
            'wInitScale': 0.01,
            'bInitScale': 0.01,

            # L2 regularization term
            'regularization_strength': 0.0,

            # shuffle training samples for every epopch (needed together with batch)
            'shuffle': True,

            # use systematic uncertainties to smear training sample. scale up by factor
            # and add additional noise
            'systematics_scaling_factor': 1.0,
            'additional_noise': 0.00,

            # does not really help, to be removed
            'adaptiveRate': False,
            'rateGamma': 1.0,

            # scale layer outputs (inputs are already scaled before feeding to the network)
            # for later layers it can be given as a list of layer indices, starting from 1 for the first hidden layer, e.g.: [3,4]
            'batchNormalization': [1, 2, 3, 4],

            # loss function, possible values:
            #   cross_entropy:      standard cross entropy loss
            #   significance:       (smoothly) binned S/sqrt(S+B) estimate
            #   significance_qA     Asimov
            #   significance_qAS    Asimov with background uncertainty
            'loss': 'cross_entropy',

            # parameters for 'significance' loss function
            'signif_loss_smoothness': 500.0,
            'signif_loss_nbins': 15,

            # prevent loss function from optimizing for "empty" bins for sqrt(q_A)
            'signif_loss_low_b_threshold': 3.0,
            'signif_loss_low_b_threshold_width': 1.0,

            # approximated systematic uncertainty on background
            'signif_loss_sysApprox_constant': 1.5,
            'signif_loss_sysApprox_linear': 0.1,

            # number of toys to evaluate background variance for systematics given as per event weights
            'weight_sys_ntoys': -1,

            # constant which is added to the background variance computed by toys
            'signif_loss_sys_variance_offset': 0.1,

            # print message if value higher than this is found
            'best_expected_significance': 99999,
        }

        self.parameters.update(parameters)
        if nFeatures:
            self.nFeatures = nFeatures
        self.nClasses = 2

        # when training with significance loss function, don't reweight signal and background to have same weights
        if self.parameters['loss'].startswith('significance') and self.parameters['balanceSignalBackground']:
            self.parameters['balanceSignalBackground'] = False
            print("\x1b[97m\x1b[41mwarning: balanceSignalBackground was disabled due to significance loss function! \x1b[0m")

        self.projectName = project if project else "dnn"
        self.inputFileName = input

        self.systematics = ['JER_Up', 'JER_Down', 'PileUpDataMC_Up', 'PileUpDataMC_Down', 'PileUpPtRef_Up', 'PileUpPtRef_Down',
                       'PileUpPtBB_Up', 'PileUpPtBB_Down', 'PileUpPtEC1_Up', 'PileUpPtEC1_Down', 'RelativeJEREC1_Up',
                       'RelativeJEREC1_Down', 'RelativeFSR_Up', 'RelativeFSR_Down', 'RelativeStatFSR_Up',
                       'RelativeStatFSR_Down',
                       'RelativeStatEC_Up', 'RelativeStatEC_Down', 'RelativePtBB_Up', 'RelativePtBB_Down',
                       'RelativePtEC1_Up',
                       'RelativePtEC1_Down', 'AbsoluteScale_Up', 'AbsoluteScale_Down', 'AbsoluteMPFBias_Up',
                       'AbsoluteMPFBias_Down',
                       'AbsoluteStat_Up', 'AbsoluteStat_Down', 'SinglePionECAL_Up', 'SinglePionECAL_Down',
                       'SinglePionHCAL_Up',
                       'SinglePionHCAL_Down', 'Fragmentation_Up', 'Fragmentation_Down']
        self.systematics_withoutUD = list(set(['_'.join(syst.split('_')[:-1]) for syst in self.systematics]))

        if input:
            self.outputBaseDir = "results/"
            self.outputName = "rnd_"
            self.outputNumber = 1

            self.directories = {'output': self.getOutputFolder()}
            self.setOutputFolders(self.directories['output'])

            print("saving model to \x1b[34m", self.directories['output'], "\x1b[0m")
            try:
                os.makedirs(self.directories['output'])
            except:
                pass

            self.logfileName = self.initLogger()
            print("saving logfile to \x1b[34m", self.logfileName, "\x1b[0m")

            print("-" * 50)
            print("start tensorboard with:")
            print("   \x1b[34mtensorboard  --logdir %s\x1b[0m" % (os.environ['PWD'] + '/' + self.directories['output']))
            print("-" * 50)
        else:
            self.directories = {'output': 'eval_output'}
            self.setOutputFolders(self.directories['output'])

    # ------------------------------------------
    # set output folder
    # ------------------------------------------
    def setOutputFolders(self, folder):
        self.directories['output'] = folder
        self.directories['checkpoint'] = self.directories['output'] + '/checkpoints/model.ckpt'
        self.directories['checkpoint_best'] = self.directories['output'] + '/checkpoints/best/model.ckpt'

    # ------------------------------------------
    # create folder name from network parameters
    # ------------------------------------------
    def getOutputFolder(self):
        # create always a new output directory
        nodesStr = []
        dropoutStr = []
        for i in range(len(self.parameters['nNodes'])):
            if self.parameters['nNodes'][i] < 1:
                break
            nodesStr += ["%d" % self.parameters['nNodes'][i]]
            dropoutStr += ["%.2f" % self.parameters['pDropout'][i]]
        inputFileName = self.inputFileName.split('/')[-1]
        if len(nodesStr) < 1:
            nodesStr = ["reg"]
        output_dir = self.outputBaseDir + '/%s/%s/%s/%s/%1.3e/%s%d' % (
            self.projectName, inputFileName, '-'.join(nodesStr), '-'.join(dropoutStr),  self.parameters['learning_rate_adam_start'], self.outputName,
            self.outputNumber)
        while os.path.isdir(output_dir) and self.outputNumber < 9999:
            self.outputNumber += 1
            output_dir = self.outputBaseDir + '/%s/%s/%s/%s/%1.3e/%s%d' % (
                self.projectName, inputFileName, '-'.join(nodesStr), '-'.join(dropoutStr), self.parameters['learning_rate_adam_start'],
                self.outputName, self.outputNumber)
        return output_dir

    # ------------------------------------------
    # initLogger
    # ------------------------------------------
    def initLogger(self):
        logfileName = self.directories['output'] + "/output.txt"
        sys.stdout = Logger(logfileName)
        return logfileName

    def unpickle(self, fileName):
        if '.dmpz' in fileName:
            inputFile = gzip.open(fileName,'rb')
        else:
            inputFile = open(fileName, "rb")
        try:
            data = pickle.load(inputFile, encoding='latin1')
        except:
            data = pickle.load(inputFile)
        inputFile.close()
        return data

    # ------------------------------------------
    # load pickled numpy arrays as input data
    # ------------------------------------------
    def loadInput(self):
        # read input files, which contain pickled numpy arrays
        self.data_original = self.unpickle(self.inputFileName)
        self.data = deepcopy(self.data_original)

        self.data_original['test'] = None


        if 'limitTestSamples' in self.parameters and self.parameters['limitTestSamples'] > 0:

            # take a random subset with maximum length
            shuffleIndices = list(range(0, len(self.data['test']['y'])))
            random.shuffle(shuffleIndices)
            reweight = 1.0*len(self.data['test']['y']) / min(len(self.data['test']['y']),self.parameters['limitTestSamples'])
            shuffleIndices = shuffleIndices[:min(len(self.data['test']['y']),self.parameters['limitTestSamples'])]

            self.data['test']['X'] = self.data['test']['X'][shuffleIndices]
            self.data['test']['y'] = self.data['test']['y'][shuffleIndices]
            self.data['test']['sample_weight'] = self.data['test']['sample_weight'][shuffleIndices]

            # reweight samples
            for i in range(len(self.data['test']['sample_weight'])):
                self.data['test']['sample_weight'][i] *= reweight

        #load systematics list from data
        if 'systematics' in self.data['meta']:
            self.systematics = self.data['meta']['systematics']
            self.systematics_withoutUD = list(set(['_'.join(syst.split('_')[:-1]) for syst in self.systematics]))

        # remove columns
        variableNames = self.data['meta']['variables'].split(' ')
        for removeFeature in sorted(self.parameters['removeFeature'], reverse=True):
            self.data['train']['X'] = np.delete(self.data['train']['X'], removeFeature, 1)
            for syst in self.systematics:
                self.data['train']['X_'+syst] = np.delete(self.data['train']['X_'+syst], removeFeature, 1)
            self.data['test']['X'] = np.delete(self.data['test']['X'], removeFeature, 1)
            self.data_original['train']['X'] = np.delete(self.data_original['train']['X'], removeFeature, 1)
            del variableNames[removeFeature]
        self.data['meta']['variables'] = ' '.join(variableNames)

        self.nFeatures = self.data['train']['X'].shape[1]
        print("nFeatures = ", self.nFeatures)

        if 'category_labels' not in self.data:
            self.data['category_labels'] = {0: 'BKG', 1: 'SIG'}

        # count test set
        counts = [0] * len(self.data['category_labels'])
        for s in self.data['test']['y']:
            counts[int(s)] += 1
        print("test set:")
        for i,k in enumerate(counts):
            print(self.data['category_labels'][i],":",k)

        # count training set
        counts = [0] * len(self.data['category_labels'])
        for s in self.data['train']['y']:
            counts[int(s)] += 1
        print("train set:")
        for i,k in enumerate(counts):
            print(self.data['category_labels'][i],":",k)

        # multiclass handling
        self.classLabels = self.data['category_labels']
        self.nClasses = len(list(self.classLabels.keys()))
        self.signalClassIDs = [x for x, y in self.classLabels.items() if 'SIG_' in y.upper() or y.upper() == 'SIG']
        if len(self.signalClassIDs) < 1 or len(self.signalClassIDs) >= self.nClasses:
            print("ERROR: no signal or no background defined!")
            exit()
        print("list of classes: (signals in \x1b[32mgreen\x1b[0m, backgrounds in \x1b[31mred\x1b[0m)")
        for classID, classLabel in self.classLabels.items():
            if classID in self.signalClassIDs:
                print("\x1b[32mclass", classID, "=>", classLabel, "\x1b[0m")
            else:
                print("\x1b[31mclass", classID, "=>", classLabel, "\x1b[0m")

        # for multiclass problems precompute additional vector to indicate if signal or background for loss computation
        self.data['train']['c'] = deepcopy(self.data['train']['y'])
        for i in range(len(self.data['train']['c'])):
            self.data['train']['c'][i] = 1 if int(self.data['train']['y'][i]) in self.signalClassIDs else 0

        self.data_original['train']['c'] = deepcopy(self.data_original['train']['y'])
        for i in range(len(self.data_original['train']['c'])):
            self.data_original['train']['c'][i] = 1 if int(self.data_original['train']['y'][i]) in self.signalClassIDs else 0

        self.data['test']['c'] = deepcopy(self.data['test']['y'])
        for i in range(len(self.data['test']['c'])):
            self.data['test']['c'][i] = 1 if int(self.data['test']['y'][i]) in self.signalClassIDs else 0

        # systematics applied as weights
        self.sysWeights_withoutUD = []
        if 'meta' in self.data:
            print("train-cut:", self.data['meta']['trainCut'])
            print("test-cut:", self.data['meta']['testCut'])
            if 'weightSYS' in self.data['meta']:
                print("systematics applied as per event weights:")
                for w in self.data['meta']['weightSYS']:
                     self.data_original['sample_weight_' + w] = None
                     print(" -> ", w)
                     wR = w.replace('Up','').replace('Down','')
                     if wR not in self.sysWeights_withoutUD:
                         self.sysWeights_withoutUD.append(wR)
                print("list w/o U/D:", self.sysWeights_withoutUD)

    # ------------------------------------------
    # feature scaling ----> mean 0.0 stddev 1.0
    # ------------------------------------------
    def scaleFeatures(self):
        self.scaler = StandardScaler()
        #self.scaler = QuantileTransformer()
        variableNames = self.data['meta']['variables'].split(' ')
        print('weights:')
        print('train', ' '.join(['%r'%x for x in list(self.data['train']['sample_weight'][0:10])]))
        print('test ', ' '.join(['%r'%x for x in list(self.data['test']['sample_weight'][0:10])]))
        for i,variableName in enumerate(variableNames):
            if len(variableName) > 49:
                variableName = variableName[:46] + '...'
            print(variableName.ljust(50), 'train', ' '.join(['%r'%x for x in list(self.data['train']['X'][0:5,i])]))
            print(variableName.ljust(50), 'test ', ' '.join(['%r'%x for x in list(self.data['test']['X'][0:5,i])]))

        # load scaler parameters file if given, e.g. for evaluation on a different data set or continuing training on a
        # different dataset
        if 'scaler' in self.parameters:
            print("scaler loaded from \x1b[34m", self.parameters['scaler'], " \x1b[0m")
            with open(self.parameters['scaler'], "r") as inputFile:
                self.scaler = pickle.load(inputFile)
            self.data['train']['X'] = self.scaler.transform(self.data['train']['X'])
        else:
            # initialize a new scaler
            self.data['train']['X'] = self.scaler.fit_transform(self.data['train']['X'])

        # apply scaling of input features
        if self.parameters['systematics_scaling_factor'] > 0.0:
            for syst in self.systematics:
                self.data['train']['X_' + syst] = self.scaler.transform(self.data['train']['X_' + syst])
        self.data['test']['X'] = self.scaler.transform(self.data['test']['X'])
        self.data_original['train']['X'] = self.scaler.transform(self.data_original['train']['X'])

        try:
            print("scaler params:", self.scaler, self.scaler.get_params(), self.scaler.scale_, self.scaler.mean_)
        except:
            print("using a different scaler than StandardScaler()!")

        ### write scaler parameters to disk
        ##scalerFileName = self.directories['output'] + "/scaler.dmp"
        ##with open(scalerFileName, "wb") as outputFile:
        ##    pickle.dump(self.scaler, outputFile, protocol=2)
        ##print("scaler dumped to \x1b[34m", scalerFileName, " \x1b[0m")

        # print list of scaled variables
        print("Variable".ljust(50),"Mean".ljust(10), "stddev".ljust(10))
        variableNames = self.data['meta']['variables'].split(' ')
        for i,variableName in enumerate(variableNames):
            if len(variableName) > 49:
                variableName = variableName[:46] + '...'
            print(variableName.ljust(50), ("%1.3f"%self.scaler.mean_[i]).ljust(10), ("%1.3f"%self.scaler.scale_[i]).ljust(10))

    # ------------------------------------------
    # further pre-processing, e.g. class weights
    # ------------------------------------------
    def preprocess(self):
        # balance the total normalization of signal/background
        if self.parameters['balanceSignalBackground']:
            totalWeight = {0: 0.0, 1: 0.0}
            for i in range(len(self.data['train']['sample_weight'])):
                totalWeight[1 if int(self.data['train']['y'][i]) in self.signalClassIDs else 0] += self.data['train']['sample_weight'][
                    i]
            totalWeightAllClasses = sum([y for x, y in totalWeight.items()])
            reweight = {0: 1.0, 1: 1.0}
            for i in [0, 1]:
                reweight[i] = totalWeightAllClasses / totalWeight[i] if totalWeight[i] > 0.0 else 1.0
                print("balancing classes, reweight group", i, "by", reweight[i])
            reweight[1] *= self.parameters['addWeightSignal']

            # apply balancing, original weights are kept in data_original
            for i in range(len(self.data['train']['sample_weight'])):
                self.data['train']['sample_weight'][i] *= reweight[1 if int(self.data['train']['y'][i]) in self.signalClassIDs else 0]

        # balance the total normalization of each class (multiclass output!)
        # will yield to imbalanced signal/background if different numbers of signal/background classes exist
        elif self.parameters['balanceClasses']:
            totalWeight = {x: 0.0 for x in list(self.classLabels.keys())}
            for i in range(len(self.data['train']['sample_weight'])):
                totalWeight[int(self.data['train']['y'][i])] += self.data['train']['sample_weight'][i]
            totalWeightAllClasses = sum([y for x, y in totalWeight.items()])
            reweight = {x: 1.0 for x in list(self.classLabels.keys())}
            for classID in list(self.classLabels.keys()):
                reweight[classID] = totalWeightAllClasses / totalWeight[classID] if totalWeight[classID] > 0.0 else 1.0
                print("balancing classes, reweight", self.classLabels[classID], "by", reweight[classID])

            # apply balancing, original weights are kept in data_original
            for i in range(len(self.data['train']['sample_weight'])):
                self.data['train']['sample_weight'][i] *= reweight[self.data['train']['y'][i]]

        print("shape train:", self.data['train']['X'].shape)
        print("shape test: ", self.data['test']['X'].shape)
        self.truth = self.data['test']['y'].astype(int)
        self.truth_train = self.data_original['train']['y'].astype(int)

    # ------------------------------------------
    # TensorFlow Implementation
    # ------------------------------------------
    def buildModel(self):
        print("building tensorflow graph with parameters")
        for k,v in self.parameters.items():
            print(" "+k.ljust(40)+" "+"%r"%v)
        print("initialize session...")
        if self.limitResources:
            self.session = tf.Session(config=tf.ConfigProto(device_count={"CPU":1}, inter_op_parallelism_threads=1, intra_op_parallelism_threads=1))
        else:
            self.session = tf.Session()
        print("initialized session!")

        # send expected significance for training and test set to tensorboard
        self.tensorboardVar = tf.Variable(0, "tensorboardVar", dtype=tf.float32)
        self.tensorboardVar_loss = tf.Variable(0, "tensorboardVar_loss", dtype=tf.float32)
        self.pythonVar = tf.placeholder("float32", [], name="expected_signif")
        self.pythonVar_loss = tf.placeholder("float32", [], name="loss_value")
        self.update_tensorboardVar = self.tensorboardVar.assign(self.pythonVar)
        self.update_tensorboardVar_loss = self.tensorboardVar_loss.assign(self.pythonVar_loss)
        tf.summary.scalar("expected_significance", self.update_tensorboardVar)
        tf.summary.scalar("loss", self.update_tensorboardVar_loss)
        print("....>>>>")

        self.merged = tf.summary.merge_all()
        
        # CP don't need this output
        #self.sum_writer_test = tf.summary.FileWriter(self.directories['output'] + '/plot_signif_test/')
        #self.sum_writer_train = tf.summary.FileWriter(self.directories['output'] + '/plot_signif_train/')
        #self.sum_writer_hist = tf.summary.FileWriter(self.directories['output'] + '/plot_signif_hist/')

        # to be implemented: fill global_step variable to be able to continue from
        # snapshot and have right step number in tensorboard
        self.global_step_tensor = tf.Variable(0, trainable=False, name='global_step')

        # helper functions to initialize weights for layers
        def weight_variable(shape, ws=0.1, name=None):
            #initial = tf.truncated_normal(shape, stddev=0.1)
            #initial = tf.sqrt(2.0/shape[0]) * tf.random_uniform(shape, minval=0, maxval=shape[0]) #arxiv:1502.01852
            initial = tf.truncated_normal(shape, stddev=math.sqrt(2.0/(shape[0]+shape[1]))) #glorot
            return tf.Variable(initial, name=name)

        # always use a small positive bias for initialization!
        def bias_variable(shape, ws=0.1, name=None):
            initial = tf.constant(ws, shape=shape)
            return tf.Variable(initial, name=name)

        # ------------------------------------------
        # INPUT variables
        # ------------------------------------------
        self.x = tf.placeholder(tf.float32, shape=[None, self.nFeatures], name='input_data')
        self.y_ = tf.placeholder(tf.int32, shape=[None], name='labels')
        self.c_ = tf.placeholder(tf.int32, shape=[None], name='class_labels')
        self.w = tf.placeholder(tf.float32, shape=[None], name='sample_weights')
        self.w_toys = tf.placeholder(tf.float32, shape=[None, None], name='sample_weights_toys')

        # this allows multi-class input files be used for normal signal/background training
        self.cl = tf.placeholder(tf.float32, shape=[None, 2], name='classes')

        # ------------------------------------------
        # DNN implementation
        # ------------------------------------------
        self.learning_rate_adam = tf.placeholder(tf.float32, shape=[], name="learning_rate_adam")

        self.is_training = tf.placeholder(tf.bool, name="is_training")
        self.keep_prob = tf.placeholder(tf.float32, name="keep_prob")
        self.pKeep = np.array([1.0 - pDrop for pDrop in self.parameters['pDropout']], dtype=np.float32)
        self.h_fc_bn = [self.x]
        self.h_fc_drop = [self.x]
        self.W_fc = [None]
        self.b_fc = [None]
        self.h_fc = [None]
        self.shortConnections = []

        # contruct DNN layer by layer
        lastLayer = 0
        print("add layers...")
        for lIdx,k in enumerate(self.parameters['nNodes'], 1):
            if self.parameters['nNodes'][lIdx-1] > 0:

                # add a new layer with weights and bias nodes
                weightScale = self.parameters['wInitScale'] / self.parameters['nNodes'][lIdx-1]
                biasScale = self.parameters['bInitScale'] / self.parameters['nNodes'][lIdx-1]
                self.W_fc.append(weight_variable([self.parameters['nNodes'][lIdx-2] if lIdx >=2 else self.nFeatures, self.parameters['nNodes'][lIdx-1]], ws=weightScale, name=("hidden_%d"%lIdx) if lIdx == 1 else None))
                self.b_fc.append(bias_variable([self.parameters['nNodes'][lIdx-1]], ws=biasScale))

                # activation function
                activationInput = tf.matmul(self.h_fc_bn[-1], self.W_fc[-1]) + self.b_fc[-1]
                if 'skipConnections' in self.parameters and lIdx in self.parameters['skipConnections']:
                    for fromLayerIdx in self.parameters['skipConnections'][lIdx]:
                        scWeights = weight_variable([self.parameters['nNodes'][fromLayerIdx-1],  self.parameters['nNodes'][lIdx-1]], ws=weightScale, name='skipConn')
                        scBias = bias_variable([self.parameters['nNodes'][lIdx-1]], ws=biasScale)
                        activationInput = activationInput + tf.matmul(self.h_fc_bn[fromLayerIdx], scWeights) + scBias
                self.h_fc.append(tf.nn.leaky_relu(activationInput, name="relu_%d"%lIdx))

                # dropout for regularization
                self.h_fc_drop.append(tf.nn.dropout(self.h_fc[-1], self.keep_prob[lIdx-1], name='dropout_%d'%lIdx))

                # normalize layer outputs (needs a lot of memory)
                if self.parameters['batchNormalization'] and lIdx in self.parameters['batchNormalization']:
                    self.h_fc_bn.append(tf.layers.batch_normalization(inputs=self.h_fc_drop[-1], axis=-1, momentum=0.99, epsilon=0.001, center=False, scale=False, training=self.is_training, name="batch_normalization_%d"%lIdx))
                else:
                    self.h_fc_bn.append(self.h_fc_drop[-1])

                lastLayer = lIdx
            else:
                break


        # output
        if lastLayer > 0:
            self.W_out = weight_variable([self.parameters['nNodes'][lastLayer-1], self.nClasses])
            self.b_out = bias_variable([self.nClasses])
            self.y = tf.matmul(self.h_fc_bn[-1], self.W_out) + self.b_out
        else:
            self.W_out = weight_variable([self.nFeatures, self.nClasses], name='output_weights')
            self.b_out = bias_variable([self.nClasses], name='output_bias')
            self.y = tf.matmul(self.x, self.W_out) + self.b_out

        # L2 regularization
        if self.parameters['regularization_strength'] > 0:
            self.regularizer = tf.contrib.layers.l2_regularizer(scale=self.parameters['regularization_strength'], scope=None)
            self.allweights = tf.trainable_variables()
            self.regularization_penalty = tf.contrib.layers.apply_regularization(self.regularizer, self.allweights)

        # probabilities for each of the two classes
        self.predictions = tf.nn.softmax(self.y, name='label_predictions')

        # cross-entropy loss function
        self.cross_entropy = tf.reduce_mean(self.w * tf.nn.sparse_softmax_cross_entropy_with_logits(logits=self.y, labels=self.y_), name='cross_entropy')

        # LOSS function
        if self.parameters['loss'].startswith('significance'):
            # use approximated expected significance as loss function
            self.histMin = tf.reduce_min(self.predictions[:,0])
            self.histMax = tf.reduce_max(self.predictions[:,0])
            self.histR = self.histMax - self.histMin + 1.0e-12
            # bin width
            self.histDR = 1.0 / self.parameters['signif_loss_nbins']
            self.predictions_S = ((1.0 - tf.cast(self.y_, dtype=tf.float32)) * self.predictions[:,0] - self.histMin) / self.histR
            self.hist_S_bin = [None] * self.parameters['signif_loss_nbins']
            self.predictions_B = ((tf.cast(self.y_, dtype=tf.float32)) * self.predictions[:,0] - self.histMin) / self.histR
            self.hist_B_bin = [None] * self.parameters['signif_loss_nbins']
            self.hist_SSB_bin = [None] * self.parameters['signif_loss_nbins']

            # calculate S and B per "smooth bin"
            self.hist_S_bin[0] = tf.constant(0.0)
            self.hist_B_bin[0] = tf.constant(0.1)
            for iBin in range(1,self.parameters['signif_loss_nbins']-1):
                self.hist_S_bin[iBin] = tf.reduce_sum(self.w * (0.5 - 0.5 * tf.erf(self.parameters['signif_loss_smoothness'] * (self.predictions_S - (iBin) * self.histDR) ) * tf.erf(self.parameters['signif_loss_smoothness'] * (self.predictions_S - (iBin+1) * self.histDR))))
                self.hist_B_bin[iBin] = tf.reduce_sum(self.w * (0.5 - 0.5 * tf.erf(self.parameters['signif_loss_smoothness'] * (self.predictions_B - (iBin) * self.histDR) ) * tf.erf(self.parameters['signif_loss_smoothness'] * (self.predictions_B - (iBin+1) * self.histDR))))
            self.hist_S_bin[self.parameters['signif_loss_nbins']-1] = tf.reduce_sum(self.w * (0.5 - 0.5 * tf.erf(self.parameters['signif_loss_smoothness'] * (self.predictions_S - (self.parameters['signif_loss_nbins']-1.0) * self.histDR) ) * tf.erf(500.0 * (self.predictions_S - (self.parameters['signif_loss_nbins']+0.01) * self.histDR))))
            self.hist_B_bin[self.parameters['signif_loss_nbins']-1] = tf.reduce_sum(self.w * (0.5 - 0.5 * tf.erf(self.parameters['signif_loss_smoothness'] * (self.predictions_B - (self.parameters['signif_loss_nbins']-1.0) * self.histDR) ) * tf.erf(500.0 * (self.predictions_B - (self.parameters['signif_loss_nbins']+0.01) * self.histDR))))

            # S^2/(S+B) per bin
            if self.parameters['loss'] == 'significance_qA':
                for iBin in range(self.parameters['signif_loss_nbins']):
                    # require minimum number of events per bin, otherwise decrease this bin's significance
                    minSBperBin = 0.5*(1.0 + tf.erf((1.0/self.parameters['signif_loss_low_b_threshold_width']) * (self.hist_B_bin[iBin] - self.parameters['signif_loss_low_b_threshold'])))
                    # asimov expected significance
                    self.hist_SSB_bin[iBin] = minSBperBin * 2.0 * ((self.hist_S_bin[iBin] + self.hist_B_bin[iBin]) * tf.log(1.0 + (self.hist_S_bin[iBin] / (self.hist_B_bin[iBin] + 0.001))) - self.hist_S_bin[iBin])

            elif self.parameters['loss'] == 'significance_qASapprox':
                for iBin in range(self.parameters['signif_loss_nbins']):
                    # require minimum number of events per bin, otherwise decrease this bin's significance
                    minSBperBin = 0.5*(1.0 + tf.erf((1.0/self.parameters['signif_loss_low_b_threshold_width']) * (self.hist_B_bin[iBin] - self.parameters['signif_loss_low_b_threshold'])))
                    # asimov expected significance with some approximate uncertainty estimate for background
                    sigma_b = tf.sqrt(self.hist_B_bin[iBin]) + 0.1*self.hist_B_bin[iBin] + 1.5
                    s = self.hist_S_bin[iBin]
                    b = self.hist_B_bin[iBin]
                    self.hist_SSB_bin[iBin] = minSBperBin * 2.0 * ((s+b)*tf.log(((s+b)*(b+sigma_b*sigma_b))/(b*b + (s+b)*sigma_b*sigma_b + 0.0001)) - (b*b)/(sigma_b*sigma_b)*tf.log(1.0+sigma_b*sigma_b*s/(b*(b+sigma_b*sigma_b))))

            elif self.parameters['loss'] == 'significance_qAS':

                nBins = self.parameters['signif_loss_nbins']
                nToys = self.parameters['weight_sys_ntoys']
                alpha = self.parameters['signif_loss_smoothness']

                # compute background prediction for all the toys
                # TODO: this could need some optimization...
                self.hist_B_bin_toy = [[tf.constant(0.0)] * nToys]
                for iBin in range(1, nBins-1):
                    self.hist_B_bin_toy.append([None] * nToys)
                    for iToy in range(nToys):
                        self.hist_B_bin_toy[iBin][iToy] = tf.reduce_sum(self.w_toys[iToy] * (0.5 - 0.5 * tf.erf(alpha * (self.predictions_B - (iBin) * self.histDR) ) * tf.erf(alpha * (self.predictions_B - (iBin+1) * self.histDR))))

                self.hist_B_bin_toy.append([None] * nToys)
                for iToy in range(nToys):
                    self.hist_B_bin_toy[nBins-1][iToy] = tf.reduce_sum(self.w_toys[iToy] * (0.5 - 0.5 * tf.erf(alpha * (self.predictions_B - (nBins-1.0) * self.histDR) ) * tf.erf(500.0 * (self.predictions_B - (nBins+0.1) * self.histDR))))

                # compute variance of predicted background events
                self.hist_sigmaB_bin = []
                for iBin in range(nBins):
                    variance = sum([ (self.hist_B_bin_toy[iBin][iToy]-self.hist_B_bin[iBin])*(self.hist_B_bin_toy[iBin][iToy]-self.hist_B_bin[iBin]) for iToy in range(nToys)])/(1.0*nToys)
                    self.hist_sigmaB_bin.append(variance)

                for iBin in range(nBins):
                    # require minimum number of events per bin, otherwise decrease this bin's significance
                    minSBperBin = 0.5*(1.0 + tf.erf((1.0/self.parameters['signif_loss_low_b_threshold_width']) * (self.hist_B_bin[iBin] - self.parameters['signif_loss_low_b_threshold'])))

                    s = self.hist_S_bin[iBin]
                    b = self.hist_B_bin[iBin]
                    var_b = self.parameters['signif_loss_sys_variance_offset'] + self.hist_sigmaB_bin[iBin]

                    # asimov expected significance with some approximate uncertainty estimate for background
                    self.hist_SSB_bin[iBin] = minSBperBin * 2.0 * ((s+b)*tf.log(((s+b)*(b+var_b))/(b*b + (s+b)*var_b + 0.01)) - (b*b)/(var_b + 0.01)*tf.log(1.0+var_b*s/(b*(b+var_b))))

            elif self.parameters['loss'] == 'significance_sysApprox':
                for iBin in range(self.parameters['signif_loss_nbins']):
                    #self.hist_SSB_bin[iBin] = self.hist_S_bin[iBin] * self.hist_S_bin[iBin] / tf.square(self.parameters['signif_loss_sysApprox_constant'] + tf.sqrt(tf.abs(self.hist_B_bin[iBin])) + tf.abs(self.parameters['signif_loss_sysApprox_linear']*self.hist_B_bin[iBin]))
                    self.hist_SSB_bin[iBin] = self.hist_S_bin[iBin] * self.hist_S_bin[iBin] / tf.square(1.0 + tf.sqrt(tf.abs(self.hist_B_bin[iBin])))

            else:
                for iBin in range(self.parameters['signif_loss_nbins']):
                    self.hist_SSB_bin[iBin] = self.hist_S_bin[iBin] * self.hist_S_bin[iBin] / (self.hist_S_bin[iBin] + self.hist_B_bin[iBin] + 1.0e-8)

            # total S/sqrt(S+B)
            self.hist_SSB = tf.sqrt(sum(self.hist_SSB_bin))

            self.total_loss = - self.hist_SSB
            #+ tf.cond(self.hist_B_bin[1]<self.hist_B_bin[self.parameters['signif_loss_nbins']-2], lambda: self.hist_B_bin[self.parameters['signif_loss_nbins']-2] - self.hist_B_bin[1], lambda: tf.constant(0.0))
        else:
            # cross entropy
            if self.parameters['regularization_strength'] > 0:
                self.total_loss = self.cross_entropy + self.regularization_penalty
            else:
                self.total_loss = self.cross_entropy

        self.update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
        with tf.control_dependencies(self.update_ops):
            self.train_step = tf.train.AdamOptimizer(
                learning_rate=self.parameters['learning_rate_adam_start'],
                epsilon=self.parameters['adam_epsilon'],
            ).minimize(self.total_loss)

        ## CP don't need this
        #file_writer = tf.summary.FileWriter(self.directories['output'], self.session.graph)

        # initialize variables
        self.session.run(tf.global_variables_initializer())
        self.saver = tf.train.Saver(restore_sequentially=True)

        # copy this script to every output folder!
        # shutil.copy('train.py', self.directories['output'] + '/train.py')

    # ------------------------------------------
    # load checkpoint from file
    # ------------------------------------------
    def restore(self, savePath):
        print("going to restore from \x1b[34m",savePath,"\x1b[0m")
        try:
            self.saver.restore(self.session, savePath)
            return True
        except Exception as e:
            print("EXCEPTION:",e)
            print("\x1b[31merror: could not restore checkpoint \x1b[0m")
            return False

    # ------------------------------------------
    # train the classifier
    # ------------------------------------------
    def fit(self):
        expected_significance_test = 0
        expected_significance_train = 0
        recordForEarlyStop = 0
        recordForEarlyStopTime = 0
        firstCheckpoint = True
        learningRate = self.parameters['learning_rate_adam_start']

        # make all backgrounds 0 and all signals 1
        classToSB0 = []
        for i in range(self.nClasses):
            if i in self.signalClassIDs:
                classToSB0.append([0.0, 1.0])
            else:
                classToSB0.append([1.0, 0.0])
        classToSB = np.array(classToSB0, dtype=np.float32)

        logfile = open(self.directories['output'] + '/log.txt', 'w')

        nSamples = len(self.data['train']['y'])

        # epochs here consist of several 'subepochs'
        for i in range(self.parameters['nEpochs']):

            # shuffle samples every few epochs
            if self.parameters['shuffle']:
                shuffleIndices = list(range(0, nSamples))
                random.shuffle(shuffleIndices)
                self.data['train']['X'] = self.data['train']['X'][shuffleIndices]
                if self.parameters['systematics_scaling_factor'] > 0.0:
                    for syst in self.systematics:
                        self.data['train']['X_' + syst] = self.data['train']['X_' + syst][shuffleIndices]
                    if 'weightSYS' in self.data['meta']:
                        for syst in self.data['meta']['weightSYS']:
                            self.data['train']['sample_weight_' + syst] = self.data['train']['sample_weight_' + syst][shuffleIndices]
                self.data['train']['y'] = self.data['train']['y'][shuffleIndices]
                self.data['train']['c'] = self.data['train']['c'][shuffleIndices]
                self.data['train']['sample_weight'] = self.data['train']['sample_weight'][shuffleIndices]

            # draw a new random sample for the inputs
            self.data['train']['X_SYS'] = deepcopy(self.data['train']['X'])
            if self.parameters['systematics_scaling_factor'] > 0:
                for syst in self.systematics_withoutUD:
                    v_gauss = np.random.normal(0, self.parameters['systematics_scaling_factor'], [nSamples])
                    v_p = np.array([x1 if x1 > 0.0 else 0.0 for x1 in v_gauss], dtype=np.float32)
                    v_n = np.array([x1 if x1 < 0.0 else 0.0 for x1 in v_gauss], dtype=np.float32)
                    self.data['train']['X_SYS'] = self.data['train']['X_SYS'] + np.multiply((self.data['train']['X_' + syst + '_Up']-self.data['train']['X']), v_p[:, np.newaxis]) + np.multiply((self.data['train']['X_' + syst + '_Down']-self.data['train']['X']), v_n[:, np.newaxis])
            if self.parameters['additional_noise'] > 0.0:
                v_noise = np.random.normal(0, self.parameters['additional_noise'], [nSamples, self.nFeatures])
                self.data['train']['X_SYS'] = self.data['train']['X_SYS'] + v_noise
            batchSize = self.parameters['batchSize']
            if batchSize < 0:
                batchSize = len(self.data['train']['X_SYS'])

            # draw a new random sample for the weightS
            #self.data['train']['sample_weight_SYS'] = deepcopy(self.data['train']['sample_weight'])
            #if self.parameters['systematics_scaling_factor'] > 0:
            #    for syst in self.sysWeights_withoutUD:
            #        w_gauss = np.random.normal(0, self.parameters['systematics_scaling_factor'], [nSamples])
            #        w_p = np.array([x1 if x1 > 0.0 else 0.0 for x1 in w_gauss], dtype=np.float32)
            #        w_n = np.array([x1 if x1 < 0.0 else 0.0 for x1 in w_gauss], dtype=np.float32)
            #        self.data['train']['sample_weight_SYS'] = self.data['train']['sample_weight_SYS'] + np.multiply((self.data['train']['sample_weight_' + syst + 'Up']-self.data['train']['sample_weight']), w_p) + np.multiply((self.data['train']['sample_weight_' + syst + 'Down']-self.data['train']['sample_weight']), w_n)
            self.data['train']['sample_weight_SYS'] = self.data['train']['sample_weight']

            # create toys for systematic variations (weight systematics)
            #print("weights before:")
            #print(self.data['train']['sample_weight_toys'])

            self.data['meta']['weightSYS']
            for dataset in ['train', 'test']:
                self.data[dataset]['sample_weight_toys'] = np.vstack([self.data[dataset]['sample_weight']] * (self.parameters['weight_sys_ntoys'] if self.parameters['weight_sys_ntoys'] > 0 else 1))
                if dataset == 'train':
                    if self.parameters['weight_sys_ntoys']>0:
                        print("generating toys...")
                        for j in range(self.parameters['weight_sys_ntoys']):
                            for syst in self.sysWeights_withoutUD:
                                # systematics given as two half-sided gaussians
                                w_gauss = np.random.normal(0, self.parameters['systematics_scaling_factor'], [nSamples])
                                w_p = np.array([x1 if x1 > 0.0 else 0.0 for x1 in w_gauss], dtype=np.float32)
                                w_n = np.array([x1 if x1 < 0.0 else 0.0 for x1 in w_gauss], dtype=np.float32)
                                # add a random variation from this systematic uncertainty
                                self.data[dataset]['sample_weight_toys'][j] = (self.data[dataset]['sample_weight_toys'][j] +
                                    np.multiply((self.data[dataset]['sample_weight_' + syst + 'Up']-self.data[dataset]['sample_weight']), w_p) +
                                    np.multiply((self.data[dataset]['sample_weight_' + syst + 'Down']-self.data[dataset]['sample_weight']), w_n))
                        print("generating toys...done. ->", self.parameters['weight_sys_ntoys'], "toys generated")

            #print("weights after:")
            #print(self.data['train']['sample_weight_toys'])
            numTrainingEventsTotal = (i+1)*(self.parameters['nStepsPerEpoch'])

            # run a few steps
            for j in range(self.parameters['nStepsPerEpoch']):
                self.train_step.run(feed_dict={
                                self.x: self.data['train']['X_SYS'][:batchSize],
                                self.y_: self.data['train']['y'][:batchSize],
                                #self.c_: self.data['train']['c'][:batchSize],
                                self.w: self.data['train']['sample_weight_SYS'][:batchSize],
                                self.w_toys: self.data['train']['sample_weight_toys'],
                                self.keep_prob: self.pKeep,
                                self.learning_rate_adam: learningRate,
                                self.is_training: True,
                                self.cl: classToSB,
                            }, session=self.session)

            if 'rateGamma' in self.parameters:
                learningRate *= self.parameters['rateGamma']

            if 'dropoutDecay' in self.parameters:
                if self.parameters['dropoutDecay'] > 1.0:
                    self.parameters['pDropout'] = [min(1.0, 1.0-((1.0-x) * math.pow((1.0/self.parameters['dropoutDecay']), self.parameters['nStepsPerEpoch']))) for x in self.parameters['pDropout']]
                else:
                    self.parameters['pDropout'] = [x * math.pow(self.parameters['dropoutDecay'], self.parameters['nStepsPerEpoch']) for x in self.parameters['pDropout']]

                self.pKeep = np.array([1.0 - pDrop for pDrop in self.parameters['pDropout']], dtype=np.float32)
                if self.parameters['dropoutDecay'] != 1.0:
                    print("decay dropout ---> keep:", self.pKeep)

            # ------------------------------------------
            # compute predictions and loss
            # ------------------------------------------
            feed_dict_test = {
                            self.x: self.data['test']['X'],
                            self.y_: self.data['test']['y'].astype(int),
                            #self.c_: self.data['test']['c'],
                            self.w: self.data['test']['sample_weight'],
                            self.w_toys: self.data['test']['sample_weight_toys'],
                            self.keep_prob: np.array([1.0]*len(self.pKeep), dtype=np.float32),
                            self.learning_rate_adam: learningRate,
                            self.is_training: False,
                            self.cl: classToSB,
                        }
            feed_dict_train = {
                            self.x: self.data['train']['X'],
                            self.y_: self.data['train']['y'].astype(int),
                            #self.c_: self.data['train']['c'],
                            self.w: self.data['train']['sample_weight'],
                            self.w_toys: self.data['train']['sample_weight_toys'],
                            self.keep_prob: np.array([1.0]*len(self.pKeep), dtype=np.float32),
                            self.learning_rate_adam: learningRate,
                            self.is_training: False,
                            self.cl: classToSB,
                        }
            feed_dict_train_org = {
                            self.x: self.data_original['train']['X'],
                            self.y_: self.data_original['train']['y'].astype(int),
                            #self.c_: self.data_original['train']['c'],
                            self.w: self.data_original['train']['sample_weight'],
                            self.w_toys: self.data['train']['sample_weight_toys'],
                            self.keep_prob: np.array([1.0]*len(self.pKeep), dtype=np.float32),
                            self.learning_rate_adam: learningRate,
                            self.is_training: False,
                            self.cl: classToSB,
                        }

            probabilities = self.session.run(self.predictions, feed_dict=feed_dict_test)
            probabilities_train = self.session.run(self.predictions, feed_dict=feed_dict_train_org)

            # print loss function
            loss_test = self.total_loss.eval(feed_dict=feed_dict_test, session=self.session)
            loss_train = self.total_loss.eval(feed_dict=feed_dict_train, session=self.session)

            self.computeScoreRescaling({'test':probabilities, 'train':probabilities_train})

            # ------------------------------------------
            # evaluate expected significances
            # ------------------------------------------
            logfile.write("---- test data -----\n")
            resultsTest = self.expected_significance(probabilities, self.truth, self.data['test']['sample_weight'], self.signalClassIDs, classLabels=self.classLabels)
            expected_significance_test = resultsTest['significance']['ssb']
            ssbs = resultsTest['histograms']['ssb']

            logfile.write(resultsTest['log'])
            if (i % 10 == 0):
               print("test data scores:")
               print(resultsTest['log'])

            if i%40 == 0:
                self.saveCheckpoint()

            if expected_significance_test > self.parameters['best_expected_significance']:
                self.parameters['best_expected_significance'] = expected_significance_test
                save_path = self.saver.save(self.session, self.directories['checkpoint_best'])
                print("saved \x1b[31mbest\x1b[0m checkpoint to \x1b[34m", save_path, "\x1b[0m")
                with open('best_checkpoint.txt','w') as outputFile:
                    outputFile.write(save_path)
                    outputFile.write("\n%1.5f"%self.parameters['best_expected_significance'])
                self.saveParameters(self.directories['checkpoint_best'] + '.parameters')

            if expected_significance_test > recordForEarlyStop:
                recordForEarlyStop = expected_significance_test
                recordForEarlyStopTime = (i+1)*self.parameters['nStepsPerEpoch']

            earlyStopEpochsLeft = recordForEarlyStopTime + self.parameters['earlyStopEpochs'] - (i+1)*self.parameters['nStepsPerEpoch']
            if earlyStopEpochsLeft < 0:
                print("early stop triggered!")
                break

            #try:
            #    print ("loss shape S: %1.6f"%self.shape_lossS.eval(feed_dict=feed_dict_train, session=self.session))
            #    print ("loss shape B: %1.6f"%self.shape_loss.eval(feed_dict=feed_dict_train, session=self.session))
            #except:
            #    pass

            # ------------------------------------------
            # produce HISTOGRAMS for tensorboard
            # ------------------------------------------
            hist = tf.HistogramProto()
            hist.min = 0
            hist.max = 15
            hist.num = 15
            hist.sum = sum(ssbs)
            hist.sum_squares = expected_significance_test*expected_significance_test

            bin_edges = list(range(1,len(ssbs)+1))
            # Add bin edges and counts
            for edge in bin_edges:
                hist.bucket_limit.append(edge)
            for c in ssbs:
                hist.bucket.append(c)

            hist2 = tf.HistogramProto()
            hist2.min = 0
            hist2.max = 15
            hist2.num = 15
            # Add bin edges and counts
            for edge in bin_edges:
                hist2.bucket_limit.append(edge)
            for c in resultsTest['histograms']['signal']:
                hist2.bucket.append(c)

            hist3 = tf.HistogramProto()
            hist3.min = 0
            hist3.max = 15
            hist3.num = 15
            # Add bin edges and counts
            for edge in bin_edges:
                hist3.bucket_limit.append(edge)
            for c in resultsTest['histograms']['background']:
                hist3.bucket.append(c)

            # Create and write Summary
            hsummary = tf.Summary(value=[tf.Summary.Value(tag="binSignificance", histo=hist)])
            hsummary2 = tf.Summary(value=[tf.Summary.Value(tag="signal", histo=hist2)])
            hsummary3 = tf.Summary(value=[tf.Summary.Value(tag="bckground", histo=hist3)])

            _, _, result = self.session.run([self.update_tensorboardVar, self.update_tensorboardVar_loss, self.merged], feed_dict={self.pythonVar: expected_significance_test, self.pythonVar_loss: loss_test})
            self.sum_writer_test.add_summary(result, numTrainingEventsTotal)
            self.sum_writer_test.add_summary(hsummary, numTrainingEventsTotal)
            self.sum_writer_test.add_summary(hsummary2, numTrainingEventsTotal)
            self.sum_writer_test.add_summary(hsummary3, numTrainingEventsTotal)
            self.sum_writer_test.flush()

            logfile.write("---- training data -----\n")
            resultsTrain = self.expected_significance(probabilities_train, self.truth_train, self.data_original['train']['sample_weight'], self.signalClassIDs, classLabels=self.classLabels)
            expected_significance_train = resultsTrain['significance']['ssb']

            logfile.write(resultsTrain['log'])
            _, _, result = self.session.run([self.update_tensorboardVar, self.update_tensorboardVar_loss, self.merged], feed_dict={self.pythonVar: expected_significance_train, self.pythonVar_loss: loss_train})
            self.sum_writer_train.add_summary(result, numTrainingEventsTotal)
            self.sum_writer_train.flush()

            nEpochsPassed = (i+1)*self.parameters['nStepsPerEpoch']
            print(("%1.4f"%expected_significance_test).ljust(8), ("%1.4f"%expected_significance_train).ljust(8),("%d"%(nEpochsPassed)).ljust(12), "%d"%numTrainingEventsTotal, " (l:", earlyStopEpochsLeft, ")", "Asimov (test/train): %1.4f %1.4f, overtr = %1.4f"%(resultsTest['significance']['qA'], resultsTrain['significance']['qA'], resultsTrain['significance']['qA']/resultsTest['significance']['qA']))

    def saveCheckpoint(self):
        save_path = self.saver.save(self.session, self.directories['checkpoint'])
        print("saved checkpoint to \x1b[34m", save_path, "\x1b[0m")
        with open('last_checkpoint.txt','w') as outputFile:
            outputFile.write(save_path)
            outputFile.write("\n%1.5f"%self.parameters['best_expected_significance'])
        self.saveParameters()

    def saveParameters(self, parametersFileName=None):
        status = True
        try:
            if not parametersFileName:
                parametersFileName = self.directories['checkpoint'] + '.parameters'
            with open(parametersFileName, 'w') as parametersFile:
                json.dump(self.parameters, parametersFile)
        except Exception as e:
            print("\x1bp31merror: could not save json parameter file!", e, "\x1b[0m")
            status = False
        return status

    def rescaleMVAscore(self, x):
        scoreScaled = (x - self.scoreMin) * self.scoreScaling
        if scoreScaled >= 1:
            scoreScaled = 1 - 1.0e-8
        if scoreScaled < 0:
            scoreScaled = 0
        return scoreScaled

    def computeScoreRescaling(self, probabilities):
        pTest = sorted(probabilities['test'][:, 0])
        self.scoreMin = pTest[int(self.parameters['mvaScoreRescalingPercentileLow'] * len(pTest))]
        self.scoreMax = pTest[int(self.parameters['mvaScoreRescalingPercentileHigh'] * len(pTest))]
        if self.scoreMin < 0:
            self.scoreMin = 0
        self.scoreScaling = 1.0 / (self.scoreMax-self.scoreMin) if (self.scoreMax-self.scoreMin) > 0 else 1.0

    def validate(self, verbose=False):
        try:
            self.createValidationPlots(verbose)
        except:
            pass

    def createValidationPlots(self, verbose=False):
        self.controlPlots = ['SB', 'SBstacked','SBSeff','ROC']
        if verbose:
            self.controlPlots.append('ROC')
        fileTypes = ['png', 'root', 'pdf']

        probabilities = {}
        probabilities['test'] = self.session.run(self.predictions, feed_dict={
                            self.x: self.data['test']['X'],
                            self.y_: self.data['test']['y'].astype(int),
                            self.w: self.data['test']['sample_weight'],
                            self.keep_prob: np.array([1.0]*len(self.pKeep), dtype=np.float32),
                            self.learning_rate_adam: self.parameters['learning_rate_adam_start'],
                            self.is_training: False,
                        })
        probabilities['train'] = self.session.run(self.predictions, feed_dict={
                        self.x: self.data_original['train']['X'],
                        self.y_: self.data_original['train']['y'].astype(int),
                        self.w: self.data_original['train']['sample_weight'],
                        self.keep_prob: np.array([1.0]*len(self.pKeep), dtype=np.float32),
                        self.learning_rate_adam: self.parameters['learning_rate_adam_start'],
                        self.is_training: False,
                    })

        histograms = {
            'test': {0: None, 1: None},
            'train': {0: None, 1: None},
        }

        self.computeScoreRescaling(probabilities)
        self.data['train'] = self.data_original['train']

        if 'SB' in self.controlPlots:
            for ds in ['test', 'train']:
                for cl in [0, 1]:
                    histograms[ds][cl] = ROOT.TH1D(ds + '%d'%cl, ds + '%d'%cl, 30, 0.0, 1.0)
                    histograms[ds][cl].Sumw2()
                    histograms[ds][cl].SetStats(0)

            for i in range(len(self.data['test']['y'])):
                signalBackground = 1 if int(self.data['test']['y'][i]) in self.signalClassIDs else 0
                histograms['test'][signalBackground].Fill(self.rescaleMVAscore(probabilities['test'][i][0]), self.data['test']['sample_weight'][i])

            for i in range(len(self.data_original['train']['y'])):
                signalBackground = 1 if int(self.data_original['train']['y'][i]) in self.signalClassIDs else 0
                histograms['train'][signalBackground].Fill(self.rescaleMVAscore(probabilities['train'][i][0]), self.data_original['train']['sample_weight'][i])

            maxBinContent = 0
            for ds in ['test', 'train']:
                for cl in [0, 1]:
                    histograms[ds][cl].Scale(1.0/histograms[ds][cl].Integral())
                    maxBinContent = max(maxBinContent, histograms[ds][cl].GetBinContent(histograms[ds][cl].GetMaximumBin()))

            for ds in ['test', 'train']:
                for cl in [0, 1]:
                    histograms[ds][cl].GetYaxis().SetRangeUser(0, maxBinContent * 1.5)

            histograms['test'][0].SetLineColor(ROOT.kRed)
            histograms['test'][1].SetLineColor(ROOT.kBlue)
            histograms['test'][0].SetMarkerStyle(21)
            histograms['test'][1].SetMarkerStyle(21)
            histograms['test'][0].SetMarkerSize(0.7)
            histograms['test'][1].SetMarkerSize(0.7)
            histograms['test'][0].SetMarkerColor(ROOT.kRed)
            histograms['test'][1].SetMarkerColor(ROOT.kBlue)

            histograms['train'][0].SetLineColor(ROOT.kRed)
            histograms['train'][1].SetLineColor(ROOT.kBlue)
            histograms['train'][0].SetFillStyle(3004)
            histograms['train'][1].SetFillStyle(1001)
            histograms['train'][0].SetFillColor(ROOT.kRed)
            histograms['train'][1].SetFillColor(ROOT.kBlue-10)

            c1 = ROOT.TCanvas("c1","c1",500,500)
            histograms['train'][1].SetTitle("overtraining check")
            histograms['train'][1].Draw("HIST")
            histograms['train'][0].Draw("HIST same")
            histograms['test'][0].Draw("E same")
            histograms['test'][1].Draw("E same")

            leg = ROOT.TLegend(0.15,0.7,0.8,0.86)
            leg.AddEntry(histograms['train'][0], "BKG (training)")
            leg.AddEntry(histograms['train'][1], "SIG (training)")
            leg.AddEntry(histograms['test'][0], "BKG (test)")
            leg.AddEntry(histograms['test'][1], "SIG (test)")
            leg.SetNColumns(2)
            leg.Draw("same")
            for ext in fileTypes:
                c1.SaveAs(self.directories["output"] +"/overtraining." + ext)
            ROOT.gPad.SetLogy(1)
            ROOT.gPad.Update()

            for ds in ['test', 'train']:
                for cl in [0, 1]:
                    histograms[ds][cl].GetYaxis().SetRangeUser(0.0001, 2.0)

            for ext in fileTypes:
                c1.SaveAs(self.directories["output"] +"/overtraining_log." + ext)

        resultsTest = self.expected_significance(probabilities['test'], self.truth,
                                                                                             self.data['test'][
                                                                                                 'sample_weight'],
                                                                                             self.signalClassIDs,
                                                                                             classLabels=self.classLabels)
        print(resultsTest['log'])

        resultsTrain = self.expected_significance(probabilities['train'], self.truth_train,
                                                                         self.data_original['train']['sample_weight'],
                                                                         self.signalClassIDs,
                                                                         classLabels=self.classLabels)
        print(resultsTrain['log'])

        # ------------------------------------------
        # SBSeff
        # TODO: more efficient computation
        # ------------------------------------------
        if 'SBSeff' in self.controlPlots:
            c1 = ROOT.TCanvas("c1", "c1", 500, 500)
            leg = ROOT.TLegend(0.65,0.65,0.88,0.88)
            tGraphs = []
            first = True
            for dataset in ['train', 'test']:
                print("computing S/B vs S eff for set\x1b[34m", dataset, "\x1b[0m")
                scoreCut = 0
                curve_points_x = []
                curve_points_y = []
                scores = []
                nTotal = [0.0, 0.0]
                for i in range(len(self.data[dataset]['y'])):
                    signalBackground = 1 if int(self.data[dataset]['y'][i]) in self.signalClassIDs else 0
                    score = self.rescaleMVAscore(probabilities[dataset][i][0])
                    weight = self.data[dataset]['sample_weight'][i]
                    nTotal[signalBackground] += weight
                    scores.append([score, signalBackground, weight])
                scores.sort(key=lambda x: x[0], reverse=True)

                nAboveCut = [0.0, 0.0]
                for x in scores:
                    nAboveCut[x[1]] += x[2]
                    signalEfficiency = nAboveCut[1] / nTotal[1]
                    signalOverBackground = nAboveCut[1] / nAboveCut[0] if nAboveCut[0] > 0 else 0
                    curve_points_x.append(signalEfficiency)
                    curve_points_y.append(signalOverBackground)

                arr_x = np.array(curve_points_x, dtype=np.float64)
                arr_y = np.array(curve_points_y, dtype=np.float64)
                tg = ROOT.TGraph(len(arr_x), arr_x, arr_y)
                tg.SetLineColor(ROOT.kRed if dataset == 'train' else ROOT.kBlack)
                tg.Draw("AL" if first else "L")
                tg.GetXaxis().SetTitle("signal efficiency")
                tg.GetYaxis().SetTitle("signal/background")
                tg.GetYaxis().SetRangeUser(0, 1.0)
                leg.AddEntry(tg, 'test set' if dataset == 'test' else 'training set')
                tGraphs.append(tg)
                first = False
            leg.SetFillStyle(0)
            leg.Draw("same")
            ROOT.gPad.Update()
            for ext in fileTypes:
                c1.SaveAs(self.directories["output"] + "/sbseff." + ext)

        # ------------------------------------------
        # ROC curves
        # TODO: more efficient computation
        # ------------------------------------------
        if 'ROC' in self.controlPlots:
            c1 = ROOT.TCanvas("c1", "c1", 500, 500)
            leg = ROOT.TLegend(0.2,0.2,0.5,0.4)
            tGraphs = []
            first = True
            for dataset in ['train', 'test']:
                print("computing ROC curves for set\x1b[34m", dataset, "\x1b[0m")
                scoreCut = 0
                curve_points_x = []
                curve_points_y = []
                scores = []
                nTotal = [0.0, 0.0]
                for i in range(len(self.data[dataset]['y'])):
                    signalBackground = 1 if int(self.data[dataset]['y'][i]) in self.signalClassIDs else 0
                    score = self.rescaleMVAscore(probabilities[dataset][i][0])
                    weight = self.data[dataset]['sample_weight'][i]
                    nTotal[signalBackground] += weight
                    scores.append([score, signalBackground, weight])
                scores.sort(key=lambda x: x[0], reverse=True)

                nAboveCut = [0.0, 0.0]
                for x in scores:
                    nAboveCut[x[1]] += x[2]
                    signalEfficiency = nAboveCut[1] / nTotal[1]
                    backgroundEfficiency = nAboveCut[0] / nTotal[0] if nTotal[0] > 0 else 0
                    backgroundRejection = 1.0 - backgroundEfficiency

                    curve_points_x.append(backgroundRejection)
                    curve_points_y.append(signalEfficiency)

                arr_x = np.array(curve_points_x, dtype=np.float64)
                arr_y = np.array(curve_points_y, dtype=np.float64)
                tg = ROOT.TGraph(len(arr_x), arr_x, arr_y)
                tg.SetLineColor(ROOT.kRed if dataset == 'train' else ROOT.kBlack)
                tg.Draw("AL" if first else "L")
                tg.GetXaxis().SetTitle("background rejection")
                tg.GetYaxis().SetTitle("signal efficiency")
                leg.AddEntry(tg, 'test set' if dataset == 'test' else 'training set')
                tGraphs.append(tg)
                first = False
            leg.SetFillStyle(0)
            leg.Draw("same")
            ROOT.gPad.Update()
            for ext in fileTypes:
                c1.SaveAs(self.directories["output"] + "/roc_curve." + ext)

        # ------------------------------------------
        # S and B stacked
        # ------------------------------------------
        if 'SBstacked' in self.controlPlots:
            sbStacked = ROOT.THStack("hs", "")
            hist = [ROOT.TH1D("bhist","bhist",15,0.0,1.0), ROOT.TH1D("shist","shist",15,0.0,1.0)]
            dataset = 'test'
            for i in range(len(self.data[dataset]['y'])):
                signalBackground = 1 if int(self.data[dataset]['y'][i]) in self.signalClassIDs else 0
                score = self.rescaleMVAscore(probabilities[dataset][i][0])
                weight = self.data[dataset]['sample_weight'][i]
                hist[signalBackground].Fill(score, weight)
            hist[0].SetFillStyle(1001)
            hist[1].SetFillStyle(1001)
            hist[0].SetFillColor(ROOT.kGray)
            hist[1].SetFillColor(ROOT.kRed)
            hist[0].SetLineColor(ROOT.kGray+1)
            hist[1].SetLineColor(ROOT.kRed+1)
            sbStacked.Add(hist[0])
            sbStacked.Add(hist[1])
            c1 = ROOT.TCanvas("c1", "c1", 500, 500)
            sbStacked.Draw("hist")
            leg = ROOT.TLegend(0.5,0.7,0.85,0.88)
            leg.AddEntry(hist[0], "background")
            leg.AddEntry(hist[1], "signal")
            leg.SetFillStyle(0)
            leg.Draw("same")
            ROOT.gPad.Update()
            for ext in fileTypes:
                c1.SaveAs(self.directories["output"] + "/sb_stacked." + ext)

            ROOT.gPad.SetLogy(1)
            ROOT.gPad.Update()
            for ext in fileTypes:
                c1.SaveAs(self.directories["output"] + "/sb_stacked_log." + ext)

    # ------------------------------------------
    # estimate significance with S/sqrt(S+B)
    # ------------------------------------------
    def expected_significance(self, probabilities, truth, weights, signalClassIDs, printout=True, classLabels=None):
        resultObject = {}
        output = ""
        nBins = 15
        histogramS = [0.0] * nBins
        histogramB = [0.0] * nBins
        if classLabels:
            histogramClass = {k: [0.0] * nBins for k in list(classLabels.keys())}

        # compute min/max probabilities
        pMin = 999
        pMax = -999
        for i in range(len(truth)):
            p = sum([probabilities[i][k] for k in signalClassIDs])
            if p > pMax:
                pMax = p
            if p < pMin:
                pMin = p
        if pMax == pMin:
            pMin = 0.0
            pMax = 1.0

        # compute S/sqrt(S+B) for each of the bins
        for i in range(len(truth)):
            p = sum([probabilities[i][k] for k in signalClassIDs])
            p = self.rescaleMVAscore(p)
            binNumber = int(min(max(0, p * nBins), nBins - 1))
            if truth[i] in signalClassIDs:
                histogramS[binNumber] += weights[i]
            else:
                histogramB[binNumber] += weights[i]
            if classLabels:
                histogramClass[truth[i]][binNumber] += weights[i]
        if printout:
            classBinContents = ""
            if classLabels:
                classBinContents += "".join([classLabels[i].ljust(9) for i in list(classLabels.keys())])
            output += "bin     signal    background   S/sqrt(S+B)" + classBinContents + "\n"
        ssbSum = 0.0
        asimovSum = 0.0
        ssbs = []
        for binNumber in range(nBins):
            ssb = histogramS[binNumber] / math.sqrt(histogramS[binNumber] + histogramB[binNumber]) if histogramS[
                                                                                                          binNumber] + \
                                                                                                      histogramB[
                                                                                                          binNumber] > 0 else 0
            if printout:
                classBinContents = ""
                if classLabels:
                    classBinContents += "".join(
                        [("%1.1f" % histogramClass[i][binNumber]).ljust(9) for i in list(classLabels.keys())])
                output += ("%d" % binNumber).ljust(8) + ("%1.1f" % histogramS[binNumber]).ljust(10) + (
                    "%1.1f" % histogramB[binNumber]).ljust(13) + ("%1.3f" % ssb).ljust(13) + classBinContents + "\n"
            ssbSum += ssb * ssb
            ssbs.append(ssb)
            if histogramB[binNumber] > 0:
                asimovSum += 2.0*((histogramS[binNumber] + histogramB[binNumber]) * math.log(1.0 + histogramS[binNumber]/histogramB[binNumber]) - histogramS[binNumber])
        resultObject['histograms'] = {
                'signal': histogramS,
                'background': histogramB,
                'ssb': ssbs,
            }
        resultObject['significance'] = {
                'ssb': math.sqrt(ssbSum),
                'qA': math.sqrt(asimovSum),
            }
        if printout:
            output += "expected: %1.4f\n" % math.sqrt(ssbSum)
        resultObject['log'] = output

        return resultObject
