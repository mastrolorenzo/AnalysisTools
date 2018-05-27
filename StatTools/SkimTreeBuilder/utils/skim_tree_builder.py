import glob
import os

import concurrent.futures
import pandas
import root_pandas
import uproot

from .makedirs import makedirs


DEFAULT_BRANCHES = [
    # General Branches
    'sampleIndex', 'event', 'Pass_nominal', 'cutFlow', 'twoResolvedJets', 'controlSample', 'weight',
    # Channel Indicators
    'isZnn', 'isWenu', 'isWmunu', 'isZee', 'isZmm',
    # Matching Patterns
    'CMS_*_BDT*', 'bdtInput*',
]


def build_skim_tree(path, additional_branches=[]):
    """Build a skimmed tree from an AnalysisTools output file.

    This function deserializes a subset of branches stored in an
    AnalysisTools output file into a dataframe, filters for signal
    region events, and then splits the dataframe by parity to assign
    odd events for training and even events for testing.

    Parameters
    ----------
    path : path
        The path or XRootD url to the AnalysisTools output file.
    additional_branches : list of strings, optional
        A list of branches to keep in addition to the default
        subset.

    Returns
    -------
    tuple of pandas.DataFrame
        A tuple containing the training and test dataframes, respectively.
    """
    events = uproot.open(path)['Events']

    dataframe = (
        events.pandas.df(DEFAULT_BRANCHES + additional_branches)
        # Generic signal region selection
        .loc[lambda x: x.Pass_nominal]
        .loc[lambda x: x.twoResolvedJets]
        .loc[lambda x: x.controlSample == 0]
    )

    # Select odd events for training
    dataframe_train = dataframe.loc[lambda x: x.event % 2 == 1]

    # Select even events for testing
    dataframe_test = dataframe.loc[lambda x: x.event % 2 == 0]

    return dataframe_train, dataframe_test


class SkimTreeBuilder(object):
    """Build skimmed trees from AnalysisTools output files.

    The skimmed trees are saved in ROOT and HDF5 format.

    Parameters
    ----------
    input_pattern : string
        The path to be interpolated by a sample's input token to glob for its output
        files. This is similar to the input_pattern for PlotWithVarial.
    dst : path, optional
        The path to the output directory containing the preprocessed files.
        The default is a directory named "samples" in the current working directory.
    """
    def __init__(self, input_pattern, dst=None):
        self.input_pattern = input_pattern
        self.dst = dst or 'samples'
        makedirs(self.dst)

    def run(self, input_token, additional_branches=[]):
        """Build a skimmed tree for a sample.

        Parameters
        ----------
        input_token : string
            The sample specific pattern that is substituted into the input_pattern for
            output file globbing. This is similar to the input_token for PlotWithVarial.
        additional_branches : list of strings, optional
            A list of branches to keep in addition to the default
            subset.
        """
        # Preprocess each of the sample's output files in parallel.
        futures = []
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for path in glob.glob(self.input_pattern % input_token):
                futures.append(executor.submit(build_skim_tree, path, additional_branches))
        # As the preprocessing calls complete, collect the intermediate dataframes.
        dataframes_train = []
        dataframes_test = []
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            dataframes_train.append(result[0])
            dataframes_test.append(result[1])
        # Concatenate the intermediate dataframes together.
        dataframe_train = pandas.concat(dataframes_train)
        dataframe_test = pandas.concat(dataframes_test)
        # Save the dataframes as ROOT TTrees.
        root_path = os.path.join(self.dst, input_token + '.root')
        dataframe_train.to_root(root_path, key='train')
        dataframe_test.to_root(root_path, key='test', mode='a')
        # Save the dataframes as HDF5.
        hdf5_path = os.path.join(self.dst, input_token + '.h5')
        with pandas.HDFStore(hdf5_path) as store:
            store['train'] = dataframe_train
            store['test'] = dataframe_test

