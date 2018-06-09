import glob
import os

import concurrent.futures
import numpy
import pandas
import root_pandas
import uproot

from .makedirs import makedirs


DEFAULT_BRANCHES = [
    # General Branches
    'sampleIndex', 'event', 'Pass_nominal', 'cutFlow', 'twoResolvedJets', 'controlSample', 'weight',
    # Channel Indicators
    'isZnn', 'isWenu', 'isWmunu', 'isZee', 'isZmm',
]

# The following branch types are supported.
BRANCH_TYPES = {
    'b': numpy.uint8,
    's': numpy.uint16,
    'i': numpy.uint32,
    'l': numpy.uint64,
    'B': numpy.int8,
    'S': numpy.int16,
    'I': numpy.int32,
    'L': numpy.int64,
    'F': numpy.float32,
    'D': numpy.float64,
    'O': numpy.bool,
}


class TreeNotFoundError(KeyError):
    pass


def build_skim_tree(path, keep_branches=[], new_branches=[], scale_factors=None):
    """Build a skimmed tree from an AnalysisTools output file.

    This function deserializes a subset of branches stored in an
    AnalysisTools output file into a dataframe, filters for signal
    region events, and then splits the dataframe by parity to assign
    odd events for training and even events for testing.

    Parameters
    ----------
    path : path
        The path or XRootD url to the AnalysisTools output file.
    keep_branches : list of strings, optional
        A list of existing branches to be kept in addition to the default subset.
    new_branches : list of callables, optional
        A list of functions defining new branches for the skimmed tree.
        The functions must take a dataframe row (an event) as its only argument and
        return the value for the new branch for that event. The function's name is
        assigned as the name of the new branch. Note that if the function shares its
        name with an existing branch, the existing branch will be overwritten.
    scale_factors : callable, optional
        A function which takes a dataframe row (an event) as its
        only argument and applies scale factors to the weight branch.

    Returns
    -------
    tuple of pandas.DataFrame
        A tuple containing the training and test dataframes, respectively.
    """
    f = uproot.open(path)
    try:
        events = f['Events']
    except KeyError:
        raise TreeNotFoundError('Unable to find a TTree named "Events" in {0}'.format(path))
    except Exception:
        raise

    dataframe = (
        events.pandas.df(DEFAULT_BRANCHES + keep_branches)
        # Apply a general signal region selection.
        .loc[lambda x: x.Pass_nominal]
        .loc[lambda x: x.twoResolvedJets]
        .loc[lambda x: x.controlSample == 0]
    )

    # If no events passed the signal region selection, return a sentinel value.
    if dataframe.empty:
        return None

    # Create new branches.
    for new_branch in new_branches:
        bname, btype = new_branch.__name__.rsplit('_', 1)
        dataframe[bname] = dataframe.apply(new_branch, axis=1)
        dataframe[bname] = dataframe[bname].astype(BRANCH_TYPES[btype])

    # If provided, apply scale factors.
    if scale_factors:
        dataframe['weight'] = dataframe.apply(scale_factors, axis=1)
    dataframe['weight'] = dataframe['weight'].astype(BRANCH_TYPES['F'])

    # Select odd events for training and even events for testing.
    dataframe_train = dataframe.loc[lambda x: x.event % 2 == 1]
    dataframe_test = dataframe.loc[lambda x: x.event % 2 == 0]
    
    # Scale Monte-Carlo events by a factor of two to account for the splitting.
    if any(dataframe.sampleIndex != 0):
        dataframe_train = dataframe_train.assign(weight=lambda x: x.weight * 2)
        dataframe_test = dataframe_test.assign(weight=lambda x: x.weight * 2)

    return dataframe_train, dataframe_test


class SkimTreeBuilder(object):
    """Build skimmed trees from AnalysisTools output files for MVA training.

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

    def run(self, name, input_tokens, keep_branches=[], new_branches=[], scale_factors=None):
        """Build a skimmed tree for a sample.

        Parameters
        ----------
        name : string
            The name for the files containing the sample's skimmed trees.
        input_tokens : list of strings
            The specific patterns that are substituted into the input_pattern when globbing
            for sample output files. This is similar to the input_token for PlotWithVarial.
        keep_branches : list of strings, optional
            A list of existing branches to be kept in addition to the default subset.
        new_branches : list of callables, optional
            A list of functions defining new branches for the skimmed tree.
            The functions must take a dataframe row (an event) as its only argument and
            return the value for the new branch for that event. The function's name is
            assigned as the name of the new branch. Note that if the function shares its
            name with an existing branch, the existing branch will be overwritten.
        scale_factors : callable, optional
            A callable which takes a dataframe row (an event) as its
            only argument and applies scale factors to the weight branch.
        """
        # Preprocess each of the sample's output files in parallel.
        futures = []
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for input_token in input_tokens:
                for path in glob.glob(self.input_pattern % input_token):
                    futures.append(executor.submit(build_skim_tree, path, keep_branches, new_branches, scale_factors))
        # As the preprocessing calls complete, collect the intermediate dataframes.
        dataframes_train = []
        dataframes_test = []
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            # Catch sentinel values indicating no events passed the generic signal region selection.
            if result is not None:
                dataframes_train.append(result[0])
                dataframes_test.append(result[1])
        # Concatenate the intermediate dataframes together.
        dataframe_train = pandas.concat(dataframes_train)
        dataframe_test = pandas.concat(dataframes_test)
        # Save the dataframes as ROOT TTrees.
        root_path = os.path.join(self.dst, name + '.root')
        dataframe_train.to_root(root_path, key='train')
        dataframe_test.to_root(root_path, key='test', mode='a')
        # Save the dataframes as HDF5.
        hdf5_path = os.path.join(self.dst, name + '.h5')
        with pandas.HDFStore(hdf5_path) as store:
            store['train'] = dataframe_train
            store['test'] = dataframe_test

