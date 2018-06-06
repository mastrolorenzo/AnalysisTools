import collections
import glob
import math
import os
from xml.etree import ElementTree

import pandas


Bin = collections.namedtuple('Bin', ['index', 'significance', 'num_events', 'num_signal', 'num_background', 'sumw_signal', 'sumw_background', 'num_data_outside_signal_region', 'num_background_outside_signal_region'])


def parse_rebinning_output(output):
    """Parse the rebinner output and return a dictionary
    containing detailed binning information.

    The dictionary keys are strings of the format "Bin<index>",
    where <index> is the zero-based bin index counting from the
    right, and the values are namedtuples with the following fields:
    * index
      The bin index.
    * significance
      The significance value of the bin.
    * num_events
      The total number of unweighted events in the bin.
    * num_signal
      The number of unweighted signal events in the bin.
    * num_background
      The number of unweighted background events in the bin.
    * sumw_signal
      The number of weighted signal events in the bin.
    * sumw_background
      The number of weighted background events in the bin.
    * num_data_outside_signal_region
      The number of data events outside of the signal region.
    * num_background_outside_signal_region
      The number of background events outside of the signal region.
    """
    bins = {}
    lines = output.splitlines()
    for i, line in enumerate(reversed(lines)):
        # Because we're looping through the output in reverse order,
        # we only parse up to the first line containing "Nodes".
        if 'Nodes' in line:
            break
        else:
            # Remove all whitespace.
            line = ''.join(line.split())
            # Split and discard the leaf node splitting path
            # while keeping the detailed bin information.
            _, info = line.split(':')
            # Store the index and information of the bin in a namedtuple
            fields = [i]
            fields.extend(float(x) for x in info.split(','))
            bins['Bin{0!s}'.format(i)] = Bin(*fields)
    return bins


def parse_rebinning_tree(path):
    """Parse the rebinner XML output file and return a dictionary
    containing rebinning information.

    The dictionary contains the following items:
    * bin_edges : list of floats
      The list of bin edges determined by the rebinner.
      This list does not contain the first and last bin
      edges which define the range of the distribution.
    * total_significance : float
      The total significance of the rebinning scheme.
    """
    bin_significance_squared = []
    bin_edges = []

    def recursive_search(node):
        if len(node) > 0:
            bin_edges.append(float(node.get('splitVal')))
            for child in node:
                recursive_search(child)
        else:
            bin_significance_squared.append(float(node.get('significanceSquared')))

    xmltree = ElementTree.parse(path)
    recursive_search(xmltree.getroot())

    rebinning_info = {
        'bin_edges': sorted(set(bin_edges)),
        'total_significance': math.sqrt(sum(bin_significance_squared)),
    }

    return rebinning_info


def parse_rebinning_trees(src, bounds=None):
    """Parse a multiple rebinner XML output files for rebinning information.

    Parameters
    ----------
    src : path
        The path to the directory containing the rebinner XML output files.
    bounds : tuple of int, optional
        The minimum and maximum values determining the boundaries of the
        rebinned distribution, which are not saved by the rebinner to the
        output XML file. If these are not provided, the "bin_edges" column
        in the returned DataFrame will only contain the internal bin edges.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the following columns:
        * n_bins : int64
          The number of final bins.
        * metric : object (string)
          The name of the significance metric, either "asimov" or "poisson".
        * n_minbkg : int64
          The minimum number of background events required in each bin.
        * smooth_bkg : bool
          Whether background estimates were smoothed by averaging neighboring bins,
          either True or False.
        * unc_tol : float64
          The maximum relative statistical uncertainty for background in each bin.
        * bin_edges : object (list of floats)
          The list of bin edges determined by the rebinner. This also contains the
          boundaries of the rebinned distribution if the `bounds` argument is provided.
        * total_significance : float64
          The total significance of the rebinning scheme.
    """
    columns = {x: [] for x in ['n_bins', 'metric', 'n_minbkg', 'smooth_bkg', 'unc_tol', 'bin_edges', 'total_significance']}
    for path in glob.glob(os.path.join(src, '*.xml')):
        filename, _ = os.path.splitext(os.path.basename(path))
        # Assume that the filenames follow the convention:
        # rebin_{n_bins}_{metric}_{n_minbkg}_{smooth_bkg}_{unc_tol}_{width_tol}
        fields = filename.split('_')
        columns['n_bins'].append(int(fields[1]))
        columns['metric'].append(fields[2])
        columns['n_minbkg'].append(int(fields[3]))
        columns['smooth_bkg'].append(True if fields[4] == 'True' else False)
        columns['unc_tol'].append(float(fields[5].replace('p', '.')))
        result = parse_rebinning_tree(path)
        if bounds:
            columns['bin_edges'].append([bounds[0]] + result['bin_edges'] + [bounds[1]])
        else:
            columns['bin_edges'].append(result['bin_edges'])
        columns['total_significance'].append(result['total_significance'])
    df = pandas.DataFrame(columns)
    return df

