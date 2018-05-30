import math
import os

import wurlitzer

import autocategorizer
from .makedirs import makedirs
from .parsing import parse_rebinning_output, parse_rebinning_tree


class BinUncertaintyError(Exception):
    pass


class BinWidthError(Exception):
    pass


def rebin(x, features, n_bins, metric='poisson', n_minbkg=100, smooth_bkg=False, unc_tol=0.35, min_width=None, outdir=None):
    """Rebin a distribution.

    The initial distribution must start as a single bin.

    This generates an XML file in the output directory named with the parameters used:
    rebin_{n_bins}_{metric}_{n_minbkg}_{smooth_bkg}_{unc_tol}_{width_tol}.xml

    Parameters
    ----------
    x : autocategorizer.Events
        The training events.
    features : autocategorizer.Features
        The training features.
    n_bins : int
        The target number of bins.
    metric : string, optional
        The proxy significance metric, either:
        * asimov
        * poisson (default)
    n_minbkg : int, optional
        The minimum number of background events required in each
        target bin. The default is 100.
    smooth_bkg : bool, optional
        Smooth the background estimate in a bin by averaging with neighboring
        bins. This can help with downward fluctuations in the background caused
        by low statistics. The default is to not apply smoothing.
    unc_tol : float, optional
        The maximum relative background statistical uncertainty allowed per
        bin for the rebinning scheme to be acceptable. The default is 0.35.
    min_width : iterable of numeric, optional
        Set the minimum allowed bin width relative to the full range of the distribution.
        The argument should pass in order the minimum value of the distribution, the
        maximum value of the distribution, and a width tolerance between 0 and 1. For
        example, passing [-1, 1, 0.05] causes rebinning schemes with bin widths less
        than 5% of the BDT output range of [-1, 1] to be rejected. The default is no minimum.
    outdir : path, optional
        The output directory for the XML output files. The default is a directory named
        "rebinning_schemes" in the current working directory.

    Returns
    -------
    result : dict
        The rebinning result dictionary with the following items:
        * Namedtuples containing individual bin information.
        * The list of bin edges. If the width_tol argument was
          passed, the list also includes the boundary values.
        * The total significance for the binning scheme.
    """
    if metric == 'poisson':
        scorer = autocategorizer.PoissonSignificance(0, 0, n_minbkg, False, False, smooth_bkg)
    elif metric == 'asimov':
        scorer = autocategorizer.AsimovSignificance(0, 0, n_minbkg, False, False, smooth_bkg)
    else:
        raise ValueError("Expected 'poisson' or 'asimov', found {0!s}".format(metric))

    filename = 'rebin_{0!s}_{1}_{2!s}_{3!r}_{4!s}_{5!s}'.format(n_bins, metric, n_minbkg, smooth_bkg, unc_tol, min_width[-1] if min_width else 0)
    filename = filename.replace('.', 'p') + '.xml'
    if outdir is None:
        outdir = 'rebinning_schemes'
    makedirs(outdir)
    outfile = os.path.join(outdir, filename)

    # Train the autocategorizer.
    tree = autocategorizer.Tree(x, 1)
    tree.setFeatureNames(features)
    with wurlitzer.pipes() as (stdout, stderr):
        tree.buildTree(n_bins, scorer)
    bins = parse_rebinning_output(stdout.read())

    # Enforce the bin-wise relative background statistical uncertainty tolerance.
    bins_below_unc_tol = []
    for _, _bin in bins.iteritems():
        sumw_background = _bin.sumw_background
        if math.sqrt(sumw_background) / sumw_background > unc_tol:
            bins_below_unc_tol.append(_bin.index)

    if bins_below_unc_tol:
        failed = ', '.join(str(i) for i in sorted(bins_below_unc_tol))
        raise BinUncertaintyError('The relative background statistical uncertainty for bins {0} are above the required tolerance of {1!s}. Rebinner output:\n{2!r}'.format(failed, unc_tol, bins))

    # Persist the rebinning scheme. This is required to access the bin edges.
    tree.saveToXML(outfile)
    results = parse_rebinning_tree(outfile)
    results.update(bins)

    # Enforce the relative bin width tolerance.
    if min_width:
        min_value, max_value, width_tol = min_width
        total_range = max_value - min_value
        bin_edges = [min_value] + results['bin_edges'] + [max_value]

        bins_below_width_tol = []
        for i, j in enumerate(reversed(xrange(len(bin_edges) - 1))):
            if (bin_edges[j + 1] - bin_edges[j]) / total_range < width_tol:
                bins_below_width_tol.append(i)

        if bins_below_width_tol:
            os.remove(outfile)
            failed = ', '.join(str(i) for i in sorted(bins_below_width_tol))
            raise BinWidthError('The relative widths for bins {0} are below the required tolerance of {1!s}. Rebinner output:\n{2!r}'.format(failed, width_tol, results))

    return results

