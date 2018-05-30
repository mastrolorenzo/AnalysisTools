## 1. Setup

If this is your first time setting up the workspace, running the `setup` script to check for an install any missing Python packages as well as download and install the Autocategorizer library.

```bash
./setup
```

## 2. Usage

First, create a Rebinner configuration file. The configuration file is simply a Python module which requires the following attributes to be defined:

* `input_pattern`
   The input pattern used to load sample files. The input pattern must contain a "%s" which will be substituted by a sample's input token.
* `destination`
   The output directory containing the rebinner XML files.
* `treename`
   The name of the tree from which to load events.
* `target`
   The name and boundary values of the distribution to rebin.
* `aliases`
  A dictionary mapping the default names for branches expected by the Autocategorizer to the names of their corresponding branches in the tree. The Autocategorizer expects the following branches to be present:

  - `bin`
    The index of the target distribution bin that contains the event. If the event belongs to the signal region and falls within the range of the target distribution, its index should be 0. Otherwise, the index should be -1.
  - `is_signal`
    Whether the tree contains signal (1), background (0), or data (-1) events.flags the Monte-Carlo sample as either signal or background.
  - `weight`
    The full event weight, including everything from b-tagging scale factors and pile-up reweighting to scale factors.

  If equivalent branches are present in the tree but with different names, then this option is needed to remap the branch references. If no aliases are needed, simply pass an empty dictionary "{}".
* `settings`
  A dictionary of settings used to find an optimal rebinning scheme. The following settings are available:

  - `n_bins`
    The number of final bins.
  - `metric`
    The metric used to evaluate the signifance of a rebinning scheme, either "asimov" or "poisson".
  - `n_minbkg`
    The minimum number of **unweighted** background events required in each bin.
  - `smooth_bkg`
    When set to True, the background estimate in a bin is smoothed by averaging with neighboring bins. This can help with downward fluctuations in the background caused by low statistics.
  - `unc_tol`
    The maximum relative background statistical uncertainty allowed per bin. If any bin fails this threshold, the rebinning scheme is rejected.
  - `width_tol`
    The minimum allowable bin width relative to the full range of the distribution. This catches the behaviour where, as the number of final bins gets larger, the rebinner construct bins that are very narrow. If any bin fails this threshold, the rebinning scheme is rejected.

* `samples`
   A list of the samples from which to load events for rebinning. These correspond to the input tokens which are substituted into the input pattern when loading sample files.

For a configuration file example, please see the module `config.py`.

Feel free to use the example configuration files as a template to be copied and modified as needed.

When ready, run the Rebinner using the `run` script which takes a configuration file path as its only argument. The script will then run a grid search over the possible rebinning parameters, generate XML files in the specified output directory, and report the best binning schemes for the two different metrics.

```bash
./run config.py
```

## 3. Visualizing Results

Once the rebinner has finished running, a quick way to plot the rebinned distribution would be nice for debugging purposes since all the events and information to visualize the resulting distribution is readily available. Because each channel has its own requirements on the samples used and branch names, its easier to push the responsibility of visualization to the individual channels. To that end, please see the script `plot_0lep.py` for an example of a quick and dirty plotting procedure. Feel free to copy and adapt this script to your needs.

