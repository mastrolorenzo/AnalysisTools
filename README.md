# AnalysisTools
The basic idea of AnalysisTools is to compile C++ code that loops over events from TTrees (per systematic variation) 
and execute that code with python.  The central class of AnalysisTools is AnalysisManager from which all
analysis code (in "plugins") is inherited.  The structure of the AnalysisManager fixes how the events are 
processed.  The main methods in inherited classes (such as VHbbAnalysis) decides which events to save
(in Preselection() and Analyse()) and what additional information to be computed/saved (in FinishEvent()).
The output of the program is ROOT file with a TTree with all branches used and all new branches added 
for events where the selection is passed (among any systematic variation).

## Quickstart instructions 
### Download and compile
```
cmsrel CMSSW_10_2_0_pre4
cd CMSSW_10_2_0_pre4/src
cmsenv
git clone https://github.com/capalmer85/AnalysisTools
cd AnalysisTools
source env.sh
make -j 4
```
### Running the VHbb code
```
cd VHbbAnalysis
```

This is the most basic version that will run over all samples with standard configuration which is set by vhbb_config.txt and files linked in it.   
```
python RunAnalysis.py -c vhbb_config_2017.txt
```

For a unit test something like this would be better.  This runs over 1\% of events where the sample name contains 'Wplus'.
```
python RunAnalysis.py -c vhbb_config_2017.txt --sample WplusH125_powheg -f 0.01
```

## Key configuration files
Many things can be configured with text files.  Which files to run, which BDTs to evaluate, what selection to 
apply.  Within the main configuration file (above vhbb_config_2017.txt) there are numerous configurations files 
linked.  Some types can be loaded multiple times (e.g. BDTs), others (e.g. samples, early branches, existing 
branches and new branches) should only be loaded once.

### Samples
```
samples=cfg/samples_2017.txt
```

The samples text files linked should have two main features:
1. the common part of the path to the root files
```
prefix=/eos/uscms/store/group/lpchbb/VHbbNanoPostProc/2017/V5/
```
2. Information per sample including:  the "name" of the sample (output directory name), the rest 
of the path from the common path ("dir"), a unique number ("type") to identify the sample in the 
compiled code and the cross section ("xsec") in pb. The type is 0 for data, negative for signal and
positive for background.
```
name=WplusH125_powheg           dir=WplusH_HToBB_WToLNu_M125_13TeV_powheg_pythia8/   type=-12500   xsec=0.17202
```

### Branches
Not all branches in the input trees are needed for the analysis.  In fact only some branches are needed for the 
inital preselection.  That is the idea of early branches.  Existing branches is everything needed beyond preselection
and new branches is for addtional branches added on top of the others.
```
earlybranches=cfg/earlybranches.txt
existingbranches=cfg/existingbranches.txt
newbranches=cfg/newbranches.txt
```

Each one of these files contains of lines like these:
```
type=8  name=Jet_pt           max=250     lengthBranch=nJet
```

The type number is defines the object in this way:
```
unsigned int=0, int=1, float=2, double=3, bool=4, char=5
array +6
```

This is one of the nicer features of the code.  The setup of branches is automatic and 
done fully via configuration.


### Settings
The other main component of the configuration is a settings files where cuts can be 
adjusted without recompiling.
```
settings=cfg/settings.txt
```

These are set up with the same syntax to the branches but only floats are allowed.  
```
type=2  name=elMetDPhiCut       val=2.0
```

All of these settings are saved in a separate tree in the output file(s).  This way you 
can look to see what cuts were applied in the files themselves.


## Making stacked plots (with Varial)
In general the output from AnalysisManager inherited classes will be a series of trees
containing per event weight, event class selection, event variables and evaluated 
discriminator outputs.  A separate set of code is used to open these trees and make
stacked plots.  

###Using PlottingTools/PlotWithVarial

In PlottingTools/PlotWithVarial you can find a plotting tool that runs in parallel 
and is configured in python.  It creates a webpage on which the plots are shown.  
The plots in different categories are linked. The tool is configured with a config file
 where a couple of variables are defined. This is an example:

```
name = 'VHbbPlots'
input_pattern = '/some/path/V25_VHbb_runonskim_20180212/haddjobs/sum_%s.root'
weight = 'weight'
enable_reuse_step = True  # try to find output on disk and don't run a step if present

from main_samples import the_samples_dict, sample_colors
from main_selections import the_category_dict
```

The last two lines import the samples and category definitions. Have a look at these 
files to see how the_category_dict and the_samples_dict are defined. This exact structure 
is expected by the tool.

IMPORTANT: input_pattern must contain "%s"! This is where the input token from the sample definition is inserted. Wildcards are allowed.

Alter input_pattern in the config file an run this script.  For the first time running 
it might be good to start with demo_config.py, which is identical to main_config.py, 
but only runs two samples in the signal region:

```
./run_plot_from_tree.py <demo_config.py>
```

