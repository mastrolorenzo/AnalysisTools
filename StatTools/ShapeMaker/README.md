# Shape Maker for Datacards, VHbb Analysis

# VHbb ShapeMaker Overview
The VHbb ShapeMaker takes four input dictionaries and outputs histograms for each process, channel, and systematic specified in the dictionaries.  The workflow goes as follows: 1.) Define the binning, samples, channels, and systematics in the Dictionaries/ folder 2.) Configure condor_config to use your ntuples, which should be hadded by sample name in a folder 3.) Run ```python runJobs.py``` to submit the condor jobs, which will output many root files in the ShapeMaker directory 4.) Run ```python dohadd.py``` to combine the root files, with the final files being produced in haddlevel3/.

## I. Dictionaries
The main change from the VHcc ShapeMaker to the VHbb ShapeMaker is that before, the old-style sample splitter runCards were used to generate the dictionaries on-the-fly that the ShapeMaker used.  This method has been replaced with explicit editable dictionaries.  A side-effect of this is that the getSRbinning.py script is broken.  However, binning can still be manually specified in VHbb_binsDict.py. 

### VHbb_channelDict.py
This defines the channels.  “channelName” is the lepton channel (“Znn”, “Wln”, or “Zll”), “cutstring” defines the selection, “nBins”, “xhigh”, “xlow”, and “varName” specify the plot settings, and “doEven”, “doRebin”, “doVV”, and “drawFromNom” are options.  “doRebin” tells the plotter to use the bins specified in VHbb_binsDict for the process instead of evenly-spaced bins. “drawFromNom” turns systematics on or off for the sample, and “doEven” only uses even events. (n.b. “doVV” does nothing currently)

### VHbb_binsDict.py
This specifies the bin edges for plots, in case “doRebin” is true for the given process.  The binning should match the range stated in VHbb_channelDict.py for the given channel. 

### VHbb_fileDict.py
This defines which sample files contribute to which processes.  Note that this contains a dictionary for Znn, Wln, and Zll.  Each sample file (e.g. “DYBJets-Pt100To200.root”) gets a key (“DYBJets-Pt100To200”) and, for each process the file contributes to, a legend, cutstring, adjustable weight, and a sample name (the weight and samplename have no effect currently).  The cutstring can be a number to specify a sampleIndex, or a generic ROOT cutstring.

### VHbb_sampleDict.
This consists of three dictionaries for Znn, Wln, and Zll that specify the relevant systematics.  Each systematic is a key paired with a shape, a scale, the branch(es) containing the weight variation, and a list of processes the systematic should be applied to.  These dictionaries follow the same structure as the older systematics.txt files (see References/systematics_Wlnshapes2017STXS.txt).

*(The following is taken and lightly modified from the VHcc readme.)*

## II. Compatibility and setup

This program makes use of ROOT RDataFrame, which is available in ROOT 6.14 and above. With ROOT 6.12, this program falls back to the experimental TDataFrame which is significantly slower.

### 1. SL7/EL7/CentOS7

Assuming the Analysis Tools repo was cloned under CMSSW 10.2, running `cmsenv` will source ROOT 6.12. Therefore, it is necessary to run the following before running the ShapeMaker to ensure best performance:

```
export ORIG_DIR=$PWD
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_10_6_0/src/
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd $ORIG_DIR
```

### 2. SL6 with RDataFrame compatibility

It appears there are no readily available builds of ROOT 6.14+ or CMSSW 10.6 for SL6. Therefore, the following "hack" using Anaconda may be used:

```
cd ~
wget https://repo.anaconda.com/archive/Anaconda2-2019.03-Linux-x86_64.sh
chmod +x Anaconda2-2019.03-Linux-x86_64.sh
./Anaconda2-2019.03-Linux-x86_64.sh
```

Select a location on a filesystem that has no limit on the file count, (eg. /nfs/dust/cms/user/username/anaconda2) as installation location.
Log out and log in after installation is complete.

```
conda install -c conda-forge root
pip install rootpy
```

This will install ROOT 6.16 (or above) on a virtual conda environment (takes quite some time). From a different machine (eg. condor cores), this environment can be sourced using
`export PATH="$CONDAPATH/bin:$PATH"`
where $CONDAPATH is the Anaconda installation directory.

### 3. SL6 without RDataFrame (not recommended)

Simply run `cmsenv` inside CMSSW 10.2/10.3. Note that this makes ShapeMaker slower than usual.

## III. Running locally

### 1. Getting SR binning (Currently broken for VHbb)

The binning of the signal regions needs to be obtained first using `getSRBinning.py`. Example:

```
python getSRBinning.py -d /nfs/dust/cms/user/lmastrol/VHbbAnalysisNtuples/VHccNtuple_March16_forApproval_fixZnnHF/haddjobs/ --multithread
```

This program takes ~1 hour and needs to be run only once for a given value of V_pT cut, in principle. A different V_pT cut can be specified using `--vptcut 300`, for example.

### 2. The main ShapeMaker code

The main script is `VHcc_shapeMaker.py`. The usage and arguments to the main script are as follows:

```
usage: Produce shape histograms for datacards. [-h] [-i INPUTFILE]
                                               [-d INPUTDIR]
                                               [-samp SAMPLENAME] [-n NOUTS]
                                               [-io IOUT]
                                               [--vptcut VPTCUT]
                                               [--skipsysts]
                                               [--dochannel DOCHANNEL]
                                               [--multithread]
  -i INPUTFILE, --inputfile INPUTFILE
                        The input root file (ntuple).
  -d INPUTDIR, --inputdir INPUTDIR
                        The directory where all the ntuples are stored.
  -samp SAMPLENAME, --samplename SAMPLENAME
                        The name of the sample, eg, 'DYToLL_HT100to200'.
  -n NOUTS, --nOuts NOUTS
                        The max. number of output files to produce. Use 0 to
                        process full input file.
  -io IOUT, --iOut IOUT
                        The index of the first output file that will be
                        produced if iOut!=0. So this instance of the program
                        will produce ith to (i+n-1)th outputs.
  --vptcut VPTCUT       The cut on V_pT, required for combination with boosted
                        analysis.
  --skipsysts           Skip systematics.
  --dochannel DOCHANNEL
                        Do only specified channel.
  --multithread         Enable Multithreading for RDataFrame.
```

For example, one can locally run over all TT samples for the first 5 channels using
```
python VHcc_shapeMaker.py -d /nfs/dust/cms/user/lmastrol/VHbbAnalysisNtuples/VHccNtuple_March16_forApproval_fixZnnHF/haddjobs/ -i TT_powheg -io 0 -n 5 --multithread
```
(this command looks for all files with fileglob \*TT_powheg\* inside the input directory and processes channels 0 to 4. Note that either one of -i or -samp needs to be specified.)

## IV. Running batch jobs on HTCondor

***Warning: The condor jobs create numerous output and log files, hence it is recommended that the jobs are submitted from a disk space that is not limited by file count. Eg: On DESY, use /nfs/dust/ instead of the /afs/ partition.***

**Step 1:** Run `python getSRBinning.py [--multithread]` (change input directory as required) if not run already. This should take ~1 hour. 

**Step 2:** Edit `condor_config.py` to set the various parameters involved in condor job submission. Submit jobs to condor using `python runJobs.py`. It will submit ~2500 jobs. All jobs should be over in 3-4 hours in DESY condor.

 - **Step 2.1:** Some jobs that end up on low memory machines may fail due to... :drumrolls:  low memory. Release them with `condor_release <username>` if you don't want to wait for the automatic periodic release.
    
 - **Step 2.2:** Check for jobs that have failed due to exceeding time limit using: ```grep -Ril aborted condor/*.log```
    
 - **Step 2.3:** Check for jobs that have failed due to corrupt input file or read error: ```find condor/*.err -not -empty```
    
 - **Step 2.4:** Resubmit jobs from 2.2: Run `python resubmitAborted.py`. #TODO: Resubmit jobs from 2.3
    
**Step 3:** After ***all*** jobs are over, run `python dohadd.py` to "hadd" the outputs. The final hadd'ed histograms can be found in the `./haddlevel*` directories.
