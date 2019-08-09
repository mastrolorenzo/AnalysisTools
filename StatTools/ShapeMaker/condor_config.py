#indir = "/nfs/dust/cms/user/lmastrol/VHbbAnalysisNtuples/VHccNtuple_March16_forApproval_fixZnnHF/haddjobs/"
#indir = "/eos/user/n/nihaubri/VHbbAnalysisNtuples/vhbb_2017_v11_test/haddjobs"
#indir = "/afs/cern.ch/work/n/nihaubri/private/CMSSW_10_2_10/src/AnalysisTools/StatTools/ShapeMaker/samples"
#indir = "/eos/user/n/nihaubri/VHbbAnalysisNtuples/vhbb_2017_V11_2/haddjobs"
indir = "/eos/user/n/nihaubri/VHbbAnalysisNtuples/vhbb_2017_V11_3/haddjobs2"
vptcut = ""                             # "" for none, "300" for 300 GeV
nChannelsPerJob = 4 #was 14
condorDir = "condor"
sizeperbunch = 2*1024*1024*1024       # O(3.5 GiB) data per job
#sizeperbunch = 3.5*1024*1024*1024       # O(3.5 GiB) data per job
#isSL6 = True                            # Anaconda with ROOT and rootpy must be installed. Edit condapath.
isSL6 = False

skipList = ["sum_WW_0.root","sum_WW_1.root","sum_WW_2.root","sum_WW_3.root","sum_WZ_0.root","sum_WZ_1.root","sum_ZZ_0.root"]
skipString = '_fil_'                    # Files with this string in the filename will be skipped

subText = '''universe = vanilla
Executable     =  condor_runscript.sh
Should_Transfer_Files     = YES
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
Notification     = never
transfer_input_files = VHbb_shapeMaker.py, Dictionaries/VHbb_binsDict.py, Dictionaries/VHbb_fileDict.py, Dictionaries/VHbb_channelDict.py, Dictionaries/VHbb_sampleDict.py, Dictionaries/__init__.py
+MaxRuntime = 4*60*60
WhenToTransferOutput=On_Exit
requirements = OpSysAndVer == "CentOS7"
#request_cpus = 2
'''

#+JobFlavour = "longlunch"
#subText = '''universe = vanilla
#Executable     =  condor_runscript.sh
#Should_Transfer_Files     = YES
#on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
#Notification     = never
#transfer_input_files = VHbb_shapeMaker.py, nano_samples_2017stxs.py, VHbb_binsDict.py, makeChannelDict.py, makeFileBasedSampleList.py, makeSampleBasedSystList.py, runCards_CR_VH_bdt_Wln.sh, runCards_CR_VH_bdt_Znn.sh, runCards_SR_VH_bdt_Zll.sh, runCards_CR_VH_bdt_Zll.sh, runCards_SR_VH_bdt_Wln.sh, runCards_SR_VH_bdt_Znn.sh, systematics_Zllshapes2017STXS.txt, systematics_Znnshapes2017STXS.txt, systematics_Wlnshapes2017STXS.txt
#+JobFlavour = "workday"
#WhenToTransferOutput=On_Exit
#requirements = OpSysAndVer == "CentOS7"
##request_cpus = 2
#'''


logFileText = '''\nOutput     = CONDORDIR/log_$(Cluster)_$(Process).stdout
Error      = CONDORDIR/log_$(Cluster)_$(Process).stderr
Log        = CONDORDIR/log_$(Cluster)_$(Process).log
'''

if not isSL6:
#    # SL7 settings    
    runscrText='''#!/bin/bash
export ORIG_DIR=$PWD
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_10_6_0/src/
source     /cvmfs/cms.cern.ch/cmsset_default.sh
eval     `scramv1 runtime -sh`
cd $ORIG_DIR
mkdir Dictionaries
cp *Dict.py Dictionaries/
cp __init__.py Dictionaries/
'''

else:
    # SL6 compatibility for ROOT RDFs
    condapath = "/nfs/dust/cms/user/spmondal/anaconda2"
    runscrText="export PATH=\"CONDAPATH/bin:$PATH\"".replace("CONDAPATH",condapath)


runscrText += "\npython VHbb_shapeMaker.py -i $1 -n $2 -io $3 -d $4 VPTCUT"

