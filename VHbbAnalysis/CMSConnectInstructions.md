# CMSConnect VHbbAnalysis Jobs

## Setup Working Area

With your CMSConnect account ready, log in the CMSConnect submission node.

```bash
ssh <username>@login.uscms.org

# You'll need to copy over your GRID certificates from lxplus
# for the first time and when you update them.
copy_certificates
```

Now request a voms-proxy and set up environment variables.

```bash
voms-proxy-init -voms cms -valid 192:00
# You can find where your proxy is using voms-proxy-info
export X509_USER_PROXY=path/to/your/proxy
# Magic variable that the LPC sys admins have determined
# allows access to LPC EOS via XRootD. You can put this
# in your ~/.bashrc
export XrdSecGSISRVNAMES=cmseos.fnal.gov
```

Now create the working area.

```bash
cmsrel CMSSW_10_2_0_pre4
cd CMSSW_10_2_0_pre4/src
cmsenv
git clone git@github.com:capalmer85/AnalysisTools.git
cd AnalysisTools
make
cd VHbbAnalysis
```

## Submitting Jobs

The state of this branch already reflects the state at the time of job submission.

The following commands are what I used to create the job submission files. Although this step takes a while because the workflow wasn't designed for this exact procedure, it was necessary because different samples required different file splittings.

I executed the following commands which you could make into a bash script. I submitted only for those samples which were not commented out in samples\_2017.txt.

Thanks to @capalmer85 for providing the file splitting estimates, which are located here:
/uscms/home/capalmer/nobackup/2018/FirstHalf/CMSSW\_10\_2\_0\_pre4/src/AnalysisToolsProduction/VHbbAnalysis/list2017.txt

```bash
# Data
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.1 -s Run2017_Ele_ReMiniAOD
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.055 -s Run2017_Mu_ReMiniAOD
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.2 -s Run2017_MET_MiniAOD
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.2 -s Run2017_DoubleEle_ReMiniAOD
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.055 -s Run2017_DoubleMu_ReMiniAOD

# I processed the signal samples at LPC, so no need to submit them again.

# W+Jets
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s WJets_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.34 -s WJets-HT100To200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s WJets-HT200To400
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s WJets-HT400To600
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s WJets-HT600To800
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.25 -s WJets-HT800To1200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s WJets-HT1200To2500

# DY+Jets
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.25 -s DYToLL_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.25 -s DYToLL_HT100to200_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.25 -s DYToLL_HT200to400_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s DYToLL_HT400to600_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s DYToLL_HT600to800_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s DYToLL_HT800to1200_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.34 -s DYToLL_HT1200to2500_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s DYToLL_HT2500toInf_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.5 -s DYToLL_M4to50_HT70to100_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s DYToLL_M4to50_HT100to200_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s DYToLL_M4to50_HT200to400_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.5 -s DYToLL_M4to50_HT400to600_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s DYToLL_M4to50_HT600toInf_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s DY1JetsToLL_highM
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.17 -s DY2JetsToLL_highM
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s DY3JetsToLL_highM

# Z+Jets
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s ZJetsToNuNu_HT100To200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s ZJetsToNuNu_HT200To400
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s ZJetsToNuNu_HT400To600
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.34 -s ZJetsToNuNu_HT600To800
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.25 -s ZJetsToNuNu_HT800To1200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s ZJetsToNuNu_HT2500ToInf

# QCD
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s QCD_HT200to300
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s QCD_HT300to500
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s QCD_HT500to700
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s QCD_HT700to1000
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s QCD_HT1000to1500
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s QCD_HT1500to2000
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s QCD_HT2000toInf

# Single Top
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s ST_tW_antitop_5f_inc
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s ST_tW_top_5f_inc_PSw
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s ST_t-c_top_4f_inc
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s ST_s-c_4f_lep_PSw

# TTbar
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.2 -s TT_DiLep
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s TT_SingleLep
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s TT_AllHadronic

# Diboson
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s WW
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.5 -s WZ
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.5 -s ZZ
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s WW_1L1Nu2Q
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 1 -s WW_LNu2Q_nlo
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June2_2017V5 -c vhbb_config_2017.txt -f 0.25 -s ZZ_4L
```

If you are confident you don't need to inspect the job submission files, then remove the `--submitJobs 0` option to submit them immediately.

To submit them after the fact, you can execute the following in a bash script, because calling condor\_submit on a glob pattern doesn't work.

```bash
condor_submit CMSConnect_June2_2017V5/DY1JetsToLL_highM/DY1JetsToLL_highM.submit
condor_submit CMSConnect_June2_2017V5/DY2JetsToLL_highM/DY2JetsToLL_highM.submit
condor_submit CMSConnect_June2_2017V5/DY3JetsToLL_highM/DY3JetsToLL_highM.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_HT100to200_madgraph/DYToLL_HT100to200_madgraph.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_HT1200to2500_madgraph/DYToLL_HT1200to2500_madgraph.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_HT200to400_madgraph/DYToLL_HT200to400_madgraph.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_HT2500toInf_madgraph/DYToLL_HT2500toInf_madgraph.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_HT400to600_madgraph/DYToLL_HT400to600_madgraph.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_HT600to800_madgraph/DYToLL_HT600to800_madgraph.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_HT800to1200_madgraph/DYToLL_HT800to1200_madgraph.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_M4to50_HT100to200_madgraph/DYToLL_M4to50_HT100to200_madgraph.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_M4to50_HT200to400_madgraph/DYToLL_M4to50_HT200to400_madgraph.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_M4to50_HT400to600_madgraph/DYToLL_M4to50_HT400to600_madgraph.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_M4to50_HT600toInf_madgraph/DYToLL_M4to50_HT600toInf_madgraph.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_M4to50_HT70to100_madgraph/DYToLL_M4to50_HT70to100_madgraph.submit
condor_submit CMSConnect_June2_2017V5/DYToLL_madgraph/DYToLL_madgraph.submit
condor_submit CMSConnect_June2_2017V5/QCD_HT1000to1500/QCD_HT1000to1500.submit
condor_submit CMSConnect_June2_2017V5/QCD_HT1500to2000/QCD_HT1500to2000.submit
condor_submit CMSConnect_June2_2017V5/QCD_HT2000toInf/QCD_HT2000toInf.submit
condor_submit CMSConnect_June2_2017V5/QCD_HT200to300/QCD_HT200to300.submit
condor_submit CMSConnect_June2_2017V5/QCD_HT300to500/QCD_HT300to500.submit
condor_submit CMSConnect_June2_2017V5/QCD_HT500to700/QCD_HT500to700.submit
condor_submit CMSConnect_June2_2017V5/QCD_HT700to1000/QCD_HT700to1000.submit
condor_submit CMSConnect_June2_2017V5/Run2017_DoubleEle_ReMiniAOD/Run2017_DoubleEle_ReMiniAOD.submit
condor_submit CMSConnect_June2_2017V5/Run2017_DoubleMu_ReMiniAOD/Run2017_DoubleMu_ReMiniAOD.submit
condor_submit CMSConnect_June2_2017V5/Run2017_Ele_ReMiniAOD/Run2017_Ele_ReMiniAOD.submit
condor_submit CMSConnect_June2_2017V5/Run2017_MET_MiniAOD/Run2017_MET_MiniAOD.submit
condor_submit CMSConnect_June2_2017V5/Run2017_Mu_ReMiniAOD/Run2017_Mu_ReMiniAOD.submit
condor_submit CMSConnect_June2_2017V5/ST_s-c_4f_lep_PSw/ST_s-c_4f_lep_PSw.submit
condor_submit CMSConnect_June2_2017V5/ST_t-c_top_4f_inc/ST_t-c_top_4f_inc.submit
condor_submit CMSConnect_June2_2017V5/ST_tW_antitop_5f_inc/ST_tW_antitop_5f_inc.submit
condor_submit CMSConnect_June2_2017V5/ST_tW_top_5f_inc_PSw/ST_tW_top_5f_inc_PSw.submit
condor_submit CMSConnect_June2_2017V5/TT_AllHadronic/TT_AllHadronic.submit
condor_submit CMSConnect_June2_2017V5/TT_DiLep/TT_DiLep.submit
condor_submit CMSConnect_June2_2017V5/TT_SingleLep/TT_SingleLep.submit
condor_submit CMSConnect_June2_2017V5/WJets-HT100To200/WJets-HT100To200.submit
condor_submit CMSConnect_June2_2017V5/WJets-HT1200To2500/WJets-HT1200To2500.submit
condor_submit CMSConnect_June2_2017V5/WJets-HT200To400/WJets-HT200To400.submit
condor_submit CMSConnect_June2_2017V5/WJets-HT400To600/WJets-HT400To600.submit
condor_submit CMSConnect_June2_2017V5/WJets-HT600To800/WJets-HT600To800.submit
condor_submit CMSConnect_June2_2017V5/WJets-HT800To1200/WJets-HT800To1200.submit
condor_submit CMSConnect_June2_2017V5/WJets_madgraph/WJets_madgraph.submit
condor_submit CMSConnect_June2_2017V5/WW/WW.submit
condor_submit CMSConnect_June2_2017V5/WW_1L1Nu2Q/WW_1L1Nu2Q.submit
condor_submit CMSConnect_June2_2017V5/WW_LNu2Q_nlo/WW_LNu2Q_nlo.submit
condor_submit CMSConnect_June2_2017V5/WZ/WZ.submit
condor_submit CMSConnect_June2_2017V5/ZJetsToNuNu_HT100To200/ZJetsToNuNu_HT100To200.submit
condor_submit CMSConnect_June2_2017V5/ZJetsToNuNu_HT200To400/ZJetsToNuNu_HT200To400.submit
condor_submit CMSConnect_June2_2017V5/ZJetsToNuNu_HT2500ToInf/ZJetsToNuNu_HT2500ToInf.submit
condor_submit CMSConnect_June2_2017V5/ZJetsToNuNu_HT400To600/ZJetsToNuNu_HT400To600.submit
condor_submit CMSConnect_June2_2017V5/ZJetsToNuNu_HT600To800/ZJetsToNuNu_HT600To800.submit
condor_submit CMSConnect_June2_2017V5/ZJetsToNuNu_HT800To1200/ZJetsToNuNu_HT800To1200.submit
condor_submit CMSConnect_June2_2017V5/ZZ/ZZ.submit
condor_submit CMSConnect_June2_2017V5/ZZ_4L/ZZ_4L.submit
```

As the jobs are submitting, sometimes it will hang just after a submission. A simple `ctrl + C` will produce a traceback complaining about some dashboard script called by condor\_submit being interrupted, but then it will continue on with job submission.

Speaking of the dashboard, to monitor the status of your jobs and view rich information about them, you can visit the dashboard at dashb-cms-job-task.cern.ch and search for your name to see your jobs.

## Utility Function

Included in the AnalysisTools/scripts directory is a command line utility named `cms-connect-helper` which takes the job submission directory and a subcommand and does something useful for you. Note that the helper script must be run within the VHbbAnalysis directory.

### Job Status

If you want a text print out that you could pipe into a log file, then use the **status** subcommand.

```
[<username>@login VHbbAnalysis]$ ./../scripts/cms-connect-helper CMSConnect_June2_2017V5 status
Initiating submission status check at Sat Jun 02 19:41:22 CDT 2018

Summary for Run2017_Ele_ReMiniAOD
--------------------------------------------------
Number of Queued Jobs: 0
Number of Running Jobs: 1683
Number of Completed Jobs: 607
 |- Number of Successes (Exit Code == 0): 604
 |- Number of Failures (Exit Code != 0): 3
Number of Jobs with Unknown Status 0

Summary for Run2017_Mu_ReMiniAOD
--------------------------------------------------
Number of Queued Jobs: 4094
Number of Running Jobs: 2836
Number of Completed Jobs: 5
 |- Number of Successes (Exit Code == 0): 4
 |- Number of Failures (Exit Code != 0): 1
Number of Jobs with Unknown Status 0
```

### Job Resubmission

Failed jobs are inevitable, but can be dealt with by the **resubmit** command. The implementation parses all available log files for job status updates indicating that a job has terminated with a non-zero return value. Once those failed have been determined, their corresponding submit files are parsed to retrieve their arguments. Finally, a new job resubmission area and file is created within the original submission directory.

```bash
[<username>@login VHbbAnalysis]$ ./../scripts/cms-connect-helper CMSConnect_June2_2017V5 resubmit
Initiating resubmission check at Sat Jun 02 19:58:33 CDT 2018
A resubmission job file combining the arguments of all
currently failed jobs under CMSConnect_June2_2017V5 has been generated and
is ready for submission:

condor_submit CMSConnect_June2_2017V5/RESUBMIT/re.submit

The log files will be available under CMSConnect_June2_2017V5/RESUBMIT.

```

By design, this command does not automatically submit the resubmission job, opting to give the user a chance to check things over before submission. Also, it will only create a resubmission file if one does not already exist to prevent accidental overwrites. One last caveat, resubmitting in this manner assumes the X509\_USER\_PROXY environment variable points to your valid voms-proxy.


# Job Submission Campaigns

## June 9, 2016V4 Analysis with Factorized JECs and CMVA
```bash
# Data
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 2 -s Run2016BToG_Ele_ReMiniAOD
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 2 -s Run2016BToG_Mu_ReMiniAOD
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 2 -s Run2016BToG_MET_ReMiniAOD
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 2 -s Run2016BToG_DoubleEle_ReMiniAOD
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 2 -s Run2016BToG_DoubleMu_ReMiniAOD

# Signals
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.04 -s ZH125_ZNuNu_powheg
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s ggZH125_ZNuNu_powheg
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s WminusH125_powheg
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s WplusH125_powheg
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s ZH125_ZLL_powheg
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s ggZH125_ZLL_powheg

# W+Jets
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.06 -s WJets_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s WJets-HT100To200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s WJets-HT200To400
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s WJets-HT400To600
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.03 -s WJets-HT600To800
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.03 -s WJets-HT800To1200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s WJets-HT1200To2500
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s WJets-HT2500ToInf
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s WBJets-Pt100To200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.015 -s WBJets-Pt200ToInf
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s WJets_BGenFilter-Pt100To200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.015 -s WJets_BGenFilter-Pt200ToInf

# DY+Jets
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s DYToLL_madgraph
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.015 -s DYToLL_HT100to200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.015 -s DYToLL_HT200to400
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s DYToLL_HT400to600
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s DYToLL_HT600to800
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s DYToLL_HT800to1200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.01 -s DYToLL_HT1200to2500
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s DYToLL_HT2500toInf
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s DYBJets-Pt100To200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s DYBJets-Pt200ToInf
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s DYBJets-Pt100To200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.015 -s DYBJets-Pt200ToInf
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.02 -s DYJets-BGenFilter-Pt100To200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.015 -s DYJets-BGenFilter-Pt200ToInf

# Z+Jets
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.125 -s ZJetsToNuNu_HT100To200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.03 -s ZJetsToNuNu_HT200To400
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.03 -s ZJetsToNuNu_HT400To600
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s ZJetsToNuNu_HT600To800
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.03 -s ZJetsToNuNu_HT800To1200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s ZJetsToNuNu_HT1200To2500
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.125 -s ZJetsToNuNu_HT2500ToInf
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s ZBJetsToNuNu_Pt-100to200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s ZBJetsToNuNu_Pt-200toInf
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s ZJetsToNuNu_BGenFilter_Pt-100to200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s ZJetsToNuNu_BGenFilter_Pt-200toInf

# QCD
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.125 -s QCD_HT100to200
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.125 -s QCD_HT200to300
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.125 -s QCD_HT300to500
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.125 -s QCD_HT500to700
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.125 -s QCD_HT700to1000
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.125 -s QCD_HT1000to1500
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.125 -s QCD_HT1500to2000
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.125 -s QCD_HT2000toInf

# Single Top
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.05 -s TToLeptons_s
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.03 -s TToLeptons_t_powheg
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.05 -s TBarToLeptons_t_powheg
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.03 -s T_tW
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.05 -s Tbar_tW

# TTbar
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.03 -s TT_powheg

# Diboson
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.03 -s WW
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.04 -s WZ
python RunCMSConnect.py --submitJobs 0 --runOnSkim 0 --doKinFit 1 -b 1 -n CMSConnect_June9_2016V4 -c vhbb_config_2016.txt -f 0.025 -s ZZ
```

