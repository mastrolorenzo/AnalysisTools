import os
import sys

#samples = ["Wj0b","Wj1b","Wj2b","Zj0b","Zj1b","Zj2b"]
samples = ["ggZH_hbb","ZH_hbb","WH_hbb","s_Top","TT","Wj0b","Wj1b","Wj2b","VVHF","VVLF","Zj0b","Zj1b","Zj2b","data_obs","QCD"]

cwd = os.getcwd()
ifile = open(sys.argv[1])
jobname = ""
doEWK = False
nSystPerJob = 10

samples_to_run = []
if (len(sys.argv)>2):
    jobname = sys.argv[2]
if (len(sys.argv)>3):
    samples_to_run = sys.argv[3].split(',')
#if (len(sys.argv)>4):
#    doEWK = bool(sys.argv[4])

#if doEWK:
#    samples = ["s_Top","TT","Wjets","VV","Zjets","data_obs","EWKWJets","IntEWKWJets","QCD_data"]

if len(samples_to_run) > 0 and samples_to_run[0]!="":
    ## run on sub-set of analysis samples
    samples = []
    for sample in samples_to_run:
        samples.append(sample)
print "samples = ",samples

scram_arch = os.environ['SCRAM_ARCH']
cmssw_version = os.environ['CMSSW_BASE'].split('/')[-1]
patch_rel = 'cmssw-patch' if 'patch' in cmssw_version else 'cmssw'

systematics = {}
if (len(sys.argv)>4):
    # prepare shape systematics
    sysfile = open(sys.argv[4])
    for line in sysfile:
        if (line[0] == '#'): continue
        line = line.strip()
        params = line.split(' ')
        paramsToKeep = []
        for param in params:
            if (param != ''):
                paramsToKeep.append(param)
        #print paramsToKeep
        if (paramsToKeep[1] != "shape"): continue
        sysSamples = paramsToKeep[3].split(',')
        # map systematic name to the corresponding weight in the ntuple and possibly an extension if it has a different BDT/Mjj output than the nominal
        # In that case, the BDT name used is [bdtname]_[paramsToKeep[5]]UpDown
        if (len(paramsToKeep) > 5):
            systematics[paramsToKeep[0]] = (paramsToKeep[4], paramsToKeep[5], sysSamples)
        else:
            systematics[paramsToKeep[0]] = (paramsToKeep[4], "", sysSamples)
    
systematics["nominal"] = (1.0,"",samples)
#print systematics
 
nToRun = (len(systematics.keys())/nSystPerJob + 1)
print sys.argv[1],sys.argv[2]
print jobname
nJobs = 0
for line in ifile:
    print line
    for sample in samples:
        nSys = 0
        for i in xrange(nToRun) :
            systematics_to_run = ""
            j = 0
            for systematic in systematics:
                if sample not in systematics[systematic][2]: continue # systematic doesn't apply to this sample
                if (j < (i+1)*nSystPerJob and j >= i*nSystPerJob):
                    systematics_to_run += systematic + ","
                j = j+1 
            systematics_to_run = systematics_to_run.strip(",")
            #print systematics_to_run
            cmd = line.replace("SAMP",sample)
            #print systematics[systematic]
            cmd = cmd.replace("SYST",systematics_to_run)
            cmd = cmd.replace("ISYS",str(i))
            submit_text = '''universe = vanilla
Executable         = condor_runscript_%s_%i_%s_%i.sh
Should_Transfer_Files         = YES
Output         = stdout_%s_%i_%s_%i
Error          = stderr_%s_%i_%s_%i
Log            = log_%s_%i_%s_%i
Notification         = never
transfer_input_files     = ../../splitSamples.py,../../systematics_ewk_Wmn.txt,../../systematics_ewk_Wen.txt,../../systematics_Wmn_demo.txt,../../systematics_Wlnshapes2016.txt,../../systematics_Znnshapes2016.txt,../../systematics_Zllshapes2016.txt,../../systematics_Wlnshapes2017.txt,../../systematics_Znnshapes2017.txt,../../systematics_Zllshapes2017.txt,../../../python/nano_samples_2016.py,../../../python/nano_samples_znn2016.py,../../../python/nano_samples_znn2017.py,../../../python/nano_samples_2017.py
WhenToTransferOutput=On_Exit
Queue          1
            ''' % (jobname,nJobs,sample,nSys,jobname,nJobs,sample,nSys,jobname,nJobs,sample,nSys,jobname,nJobs,sample,nSys)
            runscript_text = '''export ORIG_DIR=$PWD
#cd         %s
cd     /cvmfs/cms.cern.ch/%s/cms/%s/%s/src
source         /cvmfs/cms.cern.ch/cmsset_default.sh
eval         `scramv1 runtime -sh`
cd     $ORIG_DIR
     \n''' % (cwd,scram_arch,patch_rel,cmssw_version)
            runscript_text += cmd
            submitfile = open("submit_%s_%i_%s_%i.txt" % (jobname,nJobs,sample,nSys), "w")
            runscriptfile = open("condor_runscript_%s_%i_%s_%i.sh" % (jobname,nJobs,sample,nSys), "w")
            submitfile.write(submit_text)
            runscriptfile.write(runscript_text)
            submitfile.close()
            runscriptfile.close()
            nSys = nSys + 1
    nJobs += 1
print "submitting %i jobs times %i samples times %i systematics" % (nJobs,len(samples),len(systematics))
for i in range(nJobs):
    for sample in samples:
        nSys = 0
        for j in xrange(nToRun):
            print "condor_submit submit_%s_%i_%s_%i.txt" %(jobname,i,sample,j)
            os.system("condor_submit submit_%s_%i_%s_%i.txt" %(jobname,i,sample,j))
