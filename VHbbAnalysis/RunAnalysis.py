import ROOT
import sys
import os
import math
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config",       dest="configFile",  default="vhbb_config_2017.txt",   help="specify config file for this analysis" )
parser.add_argument("-b", "--batch",        dest="runBatch",    default=0, type=int,              help="run analysis jobs on condor (default=0)")
parser.add_argument("-n", "--jobName",      dest="jobName",     default="condor_jobs",            help="Specify label for condor jobs. Only to be used when running batch jobs")
parser.add_argument("-f", "--nFilesPerJob", dest="nFilesPerJob",default=10, type=float,           help="Number of input files per batch job")
parser.add_argument("-o", "--outputDir",    dest="outputDir",   default="", type=str,             help="Output directory for the jobs")
parser.add_argument("-s", "--sample",       dest="sample",      default="", type=str,             help="Run on only a specific sample (can be comma-separated list)")
parser.add_argument("-d", "--doData",       dest="doData",      default=-1, type=int,             help="If -1 run all samples, if 0 run only MC, if 1 run only data")
parser.add_argument("--site",               dest="site",        default="FNAL", type=str,         help="If running on lxplus include option --site CERN, otherwise assumes you are running on FNAL")
parser.add_argument("--sgl" ,          dest="oneJobPerCluster", default=False,action="store_true",help="Old style condor submission with one job per cluster")
parser.add_argument("--dcache",             dest="dcache",      default=False,action="store_true",help="Copy files to dcache at DESY instead of locally?")
parser.add_argument("--useSGE",             dest="useSGE",      default=0, type=int,              help="If 0 (default) use condor, if 1 use SGE job submission")
parser.add_argument("--doSkim",             dest="doSkim",      default=0, type=int,              help="If 0 (default) run analysis jobs, if 1 run skimming jobs")
parser.add_argument("--runOnSkim",          dest="runOnSkim",   default=0, type=int,              help="If 0 (default) run analysis jobs assuming the input files are already skimmed by AT, if 0 assume running directly on post-processed NanoAOD. If doSkim is 1 then this variable is assumed to be always 0.")
parser.add_argument("--onlyMvaEval",        dest="onlyMvaEval", default=0, type=int,              help="If 1 (default=0) runs the kinematic fit after the analysis job (Automatically turned off when doSkim is 1). THIS ONLY GIVES MEANINGFULL RESULT WITH -f >= 1!!")
parser.add_argument("--submitJobs",         dest="submitJobs",  default=1, type=int,              help="If 1 (default) submit jobs to batch queue, if 0 then only create submission files")
parser.add_argument("--kill_if_runtime_above_minutes",dest="kill_if_runtime_above_minutes",default=-1,type=int,help="After 100 events it's estimated how long the job would take. The job exits with code 71 if it would take longer than the given amount.")
args = parser.parse_args()

#This import needs to be after the parser or the help menu is messed up.
import ReadInput
#import ReadInputWA as ReadInput

if len(sys.argv) == 1:
    print parser.format_help()
    print "exiting",len(sys.argv)
    exit(-1)

# job workers need to compile their own version
os.system('rm aux/roccor.2016.v3_mod/RoccoR_cc*')  

ROOT.gSystem.Load("AnalysisDict.so")

# reads samples, existing branches and new branches
if (args.sample != ""):
    samplesToSubmit = args.sample.split(',')
else:
    samplesToSubmit = [] # if empty run on all samples

if args.doSkim:
    args.onlyMvaEval = False
    args.runOnSkim = False

if args.onlyMvaEval:
    assert args.nFilesPerJob >= 1, 'with onlyMvaEval, nFilesPerJob must be >=1'

ReadInput.debug = 2
am=ReadInput.ReadTextFile(args.configFile, "cfg", samplesToSubmit,"", args.doSkim, args.runOnSkim)
am.debug = 2


if (args.runBatch == False):
    print "Running locally over all samples"

    if(am.debug>100):
        am.PrintBranches()

    # loop over all the samples
    # FIXME - need to add the possibility of doing a small portion of files
    am.Loop()

else:
    print "Running analysis jobs on batch queue. Site: %s" % args.site
    if args.useSGE:
        print "Using SGE batch system..."
    else:
        print "Using Condor batch system..."
    jobName = args.jobName

    if args.site not in ("FNAL", "CERN", "DESY"):
        print "unknown site: ",args.site,"exiting..."
        sys.exit(1)

    if args.site == "FNAL":  ## let's do a little autodetect for desy and cern..
        import socket
        hostname = socket.gethostname()
        if hostname.endswith('.desy.de'):
            print 'detected site: DESY'
            args.site = 'DESY'
        if hostname.endswith('.cern.ch'):
            print 'detected site: CERN'
            args.site = 'CERN'

    if args.site == "FNAL":
        #output_dir = "/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples" # parent directory in user's EOS space to store all output files
        output_dir = "/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/" # parent directory in user's EOS space to store all output files
    elif args.site == "CERN":
        output_dir = "/eos/cms/store/user/%s/VHbbAnalysisNtuples" % os.getlogin()
    elif args.site == "DESY":
        if not args.dcache:
            output_dir = "/nfs/dust/cms/user/%s/VHbbAnalysisNtuples" % os.getlogin()
        else :
            if os.path.expandvars("$CERNUSERNAME") == '$CERNUSERNAME':
                print "No username set, exiting"
                sys.exit(1)
            output_dir = os.path.expandvars("srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN=/pnfs/desy.de/cms/tier2/store/user/$CERNUSERNAME/VHbbAnalysisNtuples/")
    if args.doSkim:
        output_dir = output_dir.replace("VHbbAnalysisNtuples","SkimmedAnalysisNtuples")

    if args.outputDir:
        print "changing output_dir to",args.outputDir
        output_dir = args.outputDir


    #if os.path.exists(jobName):
    #    print "Attempting to create jobs under jobName %s, which already exists!" % jobName
    #    sys.exit(0)
    os.system("mkdir -p %s" % jobName)
    os.system("mkdir -p %s/%s" % (output_dir,jobName))

    submitFiles = []

    #for sampleName in am.ListSampleNames():
    inputs_to_transfer = [
        args.configFile,
        'RunSample.py',
        '../AnalysisDict.so',
        '../env.sh',
        'cfg',
        'aux',
        '../python/ReadInput.py',
        '../python/kinfitter.py',
        '../python/MyStandardScaler.py',
        '../python/TensorflowEvaluator.py',
        '../python/TensorflowDNNClassifier.py',
        '../python/mva_evaluator.py',
        '../AnalysisDict_rdict.pcm',
        '../HelperClasses',
        '../AnalysisManager.h',
        '../plugins',
    ]
    # copy inputs to jobdir
    if args.site == "DESY":
        os.system('cp -r %s %s' % (' '.join(inputs_to_transfer), jobName))
        # prepend paths with '../' for the condor jobs
        inputs_to_transfer = list('../'+p if not p.startswith('../') else p for p in inputs_to_transfer)
        inputs_to_transfer.remove('../python/ReadInput.py')
        inputs_to_transfer.remove('../python/kinfitter.py')
        inputs_to_transfer.remove('../python/MyStandardScaler.py')
        inputs_to_transfer.remove('../python/TensorflowEvaluator.py')
        inputs_to_transfer.remove('../python/TensorflowDNNClassifier.py')
        inputs_to_transfer.remove('../python/mva_evaluator.py')
        inputs_to_transfer.append('../ReadInput.py')
        inputs_to_transfer.append('../kinfitter.py')
        inputs_to_transfer.append('../MyStandardScaler.py')
        inputs_to_transfer.append('../TensorflowEvaluator.py')
        inputs_to_transfer.append('../TensorflowDNNClassifier.py')
        inputs_to_transfer.append('../mva_evaluator.py')
    else:
        inputs_to_transfer = list('../../'+p for p in inputs_to_transfer)

    nFilesPerJob = args.nFilesPerJob
    for sample in am.samples:
        if (args.sample != ""):
            #if (sample.sampleName != args.sample): continue
            if (sample.sampleName not in samplesToSubmit ):
                print "sample: ",sample.sampleName," not in list, skipping..."
                continue
        if (args.doData == 0 and sample.sampleName.find("Run")!=-1):
            print "skipping data sample: ",sample.sampleName
            continue
        if (args.doData == 1 and sample.sampleName.find("Run")==-1):
            print "skipping MC sample: ",sample.sampleName
            continue
        sampleName = sample.sampleName
        print sampleName
        os.system("mkdir -p %s/%s" % (jobName,sampleName))
        os.system("mkdir -p %s/%s/%s" % (output_dir,jobName,sampleName))
        nProcJobs = 0
        nFiles = len(sample.files)
        if nFilesPerJob <1 :
            if not math.ceil(1./nFilesPerJob)==1./nFilesPerJob:
                nFilesPerJob=1/(math.ceil(1./nFilesPerJob))
        nJobs = int (( nFiles // nFilesPerJob ) + 1)
        #for filename in sample.files:
        start_event_frac=0
        end_event_frac=1
        arg_string=""
        for i in range(nJobs):
            filesToRun = ""
            if nFilesPerJob >= 1:
                for j in range(int(nFilesPerJob)):
                    index = i*nFilesPerJob + j
                    if (index >= nFiles): continue
                    filesToRun += "%s," % sample.files[int(index)]
            else :
                index=math.floor(i*nFilesPerJob)
                start_event_frac = i*nFilesPerJob-math.floor(i*nFilesPerJob)
                end_event_frac = start_event_frac+nFilesPerJob
                if (index >= nFiles): continue
                filesToRun += "%s," % sample.files[int(index)]
            filesToRun = filesToRun[:-1] # remove trailing ','
            if (filesToRun == ""): continue

            nProcJobs += 1
            RunSample_args = "runOnSkim" if args.runOnSkim else ""
            RunSample_args += ",doSkim" if args.doSkim else ""
            RunSample_args += ",onlyMvaEval" if args.onlyMvaEval else ""
            RunSample_args += ",kill_if_runtime_above_minutes=%i" % args.kill_if_runtime_above_minutes if args.kill_if_runtime_above_minutes > 0 else ""
            arg_string+="%s %s %s output_%s_%i.root %f %f %s\n"% (args.configFile, sampleName, filesToRun,sampleName, nProcJobs, start_event_frac, end_event_frac, RunSample_args)
            if args.useSGE:
                fname = "%s/%s/job%iSubmit.sh" % (jobName, sampleName,nProcJobs)
                submitFile = open(fname, "w")
                content = "source %s/%s/condor_runscript.sh %s %s %s output_%s_%i.root %f %f %s\n" % (os.getcwd(),jobName,args.configFile, sampleName, filesToRun,sampleName, nProcJobs, start_event_frac, end_event_frac, RunSample_args)
                submitFile.write(content)
                submitFile.close()
                submitFiles.append(fname)
            elif args.oneJobPerCluster:
                fname = "%s/%s/job%i.submit" % (jobName, sampleName,nProcJobs)
                submitFile = open(fname, "w")
                content =  "universe = vanilla\n"
                content += "Executable = %s/condor_runscript.sh\n" % jobName
                content += "Arguments = %s %s %s output_%s_%i.root %f %f %s\n" % (args.configFile, sampleName, filesToRun,sampleName, nProcJobs, start_event_frac, end_event_frac, RunSample_args)
                #content += "Arguments = %s %s %s %i\n" % (args.configFile, sampleName, filename, nProcJobs)
                #content += "Requirements   =  OpSys == 'LINUX' && (Arch =='INTEL' || Arch =='x86_64')\n"
                content += "initialdir = %s/%s\n" % (jobName,sampleName)
                #content += "Should_Transfer_Files = YES\n"
                content += "Output = %i.stdout\n" % nProcJobs
                content += "Error  = %i.stderr\n" % nProcJobs
                content += "Log    = %i.log\n"    % nProcJobs
                content += "Notification = never\n"
                content += "on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n"
                if args.dcache:
                    content += "x509userproxy = $ENV(X509_USER_PROXY)\n"
                #content += "WhenToTransferOutput=On_Exit\n"
                #content += "transfer_input_files = ../../%s,../../cfg/samples.txt,../../cfg/earlybranches.txt,../../cfg/existingbranches.txt,../../cfg/newbranches.txt,../../cfg/bdtsettings.txt,../../cfg/reg1_settings.txt,../../cfg/reg2_settings.txt,../../cfg/settings.txt,../../aux/TMVARegression_BDTG_ttbar_Nov23.weights.xml,../../aux/TMVA_13TeV_Dec14_3000_5_H125Sig_0b1b2bWjetsTTbarBkg_Mjj_BDT.weights.xml,../../aux/MuonIso_Z_RunCD_Reco74X_Dec1.json,../../aux/SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.json,../../aux/MuonID_Z_RunCD_Reco74X_Dec1.json,../../aux/CutBasedID_TightWP.json,../../aux/CutBasedID_LooseWP.json,../../RunSample.py,../../../AnalysisDict.so,../../cfg/systematics.txt,../../cfg/scalefactors.txt\n" % args.configFile
                content += "transfer_input_files = %s\n" % ','.join(inputs_to_transfer)
                content += "Queue  1\n"
                submitFile.write(content)
                submitFile.close()
                submitFiles.append(fname)
        if not args.oneJobPerCluster:
            fname = "%s/%s/job%s.submit" % (jobName, sampleName,sampleName)
            submitFile = open(fname, "w")
            content =  "universe = vanilla\n"
            content += "Executable = %s/condor_runscript.sh\n" % jobName
            content += "initialdir = %s/%s\n" % (jobName,sampleName)
            content += "Output = log-$(ClusterId).$(ProcId).stdout\n"
            content += "Error  = log-$(ClusterId).$(ProcId).stderr\n"
            content += "Log    = log-$(ClusterId).$(ProcId).log\n"
            content += "Notification = never\n"
            content += "on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n"
            if args.dcache:
                content += "x509userproxy = $ENV(X509_USER_PROXY)\n"
            content += "transfer_input_files = %s\n" % ','.join(inputs_to_transfer)
            if args.site == "FNAL":
                # Prevent condor from transferring files just because they are new,
                # e.g. Python bytecode, intermediate kinFit and MVA evaluation files.
                content += "transfer_output_files = \"\"\n"
            content += "Queue Arguments from (\n"
            content += "%s" %arg_string
            content += ")"
            submitFile.write(content)
            submitFile.close()
            submitFiles.append(fname)


    xrdcp_string = "xrdcp -f $ORIG_DIR/$4 %s/%s/$2" %(output_dir.replace("/eos/uscms","root://cmseos.fnal.gov/"),jobName)
    if args.site == "CERN":
        xrdcp_string = "xrdcp -f $ORIG_DIR/$4 root://eoscms.cern.ch:/%s/%s/$2" %(output_dir,jobName)
    if args.site == "DESY":
        local_path = "%s/%s" % (output_dir,jobName)
        if not args.dcache :
            xrdcp_string = "mkdir -p %s/$2; cp -vf $4 %s/$2" % (local_path, local_path)
        else:
            xrdcp_string = "gfal-copy -p file:$4 %s/$2/$4" % (local_path)
        if not args.dcache:
            xrdcp_string += '\n        python -c "import sys,ROOT; f=ROOT.TFile(\'%s/$2/$4\'); sys.exit(int(f.IsZombie() and 99))"' % local_path
    copy_string = "" # not sure how to automatically transfer files to SGE job nodes, for now do it manually
    if args.useSGE:
        copy_string = ''' cp -r {0}/cfg .
        cp -r {0}/aux .
        cp {0}/RunSample.py .
        cp {0}/../AnalysisDict.so .
        cp -r {0}/*.txt . '''.format(os.getcwd())

    #Set a few variables that we will need to identify the correct CMSSW release to source from cvmfs
    scram_arch = os.environ['SCRAM_ARCH']
    cmssw_version = os.environ['CMSSW_BASE'].split('/')[-1]
    patch_rel = 'cmssw-patch' if 'patch' in cmssw_version else 'cmssw'

    condor_runscript_text = '''

        #!/bin/bash
        ##
        ## Script to run Analysis Manager jobs on Condor from the LPC
        ## Author: Stephane Cooperstein
        ##
        ## Argument 1: Analysis config file
        ## Argument 2: Sample name to run on
        ##

        export ORIG_DIR=$PWD
        # Set up environment
        echo "setting up the environment"
        #cd %s/..
        cd /cvmfs/cms.cern.ch/%s/cms/%s/%s/src
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        eval `scramv1 runtime -sh`
        #source env.sh
        #source %s/../env.sh
        #source env.sh

        echo "successfully set up the enviroment"

        cd $ORIG_DIR
        pwd
        ls

        #echo "copying over necessary files"
        #cp -r %s/cfg .
        #cp -r %s/aux .
        #cp %s/RunSample.py .
        #cp %s/../AnalysisDict.so .
        #cp -r %s/*.txt .
        %s

        echo "running RunSample.py"
        echo $ORIG_DIR/$4
        python RunSample.py $1 $2 $3 $ORIG_DIR/$4 $5 $6 $7
        rc=$?
        if [[ $rc != 0 ]]
        then
            echo "got exit code from RunSample.py: " $rc
            exit $rc
        fi
        echo "done running, now copying output to EOS"
        ###xrdcp -f $ORIG_DIR/$4 root://cmseos.fnal.gov//store/user/sbc01/VHbbAnalysisNtuples/%s/$2
        %s
        rc=$?
        if [[ $rc != 0 ]]
        then
            echo "copy failed (either bad output from cp or file is Zombie)"
            exit $rc
        fi

        rm $ORIG_DIR/$4
        echo "all done!" ''' % (os.getcwd(),scram_arch,patch_rel,cmssw_version, os.getcwd(), os.getcwd(), os.getcwd(), os.getcwd(), os.getcwd(), os.getcwd(), copy_string, jobName, xrdcp_string )

    condor_runscript_text_desy = '''

        export PATH=/afs/desy.de/common/passwd:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/cvmfs/grid.cern.ch/emi3ui-latest/bin:/cvmfs/grid.cern.ch/emi3ui-latest/sbin:/cvmfs/grid.cern.ch/emi3ui-latest/usr/bin:/cvmfs/grid.cern.ch/emi3ui-latest/usr/sbin:$PATH
        echo "echo PATH:"
        echo $PATH
        echo "arguments: " $1 $2 $3 $4 $5 $6 $7
        echo "username and group"
        id -n -u
        id -n -g

        echo "creating tempdir and copy"
        tmp_dir=$(mktemp -d)
        cp -r %s $tmp_dir

        echo "setting up the environment"
        cd /cvmfs/cms.cern.ch/%s/cms/%s/%s/src
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        eval `scramv1 runtime -sh`
        echo "echo PATH:"
        echo $PATH
        source /cvmfs/grid.cern.ch/etc/profile.d/setup-cvmfs-ui.sh

        echo "changing to tempdir"
        cd $tmp_dir
        pwd
        ls

        echo "running RunSample.py"
        echo $4
        python RunSample.py $1 $2 $3 $4 $5 $6 $7
        rc=$?
        if [[ $rc != 0 ]]
        then
            echo "got exit code from RunSample.py: " $rc
            exit $rc
        fi
        echo "done running, now copying output to EOS"

        echo "copying output"
        %s
        rc=$?
        if [[ $rc != 0 ]]
        then
            echo "copy failed (either bad output from cp or file is Zombie)"
            exit $rc
        fi

        echo "delete tmp dir"
        cd $TMP
        rm -r $tmp_dir

        echo "all done!"
    ''' % (' '.join(inputs_to_transfer), scram_arch, patch_rel, cmssw_version, xrdcp_string)

    runscript = open("%s/condor_runscript.sh" % (jobName) , "w")
    runscript.write(condor_runscript_text_desy if args.site=='DESY' else condor_runscript_text)
    runscript.close()
    os.system("chmod 755 %s/condor_runscript.sh" % (jobName))

    if args.submitJobs:
        if not args.useSGE:
            ## Send job to condor
            print "Submit files created, sending jobs to Condor..."
            for file in submitFiles:
                print("condor_submit %s" % file)
                os.system("condor_submit %s" % file)
        else:
            ## Send jobs with SGE
            print "Submit files created, sending jobs to SGE batch..."
            for file in submitFiles:
                print 'bsub -R "pool>30000" -q 1nh -J job1 < %s' %file
                os.system('bsub -R "pool>30000" -q 1nh -J job1 < %s' %file)
    else:
        print "Submit files created, but will not submit jobs to batch queue"



