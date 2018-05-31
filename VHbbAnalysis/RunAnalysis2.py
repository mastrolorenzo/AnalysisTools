#!/usr/bin/env python
import math
import optparse
import os
import socket
import sys
import time

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import ReadInput


# Additional input files to transfer to the remote execution host.
INPUTS_TO_TRANSFER = [
    'RunSample.py',
    'aux',
    'cfg',
    '../AnalysisDict.so',
    '../AnalysisDict_rdict.pcm',
    '../AnalysisManager.h',
    '../HelperClasses',
    '../env.sh',
    '../plugins',
    '../python/ReadInput.py',
    '../python/kinfitter.py',
]

# The template for the HTCondor submit description.
SUBMIT_DESCRIPTION = """\
# Automatically generated on {timestamp}
universe = vanilla
should_transfer_files = YES

executable = {executable}
arguments = {configFile} {sampleName} {filesToRun} output_{sampleName}_{nProcJobs}.root {start_event_frac} {end_event_frac} {RunSample_args}
initialdir = {jobName}/{sampleName}
transfer_input_files = {inputs_to_transfer}
transfer_output_files = ""
output = {nProcJobs}.stdout
error = {nProcJobs}.stderr
log = {nProcJobs}.log
notification = never

{extraCommands}

queue
"""

# The template for the batch exectuable at FNAL and CERN.
RUNSCRIPT = """\
#!/usr/bin/env bash
# Automatically generated on {timestamp}

echo "setting up the environment"
export SCRAM_ARCH="{scram_arch}"
source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project CMSSW {cmssw_version}
cd {cmssw_version}/src
eval "$(scramv1 runtime -sh)"
cd "$_CONDOR_SCRATCH_DIR"
echo "successfully set up the enviroment"

{copy_string}

echo "running RunSample.py"
python RunSample.py $1 $2 $3 "$_CONDOR_SCRATCH_DIR"/$4 $5 $6 $7
echo "done running, now copying output to EOS"
{xrdcp_string}
echo "all done!"
"""

# The template for the batch exectuable on CMSConnect.
RUNSCRIPT_CMSCONNECT = """\
#!/usr/bin/env bash
# Automatically generated on {timestamp}
export XrdSecGSISRVNAMES=cmseos.fnal.gov
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch

echo "setting up the environment"
export SCRAM_ARCH="{scram_arch}"
source $VO_CMS_SW_DIR/cmsset_default.sh
scram project CMSSW {cmssw_version}
cd {cmssw_version}/src
eval "$(scramv1 runtime -sh)"
cd "$_CONDOR_SCRATCH_DIR"
echo "successfully set up the enviroment"

{copy_string}

echo "running RunSample.py"
python RunSample.py $1 $2 $3 "$_CONDOR_SCRATCH_DIR"/$4 $5 $6 $7
echo "done running, now copying output to EOS"
{xrdcp_string}
echo "all done!"
"""


# The template for the batch executable at DESY.
RUNSCRIPT_DESY = """\
export PATH=/afs/desy.de/common/passwd:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/cvmfs/grid.cern.ch/emi3ui-latest/bin:/cvmfs/grid.cern.ch/emi3ui-latest/sbin:/cvmfs/grid.cern.ch/emi3ui-latest/usr/bin:/cvmfs/grid.cern.ch/emi3ui-latest/usr/sbin:$PATH
echo "echo PATH:"
echo $PATH
echo "arguments: " $1 $2 $3 $4 $5 $6 $7
echo "username and group"
id -n -u
id -n -g

echo "creating tempdir and copy"
tmp_dir=$(mktemp -d)
cp -r {inputs_to_transfer} $tmp_dir

echo "setting up the environment"
export SCRAM_ARCH="{scram_arch}"
source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project CMSSW {cmssw_version}
cd {cmssw_version}/src
eval "$(scramv1 runtime -sh)"

echo "echo PATH:"
echo $PATH

echo "changing to tempdir"
cd $tmp_dir
pwd
ls

echo "running RunSample.py"
echo $4
python RunSample.py $1 $2 $3 $4 $5 $6 $7
echo "done running, now copying output to EOS"

echo "copying output ({xrdcp_string})"
{xrdcp_string}

echo "delete tmp dir"
cd $TMP
rm -r $tmp_dir

echo "all done!"
"""


def parse_command_line(argv):
    """Parse user settings from the command line."""
    if argv is None:
        argv = sys.argv[1:]

    parser = optparse.OptionParser()
    parser.add_option("-c", "--config", dest="configFile", default="ewk_config.txt", help="specify config file for this analysis")
    parser.add_option("-b", "--batch", dest="runBatch", default=0, type=int, help="run analysis jobs on condor (default=0)")
    parser.add_option("-n", "--jobName", dest="jobName", default="condor_jobs", help="Specify label for condor jobs. Only to be used when running batch jobs")
    parser.add_option("-f", "--nFilesPerJob", dest="nFilesPerJob", default=10, type=float, help="Number of input files per batch job")
    parser.add_option("-o", "--outputDir", dest="outputDir", default="", type=str, help="Output directory for the jobs")
    parser.add_option("-s", "--sample", dest="sample", default="", type=str, help="Run on only a specific sample (can be comma-separated list)")
    parser.add_option("-d", "--doData", dest="doData", default=-1, type=int, help="If -1 run all samples, if 0 run only MC, if 1 run only data")
    parser.add_option("--site", "--site", dest="site", default="FNAL", type=str, help="If running on lxplus include option --site CERN, otherwise assumes you are running on FNAL")
    parser.add_option("--useSGE", "--useSGE", dest="useSGE", default=0, type=int, help="If 0 (default) use condor, if 1 use SGE job submission")
    parser.add_option("--doSkim", "--doSkim", dest="doSkim", default=0, type=int, help="If 0 (default) run analysis jobs, if 1 run skimming jobs")
    parser.add_option("--runOnSkim", "--runOnSkim", dest="runOnSkim", default=1, type=int, help="If 1 (default) run analysis jobs assuming the input files are already skimmed by AT, if 0 assume running directly on post-processed NanoAOD. If doSkim is 1 then this variable is assumed to be always 0.")
    parser.add_option("--doKinFit", "--doKinFit", dest="doKinFit", default=0, type=int, help="If 1 (default=0) runs the kinematic fit after the analysis job (Automatically turned off when doSkim is 1).")
    parser.add_option("--submitJobs", "--submitJobs", dest="submitJobs", default=1, type=int, help="If 1 (default) submit jobs to batch queue, if 0 then only create submission files")

    if len(argv) == 1:
        print parser.format_help()
        sys.exit(-1)

    options, args = parser.parse_args(argv)

    return options, args


def main(argv=None):
    ROOT.gROOT.SetBatch(True)

    options, _ = parse_command_line(argv)

    # reads samples, existing branches and new branches
    if options.sample != "":
        samplesToSubmit = options.sample.split(',')
    else:
        samplesToSubmit = [] # if empty run on all samples

    configFile = options.configFile
    jobName = options.jobName
    doData = options.doData
    site = options.site
    useSGE = options.useSGE
    submitJobs = options.submitJobs

    if options.doSkim:
        options.doKinFit = False
        options.runOnSkim = False

    # Setup the analysis manager.
    analysis_manager = ReadInput.ReadTextFile(
        filename=configFile,
        filetype="cfg",
        samplesToRun=samplesToSubmit,
        filesToRun="",
        isBatch=options.runBatch,
        doSkim=options.doSkim,
        runOnSkim=options.runOnSkim,
    )

    analysis_manager.debug=2

    # Run the analysis.
    if not options.runBatch:
        print "Running locally over all samples"
        if analysis_manager.debug > 100:
            analysis_manager.PrintBranches()
        # loop over all the samples
        # FIXME - need to add the possibility of doing a small portion of files
        analysis_manager.Loop()
    else:
        # If the site is the default value, try to autodetect the correct site.
        if site == "FNAL":
            hostname = socket.gethostname()
            if hostname.endswith('.desy.de'):
                print 'detected site: DESY'
                site = 'DESY'
            if hostname.endswith('.cern.ch'):
                print 'detected site: CERN'
                site = 'CERN'

        if site not in {"FNAL", "CERN", "DESY", "CMSCONNECT"}:
            print "unknown site: %s" % site
            sys.exit(-1)

        print "Running analysis jobs on batch queue. Site: %s" % site

        if useSGE:
            print "Using SGE batch system..."
        else:
            print "Using Condor batch system..."

        if site == "FNAL":
            # parent directory in user's EOS space to store all output files
            output_dir = "/eos/uscms/store/user/%s/VHbbAnalysisNtuples" % os.getlogin()
            # parent directory in group's EOS space to store all output files
            #output_dir = "/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples/"
        elif site == "CERN":
            output_dir = "/eos/cms/store/user/%s/VHbbAnalysisNtuples" % os.getlogin()
        elif site == "DESY":
            output_dir = "/nfs/dust/cms/user/%s/VHbbAnalysisNtuples" % os.getlogin()
        elif site == "CMSCONNECT":
            output_dir = "/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples"

        if options.doSkim:
            output_dir = output_dir.replace("VHbbAnalysisNtuples", "SkimmedAnalysisNtuples")

        if options.outputDir:
            print "changing output_dir to %s" % options.outputDir
            output_dir = options.outputDir

        os.system("mkdir -p %s" % jobName)
        if site == "CMSCONNECT":
            start = output_dir.find('/store')
            os.system("xrdfs root://cmseos.fnal.gov mkdir -p %s/%s" % (output_dir[start:], jobName))
        else:
            os.system("mkdir -p %s/%s" % (output_dir, jobName))

        submitFiles = []

        inputs_to_transfer = [configFile]
        inputs_to_transfer.extend(INPUTS_TO_TRANSFER)

        # copy inputs to jobdir
        if site == "DESY":
            os.system('cp -r %s %s' % (' '.join(inputs_to_transfer), jobName))
            # prepend paths with '../' for the condor jobs
            inputs_to_transfer = list('../'+p if not p.startswith('../') else p for p in inputs_to_transfer)
            inputs_to_transfer.remove('../python/ReadInput.py')
            inputs_to_transfer.remove('../python/kinfitter.py')
            inputs_to_transfer.append('../ReadInput.py')
            inputs_to_transfer.append('../kinfitter.py')
        else:
            inputs_to_transfer = list('../../'+p for p in inputs_to_transfer)

        nFilesPerJob = options.nFilesPerJob
        for sample in analysis_manager.samples:
            if (options.sample != ""):
                if (sample.sampleName not in samplesToSubmit):
                    print "sample: ",sample.sampleName," not in list, skipping..."
                    continue
            if (options.doData == 0 and sample.sampleName.find("Run") != -1):
                print "skipping data sample: ",sample.sampleName
                continue
            if (options.doData == 1 and sample.sampleName.find("Run") == -1):
                print "skipping MC sample: ",sample.sampleName
                continue
            sampleName = sample.sampleName
            print sampleName
            os.system("mkdir -p %s/%s" % (jobName, sampleName))
            if site == "CMSCONNECT":
                start = output_dir.find('/store')
                os.system("xrdfs root://cmseos.fnal.gov mkdir -p %s/%s/%s" % (output_dir[start:], jobName, sampleName))
            else:
                os.system("mkdir -p %s/%s/%s" % (output_dir, jobName, sampleName))
            nProcJobs = 0
            nFiles = len(sample.files)
            if nFilesPerJob < 1:
                if not math.ceil(1. / nFilesPerJob) == 1. / nFilesPerJob:
                    nFilesPerJob = 1 / (math.ceil(1. / nFilesPerJob))
            nJobs = int((nFiles // nFilesPerJob) + 1)
            start_event_frac = 0.0
            end_event_frac = 1.0
            for i in range(nJobs):
                filesToRun = ""
                if nFilesPerJob >= 1:
                    for j in range(int(nFilesPerJob)):
                        index = i*nFilesPerJob + j
                        if (index >= nFiles):
                            continue
                        filesToRun += "%s," % sample.files[int(index)]
                else:
                    index = math.floor(i * nFilesPerJob)
                    start_event_frac = i*nFilesPerJob - math.floor(i*nFilesPerJob)
                    end_event_frac = start_event_frac + nFilesPerJob
                    if (index >= nFiles):
                        continue
                    filesToRun += "%s," % sample.files[int(index)]
                filesToRun = filesToRun[:-1] # remove trailing ','
                if filesToRun == "":
                    continue

                nProcJobs += 1
                RunSample_args = "runOnSkim" if options.runOnSkim else ""
                RunSample_args += ",doSkim" if options.doSkim else ""
                RunSample_args += ",doKinFit" if options.doKinFit else ""

                # Additional site-based HTCondor comamnds.
                extraCommands = {}
                if site == 'CERN':
                    extraCommands['+MaxRuntime'] = '14400'

                if not useSGE:
                    submitPath = "%s/%s/job%i.submit" % (jobName, sampleName, nProcJobs)
                    if site == 'CMSCONNECT':
                        executable = os.path.abspath("%s/condor_runscript.sh" % jobName)
                    else:
                        executable = "%s/condor_runscript.sh" % jobName
                    content = SUBMIT_DESCRIPTION.format(
                        timestamp=time.strftime('%a %b %d %H:%M:%S %Z %Y'),
                        executable=executable,
                        jobName=jobName,
                        configFile=configFile,
                        sampleName=sampleName,
                        filesToRun=filesToRun,
                        nProcJobs=str(nProcJobs),
                        start_event_frac=str(start_event_frac),
                        end_event_frac=str(end_event_frac),
                        RunSample_args=RunSample_args,
                        inputs_to_transfer=','.join(inputs_to_transfer),
                        extraCommands='\n'.join('{0} = {1}'.format(*items) for items in extraCommands.iteritems()),
                    )
                else:
                    submitPath = "%s/%s/job%iSubmit.sh" % (jobName, sampleName, nProcJobs)
                    content = "source %s/%s/condor_runscript.sh %s %s %s output_%s_%i.root %f %f %s\n" % (os.getcwd(), jobName, configFile, sampleName, filesToRun, sampleName, nProcJobs, start_event_frac, end_event_frac, RunSample_args)

                with open(submitPath, "w") as submitFile:
                    submitFile.write(content)
                submitFiles.append(submitPath)

        # Create the scripts to execute on the remote host.
        copy_string = "" # not sure how to automatically transfer files to SGE job nodes, for now do it manually
        if useSGE:
            copy_string = ''' cp -r {0}/cfg .
            cp -r {0}/aux .
            cp {0}/RunSample.py .
            cp {0}/../AnalysisDict.so .
            cp -r {0}/*.txt . '''.format(os.getcwd())

        if site == "FNAL":
            xrdcp_string = "xrdcp -f $_CONDOR_SCRATCH_DIR/$4 %s/%s/$2" % (output_dir.replace("/eos/uscms", "root://cmseos.fnal.gov/"), jobName)
            runscript_template = RUNSCRIPT
        elif site == "CERN":
            if output_dir.startswith('/eos/user'):
                xrdcp_string = "xrdcp -f $_CONDOR_SCRATCH_DIR/$4 root://eosuser.cern.ch/%s/%s/$2" % (output_dir, jobName)
            else:
                xrdcp_string = "xrdcp -f $_CONDOR_SCRATCH_DIR/$4 root://eoscms.cern.ch/%s/%s/$2" % (output_dir, jobName)
            runscript_template = RUNSCRIPT
        elif site == "DESY":
            local_path = "%s/%s" % (output_dir, jobName)
            xrdcp_string = "mkdir -p %s/$2; cp -vf $4 %s/$2" % (local_path, local_path)
            runscript_template = RUNSCRIPT_DESY
        elif site == "CMSCONNECT":
            xrdcp_string = "xrdcp -f $_CONDOR_SCRATCH_DIR/$4 root://cmseos.fnal.gov/%s/%s/$2" % (output_dir, jobName)
            runscript_template = RUNSCRIPT_CMSCONNECT

        runscript_content = runscript_template.format(
            timestamp=time.strftime('%a %b %d %H:%M:%S %Z %Y'),
            scram_arch=os.environ['SCRAM_ARCH'],
            cmssw_version=os.environ['CMSSW_VERSION'],
            copy_string=copy_string,
            xrdcp_string=xrdcp_string,
            inputs_to_transfer=' '.join(inputs_to_transfer),
        )

        with open(os.path.join(jobName, 'condor_runscript.sh'), 'w') as runscript:
            runscript.write(runscript_content)
        os.system("chmod 755 %s/condor_runscript.sh" % (jobName))

        if submitJobs:
            if not useSGE:
                ## Send job to condor
                print "Submit files created, sending jobs to Condor..."
                for submitFile in submitFiles:
                    print("condor_submit %s" % submitFile)
                    os.system("condor_submit %s" % submitFile)
            else:
                ## Send jobs with SGE
                print "Submit files created, sending jobs to SGE batch..."
                for submitFile in submitFiles:
                    print 'bsub -R "pool>30000" -q 1nh -J job1 < %s' % submitFile
                    os.system('bsub -R "pool>30000" -q 1nh -J job1 < %s' % submitFile)
        else:
            print "Submit files created, but will not submit jobs to batch queue"


if __name__ == '__main__':

    status = main()
    sys.exit(status)

