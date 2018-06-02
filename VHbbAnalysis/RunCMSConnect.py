#!/usr/bin/env python
import math
import optparse
import os
import stat
import sys
import time

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from jinja2 import Template

import ReadInput


# Arguments specifying how to tar the job input files.
TAR_ARGUMENTS = [
    'RunSample.py',
    'cfg/',
    'aux/',
    '../AnalysisManager.h',
    '../AnalysisDict.so',
    '../AnalysisDict_rdict.pcm',
    '../HelperClasses',
    '../plugins',
    '--directory=../python',
    'ReadInput.py',
    'kinfitter.py',
    'mva_evaluator.py',
    'MyStandardScaler.py',
    'TensorflowDNNClassifier.py',
    'TensorflowEvaluator.py',

]

# The submit description template.
SUBMIT = Template("""\
# Automatically generated on {{ timestamp }}
universe = vanilla
should_transfer_files = YES
notification = never
x509userproxy = $ENV(X509_USER_PROXY)

executable = {{ scriptpath }}
initialdir = {{ submission }}/{{ sample }}
transfer_input_files = ../../{{ submission }}.tar.gz
transfer_output_files = ""
output = $(Process).stdout
error = $(Process).stderr
log = $(Process).log

{% for args in arguments %}
arguments = {{ args }}
queue
{% endfor %}
""")

# The executable script template.
RUNSCRIPT = Template("""\
#!/bin/bash

# Automatically generated on {{ timestamp }}

CONDOR_EXEC="$(basename $0)"

echo "$(date) - $CONDOR_EXEC - INFO - Setting up the environment"
export XrdSecGSISRVNAMES=cmseos.fnal.gov
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch

# Check if the GRID proxy was transferred.
voms-proxy-info

# Setup the CMS software environment.
export SCRAM_ARCH="{{ environ['SCRAM_ARCH'] }}"
source $VO_CMS_SW_DIR/cmsset_default.sh

# Checkout the CMSSW release and set the runtime environment. These
# commands are often invoked by their aliases "cmsrel" and "cmsenv".
scram project CMSSW {{ environ['CMSSW_VERSION'] }}
pushd {{ environ['CMSSW_VERSION'] }}/src
eval "$(scramv1 runtime -sh)"
popd

# Untar the input files.
tar -zxvf "$_CONDOR_SCRATCH_DIR/{{ jobname }}.tar.gz"

echo "$(date) - $CONDOR_EXEC - INFO - Running RunSample.py"
python RunSample.py $1 $2 $3 $4 $5 $6 $7

echo "$(date) - $CONDOR_EXEC - INFO - Transferring Output File (10 retries)"
for i in {1..10}; do
    {{ xrdcp_command }}
    XRDEXIT=$?
    if [[ $XRDEXIT -ne 0 ]]; then
        echo "xrdcp failed with exit code $XRDEXIT"
    fi
done
""")


def parse_command_line(argv):
    """Parse user settings from the command line."""
    if argv is None:
        argv = sys.argv[1:]

    parser = optparse.OptionParser()
    parser.add_option('-c', '--config', dest='configFile', default='ewk_config.txt', help='specify config file for this analysis')
    parser.add_option('-b', '--batch', dest='runBatch', default=0, type=int, help='run analysis jobs on condor (default=0)')
    parser.add_option('-n', '--jobName', dest='jobName', default='condor_jobs', help='Specify label for condor jobs. Only to be used when running batch jobs')
    parser.add_option('-f', '--nFilesPerJob', dest='nFilesPerJob', default=10, type=float, help='Number of input files per batch job')
    parser.add_option('-o', '--outputDir', dest='outputDir', default='', type=str, help='Output directory for the jobs')
    parser.add_option('-s', '--sample', dest='sample', default='', type=str, help='Run on only a specific sample (can be comma-separated list)')
    parser.add_option('-d', '--doData', dest='doData', default=-1, type=int, help='If -1 run all samples, if 0 run only MC, if 1 run only data')
    parser.add_option('--site', '--site', dest='site', default='FNAL', type=str, help='If running on lxplus include option --site CERN, otherwise assumes you are running on FNAL')
    parser.add_option('--useSGE', '--useSGE', dest='useSGE', default=0, type=int, help='If 0 (default) use condor, if 1 use SGE job submission')
    parser.add_option('--doSkim', '--doSkim', dest='doSkim', default=0, type=int, help='If 0 (default) run analysis jobs, if 1 run skimming jobs')
    parser.add_option('--runOnSkim', '--runOnSkim', dest='runOnSkim', default=1, type=int, help='If 1 (default) run analysis jobs assuming the input files are already skimmed by AT, if 0 assume running directly on post-processed NanoAOD. If doSkim is 1 then this variable is assumed to be always 0.')
    parser.add_option('--doKinFit', '--doKinFit', dest='doKinFit', default=0, type=int, help='If 1 (default=0) runs the kinematic fit after the analysis job (Automatically turned off when doSkim is 1).')
    parser.add_option('--submitJobs', '--submitJobs', dest='submitJobs', default=1, type=int, help='If 1 (default) submit jobs to batch queue, if 0 then only create submission files')

    if len(argv) == 1:
        print parser.format_help()
        sys.exit(-1)

    options, args = parser.parse_args(argv)

    return options, args


def main(argv=None):
    ROOT.gROOT.SetBatch(True)

    options, _ = parse_command_line(argv)

    # reads samples, existing branches and new branches
    if options.sample != '':
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
        filetype='cfg',
        samplesToRun=samplesToSubmit,
        filesToRun='',
        isBatch=options.runBatch,
        doSkim=options.doSkim,
        runOnSkim=options.runOnSkim,
    )

    analysis_manager.debug = 2

    # Run the analysis.
    if not options.runBatch:
        print 'Running locally over all samples'
        if analysis_manager.debug > 100:
            analysis_manager.PrintBranches()
        # loop over all the samples
        # FIXME - need to add the possibility of doing a small portion of files
        analysis_manager.Loop()
    else:
        if 'X509_USER_PROXY' not in os.environ:
            raise RuntimeError('Missing environment variable X509_USER_PROXY, please set it to the path of a valid GRID proxy')
        output_dir = '/eos/uscms/store/group/lpchbb/VHbbAnalysisNtuples'
        if options.doSkim:
            output_dir = output_dir.replace('VHbbAnalysisNtuples', 'SkimmedAnalysisNtuples')
        if options.outputDir:
            print 'changing output_dir to %s' % options.outputDir
            output_dir = options.outputDir

        os.system('mkdir -p %s' % jobName)

        # Generate tarball of input files.
        os.system('tar czf {0}.tar.gz {1} {2}'.format(jobName, configFile, ' '.join(TAR_ARGUMENTS)))

        # Generate the job executable.
        scriptFile = os.path.join(jobName, 'condor_runscript.sh')
        with open(scriptFile, 'w') as f:
            f.write(RUNSCRIPT.render(
                timestamp=time.strftime('%a %b %d %H:%M:%S %Z %Y'),
                environ=os.environ,
                jobname=jobName,
                xrdcp_command='xrdcp -f $4 root://cmseos.fnal.gov/{0}/{1}/$2'.format(output_dir, jobName)
            ))
        current_permissions = os.stat(scriptFile)
        os.chmod(scriptFile, current_permissions.st_mode | stat.S_IEXEC)

        submitFiles = []

        nFilesPerJob = options.nFilesPerJob
        for sample in analysis_manager.samples:
            if options.sample != '':
                if sample.sampleName not in samplesToSubmit:
                    print 'sample: ', sample.sampleName, ' not in list, skipping...'
                    continue
            if options.doData == 0 and sample.sampleName.find('Run') != -1:
                print 'skipping data sample: ', sample.sampleName
                continue
            if options.doData == 1 and sample.sampleName.find('Run') == -1:
                print 'skipping MC sample: ', sample.sampleName
                continue

            sampleName = sample.sampleName
            print sampleName
            os.system('mkdir -p %s/%s' % (jobName, sampleName))
            start = output_dir.find('/store')
            os.system('xrdfs root://cmseos.fnal.gov mkdir -p %s/%s/%s' % (output_dir[start:], jobName, sampleName))

            arguments = []

            nFiles = len(sample.files)
            if nFilesPerJob < 1:
                if not math.ceil(1. / nFilesPerJob) == 1. / nFilesPerJob:
                    nFilesPerJob = 1 / (math.ceil(1. / nFilesPerJob))
            nJobs = int((nFiles // nFilesPerJob) + 1)
            start_event_frac = 0.0
            end_event_frac = 1.0
            for i in range(nJobs):
                filesToRun = ''
                if nFilesPerJob >= 1:
                    for j in range(int(nFilesPerJob)):
                        index = i*nFilesPerJob + j
                        if (index >= nFiles):
                            continue
                        filesToRun += '%s,' % sample.files[int(index)]
                else:
                    index = math.floor(i * nFilesPerJob)
                    start_event_frac = i*nFilesPerJob - math.floor(i*nFilesPerJob)
                    end_event_frac = start_event_frac + nFilesPerJob
                    if (index >= nFiles):
                        continue
                    filesToRun += '%s,' % sample.files[int(index)]
                filesToRun = filesToRun[:-1] # remove trailing ','
                if filesToRun == '':
                    continue

                RunSample_args = 'runOnSkim' if options.runOnSkim else ''
                RunSample_args += ',doSkim' if options.doSkim else ''
                RunSample_args += ',doKinFit' if options.doKinFit else ''

                arguments.append('{0} {1} {2} output_{1}_$(Process).root {3} {4} {5}'.format(configFile, sampleName, filesToRun, str(start_event_frac), str(end_event_frac), RunSample_args))

            submitFile = os.path.join(jobName, sampleName, '%s.submit' % sampleName)
            with open(submitFile, 'w') as f:
                f.write(SUBMIT.render(
                    timestamp=time.strftime('%a %b %d %H:%M:%S %Z %Y'),
                    scriptpath=os.path.abspath(scriptFile),
                    submission=jobName,
                    sample=sampleName,
                    arguments=arguments,
                ))
            submitFiles.append(submitFile)

        if submitJobs:
            print "Job submission files created, submitting jobs to HTCondor..."
            for submitFile in submitFiles:
                os.system('condor_submit %s' % submitFile)
        else:
            print "Job submission files created but not submitted to batch queue"


if __name__ == '__main__':

    status = main()
    sys.exit(status)

