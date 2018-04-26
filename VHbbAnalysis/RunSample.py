import ROOT
import sys
import os
import ReadInput

if (len(sys.argv) != 3 and len(sys.argv) != 5 and len(sys.argv) !=6):
    print "Please give two arguments:  the cfg file and the sample name"
    print "Or give four arguments: the cfg file, the sample name, a comma-separated list of input files, and the name of the output root file (plus comma-separated options: doSkim,runOnSkim,doKinFit)"
    sys.exit(0)

# do stuff :)
ROOT.gSystem.Load("AnalysisDict.so")

# reads samples, existing branches and new branches
samplesToRun = []
samplesToRun.append(sys.argv[2])
filesToRun = []
if len(sys.argv) >= 5:
    for item in sys.argv[3].split(','):
        filesToRun.append(item)

print sys.argv
if len(sys.argv)==6:
    options = sys.argv[5].split(',')
    print "options:", options
    if "doSkim" in options and "runOnSkim" in options:
        raise RuntimeError("Cannot doSkim and runOnSkim at the same time.")
else:
    options = []

am=ReadInput.ReadTextFile(sys.argv[1], "cfg", samplesToRun, filesToRun, 0, "doSkim" in options, "runOnSkim" in options)
#am.debug=20000
am.debug=2

print "Read in the input files, now let's run it!"
if(am.debug>100):
    am.PrintBranches()

print "Done printing branches, now to loop"
# loop over all the samples
# FIXME - need to add the possibility of doing a small portion of files
#aim.Loop()

if (len(sys.argv) == 3):
    am.Loop(sys.argv[2])
else:
    am.Loop(sys.argv[2], ','.join(filesToRun), sys.argv[4], "doSkim" in options)

if "doKinFit" in options:
    import kinfitter

    # mv the previous output_file to a new name, so that the final output name is the same
    output_file = sys.argv[4]
    input_file = output_file.replace('.root', '_before_kinfit.root')
    os.system('mv %s %s' % (output_file, input_file))

    kinfitter.apply_kinfit(
        input_file,
        output_file,
        am.mInt("sampleIndex") == 0,
        am.m("dataYear"),
    )
