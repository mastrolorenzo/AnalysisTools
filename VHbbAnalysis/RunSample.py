from multiprocessing import Process
import sys
import os
import time

if (len(sys.argv) != 3 and len(sys.argv) != 5 and len(sys.argv)!=6 and len(sys.argv)!=7 and len(sys.argv) !=8):
    print "Please give two arguments:  the cfg file and the sample name"
    print "Or give four arguments: the cfg file, the sample name, a comma-separated list of input files, and the name of the output root file (plus comma-separated options: doSkim,runOnSkim,doKinFit)"
    print "Or give six arguments: the cfg file, the sample name, a comma-separated list of input files, the name of the output root file, startFrac, endFrac and comma-separated options: doSkim,runOnSkim,doKinFit"
    sys.exit(61)

# do stuff :)

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
elif len(sys.argv)==8:
    options = sys.argv[7].split(',')
else:
    options = []
print "options:", options
if "doSkim" in options and "runOnSkim" in options:
    raise RuntimeError("Cannot doSkim and runOnSkim at the same time.")

t0=time.time()
def loop_func():
    # this is run in the child process
    import ReadInput
    import ROOT

    ROOT.gSystem.Load("AnalysisDict.so")
    am=ReadInput.ReadTextFile(sys.argv[1], "cfg", samplesToRun, filesToRun, 0, "doSkim" in options, "runOnSkim" in options)
    #am.debug=20000
    am.debug=2

    print "dataYear",am.m("dataYear")
    print "Read in the input files, now let's run it!"
    if(am.debug>100):
        am.PrintBranches()

    print "Done printing branches, now to loop"
    # loop over all the samples
    # FIXME - need to add the possibility of doing a small portion of files
    #aim.Loop()

    if (len(sys.argv) == 3):
        am.Loop(sys.argv[2])
    elif (len(sys.argv) > 6):
        am.Loop(sys.argv[2], ','.join(filesToRun), sys.argv[4],"doSkim" in options, float(sys.argv[5]), float(sys.argv[6]))
    else :
        am.Loop(sys.argv[2], ','.join(filesToRun), sys.argv[4], "doSkim" in options)
    os.system('rm temp.root')
    

p = Process(target=loop_func, args=tuple())
p.start()
p.join()
    
tLoop=time.time()-t0

# load am in the main process only now, in order to keep the memory footprint low
import ReadInput
import ROOT
ROOT.gSystem.Load("AnalysisDict.so")
am=ReadInput.ReadTextFile(sys.argv[1], "cfg", samplesToRun, filesToRun, 0, "doSkim" in options, "runOnSkim" in options)


if "doKinFit" in options:
    import kinfitter

    # mv the previous output_file to a new name, so that the final output name is the same
    output_file = sys.argv[4]
    input_file = output_file.replace('.root', '_before_kinfit.root')
    os.system('mv %s %s' % (output_file, input_file))

    if am.systematics.size():
        # only do lepton systematics if there are other systematics as well
        event_proxies=kinfitter.lep_sys_event_proxies + kinfitter.make_sys_event_proxies(am) + [kinfitter.EventProxy()]
    else:
        event_proxies=None

    kinfitter.apply_kinfit(
        input_file,
        output_file,
        event_proxies=event_proxies,
    )

    # os.system('rm '+input_file)

tKinFit=time.time()-tLoop-t0

if am.branchInfos['postLoopMVAEval'].val > 0.5:

    # get ready for pickling (which happens when calling subprocess..)
    output_file = sys.argv[4]
    input_file = output_file.replace('.root', '_before_mva_eval.root')

    bdt_names = list(name for name, nfo in am.bdtInfos if nfo.mvaType == 'BDT')
    dnn_names = list(name for name, nfo in am.bdtInfos if nfo.mvaType == 'DNN')
    # the first block contains all bdts and the next ones the individual dnn's
    blocks = ([bdt_names] if bdt_names else []) + list([name] for name in dnn_names)


    def worker_func(block):
        import mva_evaluator
        os.system('mv %s %s' % (output_file, input_file))
        mva_evaluator.apply_mva_eval(
            input_file,
            output_file,
            am,
            allowed_names=block
        )
        os.system('rm '+input_file)

    for block in blocks:
        p = Process(target=worker_func, args=(block,))
        p.start()
        p.join()

tMVAEval=time.time()-tKinFit-tLoop-t0

print "tMVAEval,tKinFit,tLoop",tMVAEval,tKinFit,tLoop
