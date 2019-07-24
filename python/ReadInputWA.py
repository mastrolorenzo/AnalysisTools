#~ /usr/bin/python
import sys
from os import listdir
from os.path import isfile, join, isdir
import ROOT
import json
from numpy import array
from ROOT import TH2F
import subprocess
#ROOT.gSystem.Load("SampleContainer_cc.so")
#ROOT.gSystem.Load("AnalysisManager_cc.so")
#ROOT.gSystem.Load("VHbbAnalysis_cc.so")
#ROOT.gSystem.Load("BDTInfo_h.so")
#ROOT.gSystem.Load("VHbbTrigger_h.so")
ROOT.gSystem.Load("AnalysisDict.so")

debug=200

import socket
hostname = socket.gethostname()
site=""
if hostname.endswith('.cern.ch'):
    site="CERN"
    siteIP = "root://eoscms.cern.ch/"
elif hostname.endswith('.desy.de'):
    site="DESY"
    siteIP=""
elif hostname.endswith('.fnal.gov'):
    site="FNAL"
    siteIP = "root://cmseos.fnal.gov/"
else:
    site="FNAL"
    siteIP = "root://cmseos.fnal.gov/"

def ReadTextFile(filename, filetype, samplesToRun="", filesToRun=[], doSkim=False, runOnSkim=False):
    if debug > 10:
         print "filetype is ", filetype
         print "filename is ", filename
         print "samplesToRun is ", samplesToRun
         print "filesToRun is ", filesToRun
         print "doSkim is ", doSkim
         print "runOnSkim is ", runOnSkim

    runSelectedSamples = False
    if (len(samplesToRun) > 0):
        runSelectedSamples = True

    textfile=open(filename, 'r')

    filelines=textfile.readlines()
    filelines=CleanLines(filelines)

    if debug > 1000:
        print filelines
    if filetype is "cfg":

        settings=MakeConfigMap(filelines)

        if settings.has_key("analysis"):
            am=ROOT.__getattr__(settings["analysis"])()
            #debug=100000
        #print "samplesToRun",samplesToRun
        if settings.has_key("samples"):
            aminitialized=0
            samples=ReadTextFile(settings["samples"], "samplefile",samplesToRun,filesToRun,doSkim,runOnSkim)
            for name in samples:
                addedAtLeastOneFile=False
                print site,siteIP
                samples[name]["externFileName"]=samples[name]["globalPrefix"]+name+"_sampleInfo.root"
                filesInGlobalPrefix=GetFileList(samples[name]["globalPrefix"].replace("/store","/eos/uscms/store"), site)
                samples[name]["externFileExists"]=(samples[name]["externFileName"].replace("/store","/eos/uscms/store") in filesInGlobalPrefix)
                print "name externFileExists ",samples[name]["externFileName"],samples[name]["externFileExists"]
                #print "is name in samplesToRun?",name,(name in samplesToRun)
                if (runSelectedSamples and name not in samplesToRun):
                    #print "runSelectedSamples is TRUE and",name,"not in",samplesToRun
                    continue
                sample=samples[name]
                #print sample, sampledic[sample]
                samplecon = ROOT.SampleContainer()
                samplecon.externFileName=siteIP+samples[name]["externFileName"]
                samplecon.externFileExists=samples[name]["externFileExists"]
                if sample.has_key("name"):
                    samplecon.sampleName        = sample["name"]
                if sample.has_key("type"):
                    samplecon.sampleNum         = sample["type"]
                if sample.has_key("xsec"):
                    samplecon.xsec              = sample["xsec"]
                if sample.has_key("kfac"):
                    samplecon.kFactor           = sample["kfac"]
                if sample.has_key("scale"):
                    samplecon.scale             = sample["scale"]
                if sample.has_key("npro"):
                    samplecon.processedEvents   = sample["npro"]
                else:
                    samplecon.nProFromFile      = True
                if sample.has_key("doJetFlavorSplit"):
                    samplecon.doJetFlavorSplit  = bool(sample["doJetFlavorSplit"])
                if sample.has_key("procEff"):
                    samplecon.procEff = sample["procEff"]
                if sample.has_key("lepFlav"):
                    samplecon.lepFlav = sample["lepFlav"]
                if sample.has_key("puhist"):
                    if debug > 10:
                        print "found puhist name in samples"
                    samplecon.PUHistName = sample["puhist"]
                elif settings.has_key("mcpuhistname"):
                    if debug > 10:
                        print "found puhist name in settings"
                    samplecon.PUHistName = settings["mcpuhistname"]
                #print "Reading",sample["name"],"with",len(sample["files"]),"files"

                if samplecon.externFileExists:
                    print "before",samplecon.processedEvents
                    samplecon.ReadSampleInfoFile()
                    print "after",samplecon.processedEvents

                for filename in sample["files"]:
                    #Reset pu hist name for every sample
                    if sample.has_key("puhist"):
                        samplecon.PUHistName = sample["puhist"]
                    elif settings.has_key("mcpuhistname"):
                        samplecon.PUHistName = settings["mcpuhistname"]
                    #print filename,filesToRun
                    if sample["type"]==0 and len(filesToRun)!=0 and len(filesToRun[0])!=0:
                        if filename not in filesToRun:
                            continue
                        else:
                            if debug > 10:
                                print "Files match!",filename

                    # AnalysisManager needs to be initialized
                    # with one file at the beginning
                    if aminitialized == 0:
                        filename.replace("root://cmseos.fnal.gov//store","/eos/uscms/store")
                        print("Initializing with",filename)
                        try:
                            ifile = ROOT.TFile.Open(filename)
                            treeName="Events"#"tree"
                            tree = ifile.Get(treeName)
                            ifile.Close()
                        except:
                            if debug > 10:
                                print "File: %s : no good, trying with another..." % filename
                            continue
                        am.Initialize(filename)
                        if (am.fChain.GetEntries() == 0):
                            if debug > 10:
                                print "File has am.fChain.GetEntries() == 0: %s" % filename
                            continue
                        aminitialized=1
                        # FIXME can this go elsewhere?
                        if settings.has_key("outputname"):
                            am.outputTreeName=settings["outputname"]
                    # if data and fileToRun is not empty then only run that file
                    try:
                        if not runOnSkim:
                            samplecon.AddFile(filename,int(doSkim))
                        else:
                            samplecon.AddFile(filename,2)
                        addedAtLeastOneFile=True
                    except:
                        if debug > 10:
                            print "Can't add",filename
                if addedAtLeastOneFile:
                    print("Adding sample %s to sample container with %i events " % (samplecon.sampleName, samplecon.processedEvents))
                    if not samplecon.externFileExists and samplecon.sampleNum!=0 and not runOnSkim:
                        samplecon.CreateSampleInfoFile()
                        if site=="FNAL":
                            copyToOrigin_command=['xrdcp','-f',samplecon.sampleName+'_sampleInfo.root',samplecon.externFileName]
                            success=subprocess.check_output(copyToOrigin_command)
                        elif site=="DESY":
                            copyToOrigin_command=['gfal-copy','-r',samplecon.sampleName+'_sampleInfo.root','srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN='+samplecon.externFileName]
                            success=subprocess.check_output(copyToOrigin_command)
                        else:
                            copyToOrigin_command=[]
                            print "Need to fill in and execute copy command for CERN"
                    am.AddSample(samplecon)
                else:
                    print("No inputfile could be added for sample %s. Exiting here to avoid seg faults later." % samplecon.sampleName)
                    sys.exit(-1)

        else:
            print "There are no samples in the config file."
            sys.exit(0)
        if settings.has_key("earlybranches"):
            branches=ReadTextFile(settings["earlybranches"], "branchlist",list())
            if debug>10:
                print "getting early branches"
            for branch in branches:
                if debug>10:
                    print(branch,branches[branch][0], branches[branch][1], branches[branch][3], "early", branches[branch][4], branches[branch][5], branches[branch][6])
                am.SetupBranch(branch,branches[branch][0], branches[branch][1], branches[branch][3], "early", branches[branch][4], branches[branch][5],branches[branch][6])
        else:
            print "There are no early branches in the config file."

        if settings.has_key("existingbranches"):
            branches=ReadTextFile(settings["existingbranches"], "branchlist",list())
            for branch in branches:
                am.SetupBranch(branch,branches[branch][0], branches[branch][1], branches[branch][3], "existing", branches[branch][4], branches[branch][5],branches[branch][6])
        else:
            print "There are no existing branches in the config file."

        am.ConfigureOutputTree()
        if debug>10:
            print "output tree configured"

        if settings.has_key("newbranches"):
            branches=ReadTextFile(settings["newbranches"], "branchlist",list())
            for branch in branches:
                if debug>100:
                    print(branch,branches[branch][0], branches[branch][1])
                am.SetupNewBranch(branch,branches[branch][0], branches[branch][1])
        else:
            print "There are no new branches in the config file."

        if settings.has_key("settings"):
            branches=ReadTextFile(settings["settings"], "branchlist",list())
            for branch in branches:
                am.SetupNewBranch(branch,branches[branch][0], branches[branch][1], True, "settings", branches[branch][2])
        else:
            print "There are no settings branches in the config file."

        #find BDT settings
        if not doSkim:
            bdtsettings=[]
            for key in settings:
                if key.find("bdt") != -1:
                    bdtsettings.append(key)


            for bdtsetting in bdtsettings:
                if debug > 10:
                    print "Adding a BDT configuration...",bdtsetting
                bdtInfo=ReadTextFile(settings[bdtsetting], "bdt",list())
                if debug > 10:
                    print "read the BDT settings text file for BDT %s" % bdtInfo.bdtname
                # now set up any of the branches if they don't exist yet (must be floats for BDT)
                for bdtvar in bdtInfo.bdtVars:
                    if (bdtvar.isExisting):
                        am.SetupBranch(bdtvar.localVarName, 2, -1, 0, "early")
                    else:
                        am.SetupNewBranch(bdtvar.localVarName, 2)

                # create new branches for all inputs
                iBDTvar=0
                for bdtvar in bdtInfo.bdtVars:
                    bdtVarPrefix="bdtInput_"+bdtInfo.bdtname+"_"
                    am.SetupNewBranch(bdtVarPrefix+bdtvar.localVarName, 2)
                    #bdtInfo.bdtVars[iBDTvar].localVarName=bdtVarPrefix+bdtvar.localVarName
                    #bdtInfo.bdtVars[iBDTvar].isExisting=False
                    iBDTvar=iBDTvar+1
                am.SetupNewBranch(bdtInfo.bdtname, 2)
                am.SetupNewBranch(bdtsetting, 2, -1, 1, "settings", 1)
                am.AddBDT(bdtsetting, bdtInfo)
                if debug > 10:
                    print "added BDT to analysis manager"

        if settings.has_key("systematics"):
            systs = ReadTextFile(settings["systematics"], "systematics")
            for syst in systs:
                if debug > 10:
                    print "add Systematic"
                am.AddSystematic(syst)
                am.SetupNewBranch("Pass_%s" % syst.name, 4)
                if debug > 10:
                    print "added Systematic"

        if settings.has_key("scalefactors"):
            sfs = ReadTextFile(settings["scalefactors"], "scalefactors")
            for sf in sfs:
                if debug > 10:
                    print "add scale factor"
                am.AddScaleFactor(sf)
                am.SetupNewBranch(sf.branchname, 8, 10)
                am.SetupNewBranch(sf.branchname+"_err", 8, 10)
                if debug > 10:
                    print "added scale factor"

        if settings.has_key("btagscalefactors"):
            if debug > 10:
                print "adding btag scale factors"
            am.InitializeBTagSF(settings["btagscalefactors"])

        if settings.has_key("putarget"):
            try:
                tfile=ROOT.TFile.Open(settings["putarget"])
                puhist=tfile.Get("pileup")
                am.SetGlobalPUTarget(puhist)
                if debug > 10:
                    print "setting putarget"
                tfile.Close()
            except:
                print "something went wrong with",settings["putarget"]

        if settings.has_key("puUPtarget"):
            try:
                tfile=ROOT.TFile.Open(settings["puUPtarget"])
                puhist=tfile.Get("pileup")
                am.SetGlobalPUTarget(puhist,1)
                if debug > 10:
                    print "setting puUPtarget"
                tfile.Close()
            except:
                print "something went wrong with",settings["puUPtarget"]

        if settings.has_key("puDOWNtarget"):
            try:
                tfile=ROOT.TFile.Open(settings["puDOWNtarget"])
                puhist=tfile.Get("pileup")
                am.SetGlobalPUTarget(puhist,-1)
                if debug > 10:
                    print "setting puDOWNtarget"
                tfile.Close()
            except:
                print "something went wrong with",settings["puDOWNtarget"]

        return am
    elif filetype is "samplefile":
        samples=MakeSampleMap(filelines,samplesToRun,runOnSkim,filesToRun)
        #print "writing samples to pickle file"
        #import pickle
        #with open('samples.pickle', 'wb') as fp:
        #    pickle.dump(samples, fp)
        #with open ('samples2.pickle', 'rb') as fp:
        #    samples = pickle.load(fp)
        return samples
    elif filetype is "branchlist":
        branches=MakeBranchMap(filelines)
        return branches
    elif filetype is "bdt":
        bdtInfo=SetupBDT(filelines)
        return bdtInfo
    elif filetype is "systematics":
        systContainers=SetupSyst(filelines)
        return systContainers
    elif filetype is "scalefactors":
        sfContainers=SetupSF(filelines)
        return sfContainers
    else:
        print "Unknown filetype ", filetype


def CleanLines(lines):
    cleanedLines=[]
    for line in lines:
        # skip if it is 100% whitespace
        if line.isspace():
            continue
        # remove comments
        if line.find("#") is -1:        # no comment
            cleanedLines.append(line.split("\n")[0])
        elif line.find("#") is not 0:   # there is an inline comment
            cleanedLines.append(line[:line.find("#")])
        else:                           # whole line is commented
            if debug > 1000:
                print "this line is commented and will be removed"
                print line
    return cleanedLines


def MakeConfigMap(lines):
    if debug > 0: print "In MakeConfigMap"
    settings={}

    for line in lines:
        for item in line.split():
            if len(item.split("=")) == 2:
                name,value=item.split("=")
                settings[name]=value

            else:
                print "not in setting=value format\n\t",item

    return settings


def MakeSampleMap(lines,samplesToRun,runOnSkim=False,filesToRun=[]):
    if debug > 100: print "Reading samples"
    samples={}

    globalPrefix=""
    prefix_lines = dict(
        l.split('=')
        for l in lines
        if l.startswith('prefix')
    )
    if runOnSkim and 'prefix_skim' in prefix_lines:
        globalPrefix = prefix_lines['prefix_skim']
        print "Using globalPrefix for SKIM", globalPrefix
    elif site=="DESY" and 'prefix_desy' in prefix_lines:
        globalPrefix = prefix_lines['prefix_desy']
        print "Using globalPrefix for DESY", globalPrefix
    elif site=="CERN" and 'prefix_cern' in prefix_lines:
        globalPrefix = prefix_lines['prefix_cern']
        print "Using globalPrefix for CERN", globalPrefix
    elif site=="FNAL" and 'prefix_fnal' in prefix_lines:
        globalPrefix = prefix_lines['prefix_fnal']
        print "Using globalPrefix for FNAL", globalPrefix
    elif 'prefix' in prefix_lines:
        globalPrefix = prefix_lines['prefix']
        print "Found a globalPrefix", globalPrefix

    if globalPrefix:
        lastChar=globalPrefix[-1]
        if lastChar.find("/") is -1:
            print "adding a slash to",globalPrefix
            globalPrefix=globalPrefix+"/"

    for line in lines:
        #print line
        samplename=""
        samplepaths=[]
        sampletype=-1
        samplexsec=-1
        samplekfac=1
        samplescale=1
        samplepuhist=""
        sample={}
        sample["globalPrefix"]=globalPrefix

        dontRun = False
        for item in line.split():
            name,value=item.split("=")
            if name.find("name") is 0:
                if len(samplesToRun)>0 and str(value) not in samplesToRun:
                    line = ""
                    dontRun = True
        if dontRun: continue

        for item in line.split():
            name,value=item.split("=")
            if name.find("name") is 0:
                sample["name"]=str(value)
            if name.find("file") is 0:
                samplepaths.append(globalPrefix+str(value))
            if name.find("puhist") is 0:
                sample["puhist"]=str(value)
            if name.find("dir") is 0:
                if len(filesToRun) > 0:
                    # no need to find all the files again just run on what is given
                    samplepaths = filesToRun
                    continue
                if site == "CERN":
                    value = value.replace("CERN:","")
                if value.find(',') is not 0:
                    if not runOnSkim:
                        for dirname in value.split(','):
                            globalPrefix = globalPrefix.replace("/store","/eos/uscms/store/")
                            print dirname
                            print globalPrefix,dirname,site
                            samplepaths.extend(findAllRootFiles(globalPrefix+dirname,site))
                    else:
                        # after skimming the sample directory name has been changed to sample name
                        if debug > 10:
                            print "right here, ",sample["name"]
                        samplepaths.extend(findAllRootFiles(globalPrefix+sample["name"],site))
                else:
                    samplepaths = findAllRootFiles(globalPrefix+value,site)
                #print value
                #if value.find("/store") is 0:
                #    import subprocess
                #    onlyFiles = subprocess.check_output(["/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_14/external/slc6_amd64_gcc491/bin/xrdfs", "root://cmseos.fnal.gov", "ls", value]).split('\n')
                #    for filepath in onlyFiles:
                #        #filepath = "root://xrootd-cms.infn.it/" + filepath
                #        #filepath = "root://cmsxrootd.fnal.gov/" + filepath
                #        filepath = "root://cmseos.fnal.gov/" + filepath
                #        if filepath.find(".root") != -1:
                #            samplepaths.append(filepath)
                #else:
                #    from os import listdir
                #    from os.path import isfile, join
                #    onlyfiles = [ f for f in listdir(str(value)) if isfile(join(str(value),f)) ]
                #    for rootfile in onlyfiles:
                #        #print rootfile
                #        if rootfile.find(".root") != -1:
                #            samplepaths.append(str(value)+"/"+str(rootfile))
            if name.find("type") is 0:
                sample["type"]=int(value)
            if name.find("xsec") is 0:
                sample["xsec"]=float(value)
            if name.find("kfac") is 0:
                sample["kfac"]=float(value)
            if name.find("scale") is 0:
                sample["scale"]=float(value)
            if name.find("npro") is 0:
                sample["npro"]=int(value)
            if name.find("doJetFlavorSplit") is 0:
                sample["doJetFlavorSplit"]=bool(value)
            if name.find("procEff") is 0:
                sample["procEff"]=float(value)
            if name.find("lepFlav") is 0:
                sample["lepFlav"]=int(value)

        sample["files"]=samplepaths
        if sample.has_key("name"):
            samples[sample["name"]]=sample
        else:
            if debug > 10:
                print "sample name is empty",samplename,"not filling"

    return samples


def MakeBranchMap(lines):
    branches={}

    for line in lines:
        branchname=""
        branchtype=-1
        arraylength=-1
        val=-999
        onlyMC=0
        lengthBranch=""
        allowMissingBranch=0
        dropBranch=0

        for item in line.split():
            name,value = item.split("=")
            if name.find("name") is 0:
                branchname=value
            if name.find("type") is 0:
                branchtype=int(value)
            if name.find("max") is 0:
                arraylength=int(value)
            if name.find("val") is 0:
                val=float(value)
            if name.find("onlyMC") is 0:
                onlyMC=int(value)
            if name.find("lengthBranch") is 0:
                if debug>1000:
                    print "FOUND LENGTH BRANCH",value
                lengthBranch=str(value)
            if name.find("allowMissingBranch") is 0:
                allowMissingBranch=int(value)
            if name.find("dropBranch") is 0:
                dropBranch=int(value)


        branches[branchname]= [branchtype,arraylength,val,onlyMC,lengthBranch,allowMissingBranch,dropBranch]

    return branches

def SetupBDT(lines):
    bdtname = ""
    methodName = ""
    xmlFile = ""
    inputNames = []
    localVarNames = []
    mvaType = "BDT"

    vars = {}

    for line in lines:
        inputName = ""
        localVarName = ""
        isSpec = False
        order = -1
        isExisting = False
        for item in line.split():
            name,value = item.split("=")
            if name.find("bdtname") is 0:
                bdtname=value
                break
            if name.find("methodName") is 0:
                methodName=value
                break
            #if name.find("method") is 0:
            #    method=value.replace('@', ' ')
            #    print "method set to %s" % method
            #    break
            if name.find("mvaType") is 0:
                mvaType=value
                break
            if name.find("xmlFile") is 0:
                xmlFile=value
                break
            if name.find("name") is 0:
                inputName=value
            if name.find("lname") is 0:
                localVarName=value
            if name.find("isSpec") is 0:
                if (int(value) == 1): isSpec = True
            if name.find("order") is 0:
                order=int(value)
            if name.find("isEx") is 0:
                if (int(value) == 1): isExisting = True
        vars[order] = (inputName,localVarName,isExisting,isSpec)

    bdt = ROOT.BDTInfo(methodName, bdtname, xmlFile, mvaType)

    keys = vars.keys()
    keys.sort()

    for key in keys:
        if (key == -1): continue
        name, lname, isExisting, isSpec = vars[key]
        if debug > 10:
            if not isSpec:
                print "adding variable %s (%s) existing: %i " % (name,lname,int(isExisting))
            else:
                print "adding spectator variable %s (%s) existing: %i" % (name,lname, int(isExisting))
        bdt.AddVariable(name, lname, isExisting, isSpec)
    return bdt


def SetupSyst(lines):
    systs=[]
    for line in lines:
        syst=ROOT.SystematicContainer()
        for item in line.split():
            key,value=item.split("=")
            if key=="name":
                syst.name=value
            elif key=="branches":
                for brnchName in value.split(","):
                    syst.AddBranchName(brnchName)
            elif key=="scales":
                scales=[]
                for scale in value.split(","):
                    syst.AddScale(float(scale))
            elif key=="smears":
                smears=[]
                for smear in value.split(","):
                    syst.AddSmear(float(smear))
            elif key=="scaleVar":
                for scalevar in value.split(","):
                    syst.AddScaleVar(scalevar)
            elif key=="scaleVarRef":
                for scalevarref in value.split(","):
                    syst.AddScaleVarRef(scalevarref)
            else:
                print "In systematics file, what is:",item

            # FIXME need to add passing 2d histogram from filename

        systs.append(syst)
    return systs


def SetupSF(lines):
    SFs=[]
    for line in lines:
        print line
        SF=ROOT.SFContainer()
        for item in line.split():
            key,value=item.split("=")
            if key=="json":
                SF.jsonFile=value
            elif key=="branches":
                for brnchName in value.split(","):
                    SF.AddBranch(brnchName)
            elif key=="name":
                SF.name=value
            elif key=="binning":
                SF.binning=value
            elif key=="branchname":
                SF.branchname=value
            elif key=="length":
                SF.length=value
            else:
                print "In scale factor file, what is:",item
        # now parse the json and build the scale factor map
        f = open(SF.jsonFile, 'r')
        results = json.load(f)
        if SF.name not in results.keys():
            print SF.name,": not found in json file: ",SF.jsonFile
        f.close()
        stripForEta = 5
        res = results[SF.name]
        if SF.binning not in res.keys():
            print SF.binning," not in list of binnings: ",res.keys()
        if "abseta" in SF.binning:
            stripForEta = 8

        etaBins = []
        ptBins = []
        nIter = 0

        if (SF.binning.find("pt") != 0 and SF.binning.find("DATA") == -1):
            # categorized first by eta
            for etaKey, values in sorted(res[SF.binning].iteritems()):
                etaL = float(((etaKey[stripForEta:]).rstrip(']').split(',')[0]))
                etaH = float(((etaKey[stripForEta:]).rstrip(']').split(',')[1]))
                etaBins.append(etaL)
                etaBins.append(etaH)

                for ptKey, result in sorted(values.iteritems()):
                    ptL = float(((ptKey[4:]).rstrip(']').split(',')[0]))
                    ptH = float(((ptKey[4:]).rstrip(']').split(',')[1]))
                    ptBins.append(ptL)
                    ptBins.append(ptH)
        else:
            # categorized first by pt
            for ptKey, values in sorted(res[SF.binning].iteritems()):
                if debug > 10:
                    print ptKey
                    print (ptKey[stripForEta:])
                ptL = float(((ptKey[4:]).rstrip(']').split(',')[0]))
                ptH = float(((ptKey[4:]).rstrip(']').split(',')[1]))
                ptBins.append(ptL)
                ptBins.append(ptH)

                for etaKey, result in sorted(values.iteritems()):
                    #print etaKey
                    #print (etaKey[4:]).rstrip(']').split(',')
                    etaL = float(((etaKey[stripForEta:]).rstrip(']').split(',')[0]))
                    etaH = float(((etaKey[stripForEta:]).rstrip(']').split(',')[1]))
                    etaBins.append(etaL)
                    etaBins.append(etaH)

        etaBins = list(set(etaBins)) # get rid of duplicates
        ptBins = list(set(ptBins))
        etaBins = sorted(etaBins)
        ptBins = sorted(ptBins)

        #print etaBins, ptBins
        SF.scaleMap = TH2F(SF.name, SF.name, len(ptBins)-1, array(ptBins), len(etaBins)-1, array(etaBins))
        #print array(ptBins), array(etaBins)

        if (SF.binning.find("pt") != 0 and SF.binning.find("DATA") ==- 1):
            # categorized first by eta
            for etaKey, values in sorted(res[SF.binning].iteritems()):
                etaL = float(((etaKey[stripForEta:]).rstrip(']').split(',')[0]))
                etaH = float(((etaKey[stripForEta:]).rstrip(']').split(',')[1]))
                binLow1 = SF.scaleMap.GetYaxis().FindBin(etaL)
                binHigh1 = SF.scaleMap.GetYaxis().FindBin(etaH) - 1
                for ptKey, result in sorted(values.iteritems()):
                    ptL = float(((ptKey[4:]).rstrip(']').split(',')[0]))
                    ptH = float(((ptKey[4:]).rstrip(']').split(',')[1]))
                    binLow2 = SF.scaleMap.GetXaxis().FindBin(ptL)
                    binHigh2 = SF.scaleMap.GetXaxis().FindBin(ptH) - 1
                    for ibin1 in range(binLow1,binHigh1+1):
                        for ibin2 in range(binLow2,binHigh2+1):
                            #SF.scaleMap.Fill( (ptL+ptH)/2., (etaL+etaH)/2., result["value"] )
                            #SF.scaleMap.SetBinError( SF.scaleMap.GetXaxis().FindBin( (ptL+ptH)/2.), SF.scaleMap.GetYaxis().FindBin( (etaL+etaH)/2.), result["error"] )
                            SF.scaleMap.SetBinContent(ibin2,ibin1, result["value"])
                            SF.scaleMap.SetBinError(ibin2,ibin1, result["error"])
                            #print etaL,etaH,ptL,ptH,ibin2,ibin1,SF.scaleMap.GetXaxis().GetBinLowEdge(ibin1),SF.scaleMap.GetYaxis().GetBinLowEdge(ibin2), result["value"], result["error"]
        else:
            # categorized first by pt
            for ptKey, values in sorted(res[SF.binning].iteritems()):
                ptL = float(((ptKey[4:]).rstrip(']').split(',')[0]))
                ptH = float(((ptKey[4:]).rstrip(']').split(',')[1]))
                binLow1 = SF.scaleMap.GetXaxis().FindBin(ptL)
                binHigh1 = SF.scaleMap.GetXaxis().FindBin(ptH) - 1
                for etaKey, result in sorted(values.iteritems()):
                    etaL = float(((etaKey[stripForEta:]).rstrip(']').split(',')[0]))
                    etaH = float(((etaKey[stripForEta:]).rstrip(']').split(',')[1]))
                    binLow2 = SF.scaleMap.GetYaxis().FindBin(etaL)
                    binHigh2 = SF.scaleMap.GetYaxis().FindBin(etaH) - 1
                    for ibin1 in range(binLow1,binHigh1+1):
                        for ibin2 in range(binLow2,binHigh2+1):
                            #SF.scaleMap.Fill( (ptL+ptH)/2., (etaL+etaH)/2., result["value"] )
                            #SF.scaleMap.SetBinError( SF.scaleMap.GetXaxis().FindBin( (ptL+ptH)/2.), SF.scaleMap.GetYaxis().FindBin( (etaL+etaH)/2.), result["error"] )
                            SF.scaleMap.SetBinContent(ibin1,ibin2, result["value"])
                            SF.scaleMap.SetBinError(ibin1,ibin2, result["error"])
                            #print etaL,etaH,ptL,ptH,ibin2,ibin1,SF.scaleMap.GetXaxis().GetBinLowEdge(ibin1),SF.scaleMap.GetYaxis().GetBinLowEdge(ibin2), result["value"], result["error"]
        SFs.append(SF)
    for SF in SFs:
        if debug > 10:
            print "debugging scalefactor ",SF.name
        sfmap = SF.scaleMap
        nX = sfmap.GetNbinsX()
        nY = sfmap.GetNbinsY()
        if debug > 1000:
            for i in range(1,nX+1):
                print "pt: ",sfmap.GetXaxis().GetBinLowEdge(i)
                for j in range(1, nY+1):
                    print "eta: ",sfmap.GetYaxis().GetBinLowEdge(j),": ",sfmap.GetBinContent(i,j)
    return SFs

def GetFileList(value, site):
    if value.find("/store") is 0:
        onlyFiles = subprocess.check_output(["xrdfs", siteIP, "ls", value]).split('\n')
    else:
        #onlyfiles = [ f for f in listdir(str(value)) if isfile(join(str(value),f)) ]
        onlyFiles = list(join(value,f) for f in listdir(str(value)))
    return onlyFiles


def findAllRootFiles(value, site):
    samplepaths = []
    if value.find("/store") is 0:
        #onlyFiles = subprocess.check_output(["/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_14/external/slc6_amd64_gcc491/bin/xrdfs", siteIP, "ls", value]).split('\n')
        onlyFiles = subprocess.check_output(["xrdfs", siteIP, "ls", value]).split('\n')
        for filepath in onlyFiles:
            if (filepath == ""): continue
            #filepath = "root://xrootd-cms.infn.it/" + filepath
            #filepath = "root://cmsxrootd.fnal.gov/" + filepath
            if filepath.find(".root") != -1:
                filepath = siteIP + filepath
                samplepaths.append(filepath)
            elif filepath.find("/log/")==-1:
                samplepaths.extend(findAllRootFiles(filepath,site))
    else:
        #onlyfiles = [ f for f in listdir(str(value)) if isfile(join(str(value),f)) ]
        onlyFiles = listdir(str(value))
        for rootfile in onlyFiles:
            if rootfile.find(".root") != -1:
                tmp = str(value)+"/"+str(rootfile)
                tmp = tmp.replace("/eos/uscms/store","root://cmseos.fnal.gov//store")
                print tmp
                samplepaths.append(tmp)
                #samplepaths.append(str(value)+"/"+str(rootfile))
            elif (isdir(join(str(value),rootfile))):
                print "recursive, ",str(value),"/",str(rootfile)
                samplepaths.extend(findAllRootFiles(str(value)+"/"+str(rootfile),site))
    return samplepaths