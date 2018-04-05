import time
import sys,os
import getopt
import subprocess

usage = """
tier2transfer.py [options]
--input-dir, -i                    origin directory
--output-dir, -o                   destination directory
--wildcard, -w                     selects only directories with w to transfer
--source, -s                       pisa,fnal,cern,psi
--destination, -d                  cern,fnal
--nolsmode, -n                     don't check if files already exist, default off
--checkmode, -c                    check if all files are transferred, default off
--copyMissingDirsOnly              if a directory does not exist in destination, copy it, default off
--copyMissingFilesOnly             if a file does not exist in destination, copy  it, default off
--mode                             either LCG or XRD, default XRD
"""

pisapre='srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/'
#cernpre='srm://srm-eoscms.cern.ch:8443/srm/v2/server?SFN=/eos/cms/store/'
cernpre='gsiftp://eoscmsftp.cern.ch//eos/cms/store/'
fnalpre='srm://cmseos.fnal.gov:8443/srm/v2/server?SFN=/eos/uscms/store/'
psipre='srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/'
desypre="srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN=/pnfs/desy.de/cms/tier2/store/"

sourcepre=''
destpre=''

source=''
dest=''

inputDir=''
outputDir=''

wildcard=''

copyMissingDirsOnly=False
copyMissingFilesOnly=False
checkmode=False
nolsmode=False
recursive=True
mode="GFAL"

(opts, args) = getopt.getopt(sys.argv[1:], 'd:i:o:w:s:hcnr:', ['copyMissingDirsOnly=','copyMissingFilesOnly=','source=','destination=','input-dir=', 'output-dir=', 'help', 'wildcard=','nolsmode','checkmode','recursive','mode='])

for opt,argm in opts:
    #print opt,argm
    if (opt == "--help" or opt == "-h"):
        print 'Usage: %s' % (usage)
        sys.exit(0)
    elif (opt == "-s" or opt == "--source"):
        source = argm
    elif (opt == "-d" or opt == "--destination"):
        dest = argm
    elif (opt == "-i" or opt == "--input-dir"):
        inputDir = argm
    elif (opt == "-o" or opt == "--output-dir"):
        outputDir = argm
    elif (opt == "-w" or opt == "--wildcard"):
        wildcard = argm
    elif (opt == "-c" or opt == "--checkmode"):
        checkmode=True
    elif (opt == "-n" or opt == "--nolsmode"):
        nolsmode=True
    elif (opt == "-r" or opt == "--recursive"):
        recursive=argm
    elif (opt == "--mode"):
        mode = argm
    elif (opt == "--copyMissingDirsOnly"):
        copyMissingDirsOnly=argm
    elif (opt == "--copyMissingFilesOnly"):
        copyMissingFilesOnly=argm
    else:
        print 'Wrong options: %s' % (opt)
        sys.exit(3)


if source.find('pisa') != -1:
    sourcepre=pisapre
elif source.find('cern') != -1:
    sourcepre=cernpre
elif source.find('fnal') != -1:
    sourcepre=fnalpre
elif source.find('psi') != -1:
    sourcepre=psipre
elif source.find('desy') != -1:
    sourcepre=desypre
else:
    print "source",source,"is not setup"
    sys.exit(0)

if dest.find('local') != -1:
    destpre=''
elif dest.find('cern') != -1:
    destpre=cernpre
elif dest.find('fnal') != -1:
    destpre=fnalpre
elif dest.find('psi') != -1:
    destpre=psipre
elif dest.find('desy') != -1:
    destpre=desypre
else:
    print "destination",dest,"is not setup"
    sys.exit(0)

if mode=="XRD":
    sourceIPAdd=""
    destIPAdd=""
    if source.find('cern') != -1:
        sourceIPAdd="188.184.38.46"
    else:
        print "XRD can't do source",source,"exiting..."
        sys.exit(-1)
else:
    print "mode",mode,"doesn't exist"



def LS_WRAP(args):
    print "source in LS_WRAP",source
    if mode=="LCG":
        cmd=['lcg-ls','-b','-D','srmv2']
    elif mode=="GFAL":
        cmd=['gfal-ls']
    else:
        print "what is",mode,"? exiting"
        sys.exit(0)
    cmd.extend(args)
    print cmd
    checkOutput=""
    try:
        checkOutput = subprocess.check_output(cmd)
    except:
        print "problems with",cmd
    return checkOutput


def XRDFS_WRAP(ipAdd,path):
    cmd=["xrdfs",ipAdd,"ls","-u",path]
    print cmd
    checkOutput=""
    try:
        checkOutput = subprocess.check_output(cmd)
    except:
        print "problems with",cmd
    lsOut=[]
    for path in checkOutput.split(":1094/"):
        if path.find("/store")!=-1:
            lsOut.append(path.split("\n")[0])
    print "len(lsOut)",len(lsOut)
    return lsOut


def CP_WRAP(input, output):
    print "input",input
    print "output",output
    if mode=="LCG":
        cmd=['lcg-cp','-b','-D','srmv2',input,output]
    elif mode=="GFAL":
        #cmd=['gfal-copy',input,output.split("/tree_")[0]]
        cmd=['gfal-copy','-p',input,output]
        #if not os.path.isdir(output.split("r?SFN=")[1].split("/tree_")[0]):
        #    os.makedirs(output.split("r?SFN=")[1].split("/tree_")[0])
    elif mode=="XRD":
        cmd=["xrdcp","-r",input, output]
    print cmd
    subprocess.Popen(cmd)


def DoesFileExist(inpath,outpath):
    try:
        #print inpath
        outls=LS_WRAP(["-l",outpath])
        #check that file is the right size
        inls=LS_WRAP(["-l",inpath])
        #print "insize, outsize,",inls.split()[4],",",outls.split()[4]
        if inls.split()[4] == outls.split()[4]:
            #same size
            return True
        else:
            print inpath
            print "size different",inls.split()[4],outls.split()[4]
            return False
    except Exception as e:
        print "DoesFileExist failed"
        #print "ERROR",e
        return False


ignoreList=["/log","/failed",".tmp"]
def CheckIgnoreList(filePath):
    ignore=False
    if filePath is "":
        ignore=True
    else:
        for toIgnore in ignoreList:
            if filePath.find(toIgnore) !=-1:
                ignore=True
                break
    return ignore


def MakeFileList(rootpath):
    print "rootpath",rootpath
    print "mode",mode
    files=[]
    if mode=="LCG" or mode=="GFAL":
        output=LS_WRAP([rootpath])  
        list=output.split("\n")
        print "list of contents",list
        for item in list:
            if len(item)>0:
                print "raw item",item
                #if mode=="GFAL":
                #    splitPhrase=sourcepre.split("SFN=")[1]
                #    #print "splitPhrase"
                #    #print "new item",item.split(splitPhrase)
                #    item=item.split(splitPhrase)[-1]
                #    print "GFAL item",item
                item=rootpath+"/"+item
                print "final item",item
                # removing things here could be problematic
                # if there are many subdirs without datasetnames
                if wildcard is not "":
                    if item.find(wildcard) == -1:
                        print "skipping",item
                        continue
                if item.find(".root") != -1:
                    print "found root file",item
                    files.append(item)
                    #print "Adding",item,"to list"
                elif recursive == True:
                    ignore=CheckIgnoreList(item)
                    print "ignore?",ignore
                    if not ignore:
                        try:
                            #files.extend(MakeFileList(sourcepre+item))
                            print "found dir",item
                            files.extend(MakeFileList(item))
                        except:
                            print "The following item fails "+item
        return files
    elif mode=="XRD":
        lsOut=XRDFS_WRAP(sourceIPAdd,rootpath)
        #print "in makefilelist.  loop over",len(lsOut)
        for item in lsOut:
            #print "in makefilelist.  item",item
            if wildcard is not "":
                if item.find(wildcard) == -1:
                    print "skipping",item
                    continue
            if item.find(".root") != -1:
                files.append(item)
            elif recursive == True:
                ignore=CheckIgnoreList(item)
                if not ignore:
                    try:
                        files.extend(MakeFileList(item))
                    except:
                        print "Can't LS this path",lsOut

        return files


def SetMaxThreads():
    cmd=['cat','/proc/cpuinfo']
    output=subprocess.check_output(cmd)
    corelines=[line for line in output.split('\n') if line.find("cpu cores") != -1]
    cores=0
    for coreline in corelines:
        cores=cores+int(coreline.split()[-1])

    return cores


if mode=="LCG":
    cmdCheck="lcg-cp"
elif mode=="GFAL":
    cmdCheck="gfal-copy"
elif mode=="XRD":
    cmdCheck="xrdcp"
        

MaxThreads=max(SetMaxThreads()/4,3)

listOfExistingInputDirs=[]
if mode=="LCG" or mode=="GFAL":
    listOfExistingInputDirs=LS_WRAP([sourcepre+inputDir])
elif mode=="XRD":
    listOfExistingInputDirs=XRDFS_WRAP(sourceIPAdd,inputDir)
listOfExistingInputDirs=listOfExistingInputDirs.split("\n")

listOfExistingOutputDirs=[]
if copyMissingDirsOnly:
    if mode=="LCG" or mode=="GFAL":
        listOfExistingOutputDirs=LS_WRAP([destpre+outputDir])
    elif mode=="XRD":
        listOfExistingOutputDirs=XRDFS_WRAP(destIPAdd,outputDir)
    listOfExistingOutputDirs=listOfExistingOutputDirs.split("\n")
    

print "n directories",len(listOfExistingInputDirs)
for directory in listOfExistingInputDirs:
    if directory is "":
        print "blank string... skipping"
        continue

    if copyMissingDirsOnly and directory in listOfExistingOutputDirs:
        print "This directory is already present and only new directories are being copied.",directory
        continue

    if directory.find(wildcard) == -1:
        print "skipping",directory
        continue
    
    thisInputDir=inputDir+"/"+directory
    thisOutputDir=outputDir+"/"+directory
    
    if mode=="LCG" or mode=="GFAL":
        print "directory",directory
        print "sourcepre",sourcepre
        print "thisInputDir",thisInputDir
        print "thisOutputDir",thisOutputDir
        
        files=MakeFileList(sourcepre+thisInputDir)
    elif mode=="XRD":
        files=MakeFileList(thisInputDir)
    print "n files",len(files)
        
    if copyMissingFilesOnly:
        filesAlreadyThere=[]
        if mode=="LCG" or mode=="GFAL":
            filesAlreadyThere=MakeFileList(destpre+thisOutputDir)
        print "filesAlreadyThere",filesAlreadyThere
        for iOutFile in range(len(filesAlreadyThere)):
            filesAlreadyThere[iOutFile]=filesAlreadyThere[iOutFile].split(destpre+thisOutputDir)[1]
        filesToGo=[] 
        for iInFile in range(len(filesToGo)):
            thisFile=files[iInFile].split(sourcepre+thisInputDir)[1]
            if thisFile not in filesAlreadyThere:
                filesToGo.append(files)
        print "len(filesToGo),len(files)",len(filesToGo),len(files)
        files=filesToGo

    ifile=0
    isleep=0
    newsleep=1
    untransfiles=[]
    
    if len(files) is 0:
        print "There are no files to transfer."
        continue
    
    print "MaxThreads files ",MaxThreads,len(files)
    print "first file name",files[0]
    while ifile!=len(files):
        out=[line for line in subprocess.check_output(["ps"]).split("\n") if line.find(cmdCheck) !=-1 ]
       
        #if ifile>13:
        #    break
    
        if len(out)< MaxThreads:
            if mode=="LCG" or mode=="GFAL":
                outpath=files[ifile].replace(sourcepre+inputDir,destpre+outputDir)
                inpath=files[ifile]
    
                transfer=True
                if nolsmode is False:
                    transfer= not DoesFileExist(inpath, outpath)
                
                if transfer is True:
                    if checkmode is True:
                        untransfiles.append(outpath) 
                    else:
                        CP_WRAP(inpath, outpath)
                        newsleep=1
            elif mode=="XRD":
                outpath=files[ifile].replace(thisInputDir,thisOutputDir)
                outdir=outpath.split("/tree_")[0]
                if not os.path.isdir(outdir):
                    os.makedirs(outdir)
                CP_WRAP("root://"+sourceIPAdd+"/"+files[ifile],outdir)
            
            ifile=ifile+1
        else:
            isleep=isleep+1
            if newsleep == 1:
                print "sleeping",isleep
                print "ifile",ifile,"of",len(files)
                newsleep=0
            time.sleep(2)
    
    if checkmode is True:
        print "Missing:"
        for file in untransfiles:
            print file
        print len(untransfiles),"of",len(files),"to go"
          
    else:
        out=[line for line in subprocess.check_output(["ps"]).split("\n") if line.find(cmdCheck) !=-1 ]
        print "Number of on-going "+cmdCheck+"s:",len(out)
