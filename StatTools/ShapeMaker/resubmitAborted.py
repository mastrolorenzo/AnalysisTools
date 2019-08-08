import os,sys
from condor_config import *

if "grouping.txt" not in os.listdir("."):
    print "Missing grouping.txt file. Have condor jobs been submitted yet? Exiting."
    sys.exit(1)
    
groupDict = {}
with open("grouping.txt",'r') as gf:
    for line in gf:
        line = line.strip()
        if line == "": continue
        groupDict[line.split(',')[0]] = line

os.system("grep -Ril aborted %s/*.log &> abortedJobs.log"%condorDir)

subf = open("resubmitAborted.sub",'w')
subf.write(subText)
os.system("mkdir -p %s/oldlog"%condorDir)

af = open("abortedJobs.log",'r')
for line in af:
    line = line.strip()
    os.system("mv %s %s/oldlog"%(line,condorDir))
    line = line.lstrip(condorDir+"/log_").rstrip(".log")
    ifile = '_'.join(line.split('_')[:-4])
    start = line.split('_')[-3]
    igroup = groupDict[ifile+".root"]
    logText = logFileText.replace("CONDORDIR",condorDir).replace("ROOTFILE",ifile.rstrip(".root")).replace("CHANSTART",str(start))
    subf.write(logText)
    subf.write("arguments = %s %s %s %s\nqueue\n"%(igroup,nChannelsPerJob,start,indir))
    
af.close()
subf.close()
print "Created resubmitAborted.sub"
os.system("condor_submit resubmitAborted.sub")
