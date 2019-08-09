import os, sys
histfiles = [i for i in os.listdir('.') if i.startswith("hists_") and i.endswith(".root")]

os.system("mkdir -p haddlevel1")
done = []

for fl in histfiles:
    outname = '_'.join(fl.rstrip('.root').split('_')[:-2])
    if outname in done: continue
    if 'data_obs' in fl:
        if '_Zee_' in fl: suff = "Run2016BToGDoubleEleReMiniAOD"
        elif '_Zmm_' in fl: suff = "Run2016BToGDoubleMuReMiniAOD"
        elif '_Wenu_' in fl: suff = "Run2016BToGEleReMiniAOD"
        elif '_Wmunu_' in fl: suff = "Run2016BToGMuReMiniAOD"
        elif '_Znn_' in fl: suff = "Run2016BToGMETReMiniAOD"
        else:
            print "Matching dataset not found."
            sys.exit(1)
        cmd = "hadd haddlevel1/%s.root %s_*%s*.root"%(outname,outname,suff)
    else:
        cmd = "hadd haddlevel1/%s.root %s_*.root"%(outname,outname)
    print "Running command:", cmd
    os.system(cmd)
    done.append(outname)


os.system("mkdir -p haddlevel3")
os.system("hadd haddlevel3/vhcc_Zee-2016.root haddlevel1/*_Zee_*.root")
os.system("hadd haddlevel3/vhcc_Zmm-2016.root haddlevel1/*_Zmm_*.root")
os.system("hadd haddlevel3/vhcc_Wen-2016.root haddlevel1/*_Wenu_*.root")
os.system("hadd haddlevel3/vhcc_Wmn-2016.root haddlevel1/*_Wmunu_*.root")
os.system("hadd haddlevel3/vhcc_Znn-2016.root haddlevel1/*_Znn_*.root")
