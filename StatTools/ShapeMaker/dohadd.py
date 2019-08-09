import os, sys
histfiles = [i for i in os.listdir('.') if i.startswith("hists_") and i.endswith(".root")]

os.system("mkdir -p haddlevel1")
done = []

for fl in histfiles:
    outname = '_'.join(fl.rstrip('.root').split('_')[:-2])
    #hack escape characters where needed
    safe_outname = []
    for char in outname:
        if char == "=" or char == "&" or char == ")" or char == "(":
            safe_outname.append("\\" + char)
        else:
            safe_outname.append(char)
    outname = "".join(safe_outname)

    positions = [pos for pos, char in enumerate(outname) if char == "="]

    if outname in done: continue
    if 'data_obs' in fl:
        if '_Znn_' in fl: 
            suff = "Run2017BToGMETReMiniAOD"
            cmd = 'hadd haddlevel1/%s.root %s_*METMiniAOD*.root %s_*EleReMiniAOD*.root %s_*MuReMiniAOD*.root'%(outname,outname,outname,outname)
        else:
            if '_Zee_' in fl: suff = "DoubleEleReMiniAOD"
            elif '_Zmm_' in fl: suff = "DoubleMuReMiniAOD"
            elif '_Wen_' in fl: suff = "EleReMiniAOD"
            elif '_Wmn_' in fl: suff = "MuReMiniAOD"
            else:
                print "%s: Matching dataset not found."%(fl)
                sys.exit(1)
            cmd = "hadd haddlevel1/%s.root %s_*%s*.root"%(outname,outname,suff)
    else:
        cmd = "hadd haddlevel1/%s.root %s_*.root"%(outname,outname)
    print "Running command:", cmd
    os.system(cmd)
    done.append(outname)

os.system("mkdir -p haddlevel3")
os.system("hadd haddlevel3/vhbb_Zee.root haddlevel1/*_Zee_*.root")
os.system("hadd haddlevel3/vhbb_Zmm.root haddlevel1/*_Zmm_*.root")
os.system("hadd haddlevel3/vhbb_Wen.root haddlevel1/*_Wen_*.root")
os.system("hadd haddlevel3/vhbb_Wmn.root haddlevel1/*_Wmn_*.root")
os.system("hadd haddlevel3/vhbb_Znn.root haddlevel1/*_Znn_*.root")
