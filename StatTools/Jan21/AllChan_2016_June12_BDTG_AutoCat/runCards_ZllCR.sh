python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c ttbar_high_Zee -p "isZee&&controlSample==21&&V_pt>150" -s systematics_Zllshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_ttbar_high_Zee_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --drawFromNom True --systematic SYST
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c ttbar_low_Zee -p "isZee&&controlSample==21&&V_pt>50&&V_pt<=150" -s systematics_Zllshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_ttbar_low_Zee_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --drawFromNom True --systematic SYST
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c ttbar_high_Zuu -p "isZmm&&controlSample==21&&V_pt>150" -s systematics_Zllshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_ttbar_high_Zuu_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --drawFromNom True --systematic SYST
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c ttbar_low_Zuu -p "isZmm&&controlSample==21&&V_pt>50&&V_pt<=150" -s systematics_Zllshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_ttbar_low_Zuu_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --drawFromNom True --systematic SYST
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c Zlf_high_Zee -p "isZee&&controlSample==22&&V_pt>150" -s systematics_Zllshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zlf_high_Zee_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --drawFromNom True --systematic SYST
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c Zlf_low_Zee -p "isZee&&controlSample==22&&V_pt>50&&V_pt<=150" -s systematics_Zllshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zlf_low_Zee_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --drawFromNom True --systematic SYST
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c Zlf_high_Zuu -p "isZmm&&controlSample==22&&V_pt>150" -s systematics_Zllshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zlf_high_Zuu_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --drawFromNom True --systematic SYST
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c Zlf_low_Zuu -p "isZmm&&controlSample==22&&V_pt>50&&V_pt<=150" -s systematics_Zllshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zlf_low_Zuu_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --drawFromNom True --systematic SYST
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c Zhf_high_Zee -p "isZee&&controlSample==23&&V_pt>150" -s systematics_Zllshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_high_Zee_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --drawFromNom True --systematic SYST
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c Zhf_low_Zee -p "isZee&&controlSample==23&&V_pt>50&&V_pt<=150" -s systematics_Zllshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_low_Zee_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --drawFromNom True --systematic SYST
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c Zhf_high_Zuu -p "isZmm&&controlSample==23&&V_pt>150" -s systematics_Zllshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_high_Zuu_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --drawFromNom True --systematic SYST
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c Zhf_low_Zuu -p "isZmm&&controlSample==23&&V_pt>50&&V_pt<=150" -s systematics_Zllshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_low_Zuu_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --drawFromNom True --systematic SYST
