python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c ttWmn -p "V_pt>=150&&isWmunu&&controlSample==11" -s systematics_Wlnshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_ttWmn_SAMP.root -r 0 -sa SAMP --year '2016' --channel 'Wln' --drawFromNom True
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c ttWen -p "V_pt>=150&&isWenu&&controlSample==11" -s systematics_Wlnshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0 -o $ORIG_DIR/hists_ttWen_SAMP.root -r 0 -sa SAMP --year '2016' --channel 'Wln' --drawFromNom True
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c whfWmnLow -p "V_pt>=150&&isWmunu&&H_mass<90&&controlSample==13" -s systematics_Wlnshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_whfWmnLow_SAMP.root -r 0 -sa SAMP --year '2016' --channel 'Wln' --drawFromNom True
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c whfWenLow -p "V_pt>=150&&isWenu&&H_mass<90&&controlSample==13" -s systematics_Wlnshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0 -o $ORIG_DIR/hists_whfWenLow_SAMP.root -r 0 -sa SAMP --year '2016' --channel 'Wln' --drawFromNom True
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c whfWmnHigh -p "V_pt>=150&&isWmunu&&H_mass>150&&controlSample==13" -s systematics_Wlnshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_whfWmnHigh_SAMP.root -r 0 -sa SAMP --year '2016' --channel 'Wln' --drawFromNom True
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c whfWenHigh -p "V_pt>=150&&isWenu&&H_mass>150&&controlSample==13" -s systematics_Wlnshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0 -o $ORIG_DIR/hists_whfWenHigh_SAMP.root -r 0 -sa SAMP --year '2016' --channel 'Wln' --drawFromNom True
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c wlfWmn -p "V_pt>=150&&isWmunu&&controlSample==12" -s systematics_Wlnshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_wlfWmn_SAMP.root -r 0 -sa SAMP --year '2016' --channel 'Wln' --drawFromNom True
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/CMSConnect_June9_2016V4/haddjobs/ -c wlfWen -p "V_pt>=150&&isWenu&&controlSample==12" -s systematics_Wlnshapes2016.txt -v "Jet_btagCMVA[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0 -o $ORIG_DIR/hists_wlfWen_SAMP.root -r 0 -sa SAMP --year '2016' --channel 'Wln' --drawFromNom True
