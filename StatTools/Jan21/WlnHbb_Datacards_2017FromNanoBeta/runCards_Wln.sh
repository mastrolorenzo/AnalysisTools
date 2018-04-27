python splitSamples.py -i root://cmseos.fnal.gov//store/group/lpchbb/VHbbAnalysisNtuples/2017V2_March27/haddjobs/ -c ttWmn -p "Jet_bReg[hJetInd1]>30&&Jet_bReg[hJetInd2]>30&&MET_Pt<170&&isWmunu&&V_pt>100&&Vtype==2&&controlSample==11" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -tol 0.5    -o $ORIG_DIR/hists_ttWmn_SAMP.root -b $ORIG_DIR/binStats_ttWmn_SAMP.txt  -r 0 -sa SAMP
python splitSamples.py -i root://cmseos.fnal.gov//store/group/lpchbb/VHbbAnalysisNtuples/2017V2_March27/haddjobs/ -c ttWen -p "Jet_bReg[hJetInd1]>30&&Jet_bReg[hJetInd2]>30&&MET_Pt<170&&isWenu&&V_pt>100&&Vtype==3&&controlSample==11" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -tol 0.5    -o $ORIG_DIR/hists_ttWen_SAMP.root -b $ORIG_DIR/binStats_ttWen_SAMP.txt  -r 0 -sa SAMP