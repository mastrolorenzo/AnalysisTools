python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June19_unblinding/haddjobs/ -c whfWmnLow -p "usingBEnriched&&V_pt>=150&&isWmunu&&H_mass<90&&controlSample==13" -s "" -v "TMath::Power(DNN_max_-DNN_2ndmax_,0.33)+DNN_idMax_" --xlow 0.0 --xhigh 5.0 -s systematics_Wlnshapes2017.txt   -d 1 -n 50 -w "(1+(sampleIndex>=5100&&sampleIndex<5103)*(-1+0.25)+(sampleIndex>=5300&&sampleIndex<5303)*(-1+0.9))"  -o $ORIG_DIR/hists_whfWmnLow_SAMP_ISYS.root -r 0 -sa SAMP --year '2017' --systematic SYST --channel 'Wln' --drawFromNom True --friend root://cmseos.fnal.gov://store/group/lpchbb/VHbbNanoPostProc/2017/V5/DNNFRIENDS/dnnCRSR2017-trainedON2017.root  --even
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June19_unblinding/haddjobs/ -c whfWenLow -p "usingBEnriched&&V_pt>=150&&isWenu&&H_mass<90&&controlSample==13" -s "" -v "TMath::Power(DNN_max_-DNN_2ndmax_,0.33)+DNN_idMax_" --xlow 0.0 --xhigh 5.0 -s systematics_Wlnshapes2017.txt   -d 1 -n 50 -w "(1+(sampleIndex>=5100&&sampleIndex<5103)*(-1+0.25)+(sampleIndex>=5300&&sampleIndex<5303)*(-1+0.9))" -o $ORIG_DIR/hists_whfWenLow_SAMP_ISYS.root -r 0 -sa SAMP --year '2017' --systematic SYST --channel 'Wln' --drawFromNom True --friend root://cmseos.fnal.gov://store/group/lpchbb/VHbbNanoPostProc/2017/V5/DNNFRIENDS/dnnCRSR2017-trainedON2017.root --even
