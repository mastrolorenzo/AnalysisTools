python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/2016V4_June19_unblinding/haddjobs/ -c Znn_13TeV_Signal -p "hJets_btagged_0>0.9432&&hJets_btagged_1>-0.5884&&twoResolvedJets&&isZnn&&controlSample==0&&H_mass>60&&H_mass<160" -s "" -v "1-CMS_vhbb_DNN_Znn_13TeV" --xlow -1.0 --xhigh 1.0 -s systematics_Znnshapes2016.txt   -d 1 -n 15 -w 1.0 -o $ORIG_DIR/hists_Znn_13TeV_Signal_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Znn' --even --systematic SYST -bb 0.,0.632328776988211,0.7064500044284682,0.7539476292070444,0.7896482692042642,0.818542981085576,0.8429625504030548,0.8641958748480564,0.8830351639002897,0.9000042032827302,0.9154683813838962,0.9296932985159941,0.942878421178826,0.9551775742912145,0.9667120220277801,1.0
