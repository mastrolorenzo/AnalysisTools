python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/2016V4_June19_unblinding/haddjobs/ -c SR_high_Zee -p "hJets_btagged_0>-0.5884&&hJets_btagged_1>-0.5884&&twoResolvedJets&&isZee&&controlSample==0&&V_pt>150&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s systematics_Zllshapes2016.txt -v "1-CMS_vhbb_DNN_Zll_HighPT_13TeV" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_SR_high_Zee_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --even --systematic SYST -bb 0.,0.6376518270819396,0.7142350460131041,0.7633103395929345,0.8001968028007042,0.8300512739472189,0.8552819573537043,0.8772205624033668,0.8966856124007857,0.9142182909998341,0.9301961232350959,0.944893531529949,0.9585166078117847,0.9712242864080837,0.9831418593728792,1.
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/2016V4_June19_unblinding/haddjobs/ -c SR_low_Zee -p "hJets_btagged_0>-0.5884&&hJets_btagged_1>-0.5884&&twoResolvedJets&&isZee&&controlSample==0&&V_pt>50&&V_pt<=150&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s systematics_Zllshapes2016.txt -v "1-CMS_vhbb_DNN_Zll_LowPT_13TeV" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_SR_low_Zee_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --even --systematic SYST -bb 0.,0.6352621184806369,0.7109582520671945,0.759465092475109,0.7959242892218619,0.8254329474841093,0.8503713766631454,0.8720558605712276,0.89129544132322,0.9086250338878344,0.9244177903191237,0.9389449543269904,0.9524102306170428,0.964970712534956,0.9767502408410403,1.0
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/2016V4_June19_unblinding/haddjobs/ -c SR_high_Zuu -p "hJets_btagged_0>-0.5884&&hJets_btagged_1>-0.5884&&twoResolvedJets&&isZmm&&controlSample==0&&V_pt>150&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s systematics_Zllshapes2016.txt -v "1-CMS_vhbb_DNN_Zll_HighPT_13TeV" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_SR_high_Zuu_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --even --systematic SYST -bb 0.,0.6370153562430143,0.7135734294771537,0.7626326094425266,0.7995069611741649,0.8293516297621457,0.8545740288064754,0.8765054304305907,0.8959640891812504,0.9134910110175217,0.9294635970156299,0.9441561794937937,0.9577747827102044,0.9704782888070924,0.9823919486995044,1.
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/2016V4_June19_unblinding/haddjobs/ -c SR_low_Zuu -p "hJets_btagged_0>-0.5884&&hJets_btagged_1>-0.5884&&twoResolvedJets&&isZmm&&controlSample==0&&V_pt>50&&V_pt<=150&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s systematics_Zllshapes2016.txt -v "1-CMS_vhbb_DNN_Zll_HighPT_13TeV" --xlow -1.0 --xhigh 1.0    -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_SR_low_Zuu_SAMP_ISYS.root -r 0 -sa SAMP --year '2016' --channel 'Zll' --even --systematic SYST -bb 0.,0.635149916417177,0.7108309858255353,0.7593281729572101,0.795780114011885,0.8252828997986781,0.8502163660169476,0.8718965345272113,0.8911322864380304,0.9084584302654618,0.924248043803107,0.9387723167810494,0.9522349133659183,0.9647928956404691,0.9765700797196413,1.0