python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June19_unblinding/haddjobs/ -c SR_high_Zee -p "usingBEnriched&&hJets_btagged_0>0.1522&&hJets_btagged_1>0.1522&&twoResolvedJets&&isZee&&controlSample==0&&V_pt>150&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s "" -v "1-CMS_vhbb_DNN_Zll_HighPT_13TeV" --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -bb 0.,0.6388302270218088,0.7163691683600514,0.7660568984287467,0.8034036872783153,0.833630728126449,0.8591762785470664,0.8813886670275982,0.9010966316569309,0.9188481098053722,0.9350253378416581,0.9499061628514183,0.9636992486919744,0.9765655131161601,0.9886318117636758,1.0 -w 1.0  -o $ORIG_DIR/hists_SR_high_Zee_SAMP_ISYS.root -r 0 -sa SAMP --year '2017' --systematic SYST --channel 'Zll' --even
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June19_unblinding/haddjobs/ -c SR_low_Zee -p "usingBEnriched&&hJets_btagged_0>0.1522&&hJets_btagged_1>0.1522&&twoResolvedJets&&isZee&&controlSample==0&&V_pt>50&&V_pt<=150&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s "" -v "1-CMS_vhbb_DNN_Zll_LowPT_13TeV" --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -bb 0.,0.6392246969551242,0.7166787936439221,0.7663121544116362,0.8036180776631824,0.8338120434829465,0.8593296414549889,0.8815177246989334,0.9012041244819916,0.9189361786105196,0.9350957052026273,0.9499602473181461,0.963738240490261,0.9765904263924049,0.9886435218563179,1. -w 1.0  -o $ORIG_DIR/hists_SR_low_Zee_SAMP_ISYS.root -r 0 -sa SAMP --year '2017' --systematic SYST --channel 'Zll' --even
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June19_unblinding/haddjobs/ -c SR_high_Zuu -p "usingBEnriched&&hJets_btagged_0>0.1522&&hJets_btagged_1>0.1522&&twoResolvedJets&&isZmm&&controlSample==0&&V_pt>150&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s "" -v "1-CMS_vhbb_DNN_Zll_HighPT_13TeV" --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -bb 0.,0.6388302270218088,0.7163691683600514,0.7660568984287467,0.8034036872783153,0.833630728126449,0.8591762785470664,0.8813886670275982,0.9010966316569309,0.9188481098053722,0.9350253378416581,0.9499061628514183,0.9636992486919744,0.9765655131161601,0.9886318117636758,1. -w 1.0  -o $ORIG_DIR/hists_SR_high_Zuu_SAMP_ISYS.root -r 0 -sa SAMP --year '2017' --systematic SYST --channel 'Zll' --even
python splitSamples.py -i root://cmseos.fnal.gov://store/group/lpchbb/VHbbAnalysisNtuples/2017V5_June19_unblinding/haddjobs/ -c SR_low_Zuu -p "usingBEnriched&&hJets_btagged_0>0.1522&&hJets_btagged_1>0.1522&&twoResolvedJets&&isZmm&&controlSample==0&&V_pt>50&&V_pt<=150&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s "" -v "1-CMS_vhbb_DNN_Zll_LowPT_13TeV" --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -bb 0.,0.6392246969551242,0.7166787936439221,0.7663121544116362,0.8036180776631824,0.8338120434829465,0.8593296414549889,0.8815177246989334,0.9012041244819916,0.9189361786105196,0.9350957052026273,0.9499602473181461,0.963738240490261,0.9765904263924049,0.9886435218563179,1. -w 1.0  -o $ORIG_DIR/hists_SR_low_Zuu_SAMP_ISYS.root -r 0 -sa SAMP --year '2017' --systematic SYST --channel 'Zll' --even