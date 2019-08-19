python splitSamples.py -i /nfs/dust/cms/user/dewita/VHbbAnalysisNtuples/vhbb_2017_1802/haddjobs/ -c SR_high_Zee -p "usingBEnriched&&hJets_btagged_0>0.1522&&hJets_btagged_1>0.1522&&twoResolvedJets&&isZee&&controlSample==0&&V_pt>250&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s "" -v "(TMath::Sqrt(1-CMS_vhbb_DNN_Zll_HighPT_13TeV)+TMath::Power((1-CMS_vhbb_DNN_Zll_HighPT_13TeV),12))/2." --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017STXS.txt  -d 1  -w 1.0  -o $ORIG_DIR/hists_SR_high_Zee_SAMP.root -r 0 -sa SAMP --year '2017'  --channel 'Zll' --even  -n 20
python splitSamples.py -i /nfs/dust/cms/user/dewita/VHbbAnalysisNtuples/vhbb_2017_1802/haddjobs/ -c SR_low_Zee -p "usingBEnriched&&hJets_btagged_0>0.1522&&hJets_btagged_1>0.1522&&twoResolvedJets&&isZee&&controlSample==0&&V_pt>75&&V_pt<=150&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s "" -v "(TMath::Sqrt(1-CMS_vhbb_DNN_Zll_LowPT_13TeV)+TMath::Power((1-CMS_vhbb_DNN_Zll_LowPT_13TeV),12))/2." --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017STXS.txt  -d 1 -w 1.0  -o $ORIG_DIR/hists_SR_low_Zee_SAMP.root -r 0 -sa SAMP --year '2017'  --channel 'Zll' --even -n 20
python splitSamples.py -i /nfs/dust/cms/user/dewita/VHbbAnalysisNtuples/vhbb_2017_1802/haddjobs/ -c SR_high_Zmm -p "usingBEnriched&&hJets_btagged_0>0.1522&&hJets_btagged_1>0.1522&&twoResolvedJets&&isZmm&&controlSample==0&&V_pt>250&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s "" -v "(TMath::Sqrt(1-CMS_vhbb_DNN_Zll_HighPT_13TeV)+TMath::Power((1-CMS_vhbb_DNN_Zll_HighPT_13TeV),12))/2." --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017STXS.txt  -d 1  -w 1.0  -o $ORIG_DIR/hists_SR_high_Zmm_SAMP.root -r 0 -sa SAMP --year '2017'  --channel 'Zll' --even -n 20
python splitSamples.py -i /nfs/dust/cms/user/dewita/VHbbAnalysisNtuples/vhbb_2017_1802/haddjobs/ -c SR_low_Zmm -p "usingBEnriched&&hJets_btagged_0>0.1522&&hJets_btagged_1>0.1522&&twoResolvedJets&&isZmm&&controlSample==0&&V_pt>75&&V_pt<=150&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s "" -v "(TMath::Sqrt(1-CMS_vhbb_DNN_Zll_LowPT_13TeV)+TMath::Power((1-CMS_vhbb_DNN_Zll_LowPT_13TeV),12))/2." --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017STXS.txt  -d 1  -w 1.0  -o $ORIG_DIR/hists_SR_low_Zmm_SAMP.root -r 0 -sa SAMP --year '2017'  --channel 'Zll' --even -n 20
python splitSamples.py -i /nfs/dust/cms/user/dewita/VHbbAnalysisNtuples/vhbb_2017_1802/haddjobs/ -c SR_med_Zee_0j -p "usingBEnriched&&hJets_btagged_0>0.1522&&hJets_btagged_1>0.1522&&twoResolvedJets&&isZee&&controlSample==0&&V_pt>150&&V_pt<=250&&nAddJets302p5_puid<1&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s "" -v "(TMath::Sqrt(1-CMS_vhbb_DNN_Zll_HighPT_13TeV)+TMath::Power((1-CMS_vhbb_DNN_Zll_HighPT_13TeV),12))/2." --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017STXS.txt  -d 1  -w 1.0  -o $ORIG_DIR/hists_SR_med_Zee_0j_SAMP.root -r 0 -sa SAMP --year '2017'  --channel 'Zll' --even -n 20
python splitSamples.py -i /nfs/dust/cms/user/dewita/VHbbAnalysisNtuples/vhbb_2017_1802/haddjobs/ -c SR_med_Zee_ge1j -p "usingBEnriched&&hJets_btagged_0>0.1522&&hJets_btagged_1>0.1522&&twoResolvedJets&&isZee&&controlSample==0&&V_pt>150&&V_pt<=250&&nAddJets302p5_puid>0&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s "" -v "(TMath::Sqrt(1-CMS_vhbb_DNN_Zll_HighPT_13TeV)+TMath::Power((1-CMS_vhbb_DNN_Zll_HighPT_13TeV),12))/2." --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017STXS.txt  -d 1  -w 1.0  -o $ORIG_DIR/hists_SR_med_Zee_ge1j_SAMP.root -r 0 -sa SAMP --year '2017'  --channel 'Zll' --even -n 20
python splitSamples.py -i /nfs/dust/cms/user/dewita/VHbbAnalysisNtuples/vhbb_2017_1802/haddjobs/ -c SR_med_Zmm_0j -p "usingBEnriched&&hJets_btagged_0>0.1522&&hJets_btagged_1>0.1522&&twoResolvedJets&&isZmm&&controlSample==0&&V_pt>150&&V_pt<=250&&nAddJets302p5_puid<1&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s "" -v "(TMath::Sqrt(1-CMS_vhbb_DNN_Zll_HighPT_13TeV)+TMath::Power((1-CMS_vhbb_DNN_Zll_HighPT_13TeV),12))/2." --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017STXS.txt  -d 1  -w 1.0  -o $ORIG_DIR/hists_SR_med_Zmm_0j_SAMP.root -r 0 -sa SAMP --year '2017'  --channel 'Zll' --even -n 20
python splitSamples.py -i /nfs/dust/cms/user/dewita/VHbbAnalysisNtuples/vhbb_2017_1802/haddjobs/ -c SR_med_Zmm_ge1j -p "usingBEnriched&&hJets_btagged_0>0.1522&&hJets_btagged_1>0.1522&&twoResolvedJets&&isZmm&&controlSample==0&&V_pt>150&&V_pt<=250&&nAddJets302p5_puid>0&&H_mass_fit_fallback>90&&H_mass_fit_fallback<150" -s "" -v "(TMath::Sqrt(1-CMS_vhbb_DNN_Zll_HighPT_13TeV)+TMath::Power((1-CMS_vhbb_DNN_Zll_HighPT_13TeV),12))/2." --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017STXS.txt  -d 1  -w 1.0  -o $ORIG_DIR/hists_SR_med_Zmm_ge1j_SAMP.root -r 0 -sa SAMP --year '2017'  --channel 'Zll' --even -n 20
