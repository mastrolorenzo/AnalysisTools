python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c ttbar_high_Zee -p "isZee&&controlSample==21&&V_pt>150" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_ttbar_high_Zee_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c ttbar_low_Zee -p "isZee&&controlSample==21&&V_pt>50&&V_pt<=150" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_ttbar_low_Zee_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c ttbar_high_Zuu -p "isZmm&&controlSample==21&&V_pt>150" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_ttbar_high_Zuu_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c ttbar_low_Zuu -p "isZmm&&controlSample==21&&V_pt>50&&V_pt<=150" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_ttbar_low_Zuu_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zlf_high_Zee -p "isZee&&controlSample==22&&V_pt>150" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zlf_high_Zee_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zlf_low_Zee -p "isZee&&controlSample==22&&V_pt>50&&V_pt<=150" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zlf_low_Zee_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zlf_high_Zuu -p "isZmm&&controlSample==22&&V_pt>150" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zlf_high_Zuu_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zlf_low_Zuu -p "isZmm&&controlSample==22&&V_pt>50&&V_pt<=150" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zlf_low_Zuu_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zhf_high_Zee -p "isZee&&controlSample==23&&V_pt>150" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_high_Zee_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zhf_low_Zee -p "isZee&&controlSample==23&&V_pt>50&&V_pt<=150" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_low_Zee_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zhf_high_Zuu -p "isZmm&&controlSample==23&&V_pt>150" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_high_Zuu_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zhf_low_Zuu -p "isZmm&&controlSample==23&&V_pt>50&&V_pt<=150" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_low_Zuu_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zhf_high_Zee -p "isZee&&controlSample==0&&V_pt>150&&H_mass>90&&H_mass<150" -s "" -v "CMS_vhbb_BDTG_Zll_HighPT_13TeV" --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_high_Zee_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll'
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zhf_low_Zee -p "isZee&&controlSample==0&&V_pt>50&&V_pt<=150&&H_mass>90&&H_mass<150" -s "" -v "CMS_vhbb_BDTG_Zll_LowPT_13TeV" --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_low_Zee_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll'
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zhf_high_Zuu -p "isZmm&&controlSample==0&&V_pt>150&&H_mass>90&&H_mass<150" -s "" -v "CMS_vhbb_BDTG_Zll_HighPT_13TeV" --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_high_Zuu_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll'
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zhf_low_Zuu -p "isZmm&&controlSample==0&&V_pt>50&&V_pt<=150&&H_mass>90&&H_mass<150" -s "" -v "CMS_vhbb_BDTG_Zll_LowPT_13TeV" --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_low_Zuu_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll'
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zhf_high_Zee -p "isZee&&controlSample==0&&V_pt>150&&H_mass>90&&H_mass<150" -s "" -v "CMS_vhbb_DNN_Zll_HighPT_13TeV" --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_high_ZeeDNN_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll'
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zhf_low_Zee -p "isZee&&controlSample==0&&V_pt>50&&V_pt<=150&&H_mass>90&&H_mass<150" -s "" -v "CMS_vhbb_DNN_Zll_LowPT_13TeV" --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_low_ZeeDNN_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll'
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zhf_high_Zuu -p "isZmm&&controlSample==0&&V_pt>150&&H_mass>90&&H_mass<150" -s "" -v "CMS_vhbb_DNN_Zll_HighPT_13TeV" --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_high_ZuuDNN_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll'
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c Zhf_low_Zuu -p "isZmm&&controlSample==0&&V_pt>50&&V_pt<=150&&H_mass>90&&H_mass<150" -s "" -v "CMS_vhbb_DNN_Zll_LowPT_13TeV" --xlow -1.0 --xhigh 1.0  -s systematics_Zllshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_Zhf_low_ZuuDNN_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Zll'

