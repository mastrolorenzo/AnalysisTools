python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c ttWmn -p "isWmunu&&controlSample==11" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Wlnshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_ttWmn_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Wln' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c ttWen -p "isWenu&&controlSample==11" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0 -s systematics_Wlnshapes2017.txt   -d 1 -n 15 -w 1.0 -o $ORIG_DIR/hists_ttWen_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Wln' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c whfWmnLow -p "isWmunu&&H_mass<90&&controlSample==13" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0 -s systematics_Wlnshapes2017.txt   -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_whfWmnLow_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Wln' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c whfWenLow -p "isWenu&&H_mass<90&&controlSample==13" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0 -s systematics_Wlnshapes2017.txt   -d 1 -n 15 -w 1.0 -o $ORIG_DIR/hists_whfWenLow_SAMP.root -r 0 -sa SAMP --year'2017' --channel 'Wln' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c whfWmnHigh -p "isWmunu&&H_mass>150&&controlSample==13" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0  -s systematics_Wlnshapes2017.txt  -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_whfWmnHigh_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Wln' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c whfWenHigh -p "isWenu&&H_mass>150&&controlSample==13" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0 -s systematics_Wlnshapes2017.txt   -d 1 -n 15 -w 1.0 -o $ORIG_DIR/hists_whfWenHigh_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Wln' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c wlfWmn -p "isWmunu&&controlSample==12" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0 -s systematics_Wlnshapes2017.txt   -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_wlfWmn_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Wln' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c wlfWen -p "isWenu&&controlSample==12" -s "" -v "Jet_btagDeepB[hJetInd2]" --xlow 0.0 --xhigh 1.0 -s systematics_Wlnshapes2017.txt   -d 1 -n 15 -w 1.0 -o $ORIG_DIR/hists_wlfWen_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Wln' --drawFromNom True
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c WmnHighPt -p "isWmunu&&controlSample==0&&H_mass>90&&H_mass<150" -s "" -v "CMS_vhbb_BDTG_Wln_13TeV" --xlow -1.0 --xhigh 1.0 -s systematics_Wlnshapes2017.txt   -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_WmnHighPt_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Wln'
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c WenHighPt -p "isWenu&&controlSample==0&&H_mass>90&&H_mass<150" -s "" -v "CMS_vhbb_BDTG_Wln_13TeV" --xlow -1.0 --xhigh 1.0 -s systematics_Wlnshapes2017.txt   -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_WenHighPt_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Wln'
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c WmnHighPt -p "isWmunu&&controlSample==0&&H_mass>90&&H_mass<150" -s "" -v "CMS_vhbb_DNN_Wmn_13TeV" --xlow 0 --xhigh 1.0 -s systematics_Wlnshapes2017.txt   -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_WmnHighPtDNN_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Wln'
python splitSamples.py -i /nfs/dust/cms/user/dewita/hadd2016V4_fromChris_2/haddjobs/ -c WenHighPt -p "isWenu&&controlSample==0&&H_mass>90&&H_mass<150" -s "" -v "CMS_vhbb_DNN_Wen_13TeV" --xlow 0 --xhigh 1.0 -s systematics_Wlnshapes2017.txt   -d 1 -n 15 -w 1.0  -o $ORIG_DIR/hists_WenHighPtDNN_SAMP.root -r 0 -sa SAMP --year '2017' --channel 'Wln'
