for channel in Znn_13TeV_Signal Znn_13TeV_TT Znn_13TeV_Zbb Znn_13TeV_Zlight SR_high_Zee SR_low_Zee SR_high_Zuu SR_low_Zuu ttbar_high_Zee ttbar_low_Zee ttbar_high_Zuu ttbar_low_Zuu Zlf_high_Zee Zlf_low_Zee Zlf_high_Zuu Zlf_low_Zuu Zhf_high_Zee Zhf_low_Zee Zhf_high_Zuu Zhf_low_Zuu ttWmn ttWen wlfWmn wlfWen whfWmnLow whfWmnHigh whfWenLow whfWenHigh WmnHighPt WenHighPt
do
    hadd -f hists_${channel}.root hists_${channel}_*.root
done
hadd -f vhbb_Znn-2016.root hists_Znn_13TeV_Signal.root hists_Znn_13TeV_Zlight.root hists_Znn_13TeV_Zbb.root hists_Znn_13TeV_TT.root  
hadd -f vhbb_Zee-2016.root hists_SR_low_Zee.root hists_SR_high_Zee.root hists_ttbar_low_Zee.root hists_ttbar_high_Zee.root hists_Zlf_low_Zee.root hists_Zlf_high_Zee.root hists_Zhf_low_Zee.root hists_Zhf_high_Zee.root
hadd -f vhbb_Zmm-2016.root hists_SR_low_Zuu.root hists_SR_high_Zuu.root hists_ttbar_low_Zuu.root hists_ttbar_high_Zuu.root hists_Zlf_low_Zuu.root hists_Zlf_high_Zuu.root hists_Zhf_low_Zuu.root hists_Zhf_high_Zuu.root
hadd -f vhbb_Wen-2016.root hists_ttWen.root hists_wlfWen.root hists_whfWenHigh.root hists_whfWenLow.root hists_WenHighPt.root
hadd -f vhbb_Wmn-2016.root hists_ttWmn.root hists_wlfWmn.root hists_whfWmnHigh.root hists_whfWmnLow.root hists_WmnHighPt.root
