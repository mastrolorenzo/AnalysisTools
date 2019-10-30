import TensorflowEvaluatorRun2Legacy
import sys

checkPointFile=sys.argv[1]
print checkPointFile

# H_mass *      H_pt *      V_mt *      V_pt * V_pt/H_pt * abs(TVect * Jet_btagD * Jet_btagD * max(Jet_P * min(Jet_P * abs(Jet_e *    MET_Pt * dPhiLepMe * top_mass2 *       SA5 * nAddJets3
#H_mass H_pt V_mt V_pt V_pt/H_pt abs(TVector2::Phi_mpi_pi(V_phi-H_phi)) Jet_btagDeepB[hJidx[0]] Jet_btagDeepB[hJidx[1]] max(Jet_PtReg[hJidx[0]],Jet_PtReg[hJidx[1]]) min(Jet_PtReg[hJidx[0]],Jet_PtReg[hJidx[1]]) abs(Jet_eta[hJidx[0]]-Jet_eta[hJidx[1]]) MET_Pt dPhiLepMet top_mass2_05 SA5 nAddJets302p5_puid
tf_evaluator = TensorflowEvaluatorRun2Legacy.TensorflowDNNEvaluator(checkPointFile)

raw_vars=[133.42915 , 215.59137 , 61.906452 , 251.24882 , 1.1653937 , 2.7094084 , 0.9951171, 0.9433593 , 143.99345 , 71.707962 , 1.2127075 , 194.36132 , 0.5606406 , 384.10931, 1, 0]
output=tf_evaluator.EvaluateDNN(raw_vars)
print "correct inputs",output

inversion=1./1.1653937
raw_vars=[133.42915 , 215.59137 , 61.906452 , 251.24882 , inversion , 2.7094084 , 0.9951171, 0.9433593 , 143.99345 , 71.707962 , 1.2127075 , 194.36132 , 0.5606406 , 384.10931, 1, 0]
output=tf_evaluator.EvaluateDNN(raw_vars)
print "inverted jjpt/vpt",output

raw_vars=[134.8 , 215.59137 , 61.906452 , 251.24882 , inversion , 2.7094084 , 0.9951171, 0.9433593 , 143.99345 , 71.707962 , 1.2127075 , 194.36132 , 0.5606406 , 384.10931, 1, 0]
output=tf_evaluator.EvaluateDNN(raw_vars)
print "inverted jjpt/vpt, AT mjj",output

raw_vars=[134.8 , 215.59137 , 61.906452 , 251.24882 , inversion , 2.7094084 , 0.9951171, 0.9433593 , 143.99345 , 71.707962 , 1.2127075 , 194.36132 , 0.5606406 , 394.1, 1, 0]
output=tf_evaluator.EvaluateDNN(raw_vars)
print "inverted jjpt/vpt, AT mjj, AT top mass",output

