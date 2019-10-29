from kinfitter import EventProxy, prepare_io
import ctypes
import time
import ROOT
import numpy

VERBOSE = False
SYS_VARS = (
    'H_mass',
    'H_pt',
    'V_pt',
    'nAddJets252p9_puid',
    'HJ1_HJ2_dEta',
    'HJ1_HJ2_dPhi',
    'abs_HJ1_HJ2_dPhi',
    'HJ1_HJ2_dR',
    'HVdPhi',
    'HVdR',
    'MET_pt',
    'Top1_mass_fromLepton_regPT_w4MET',
    'V_mt',
    'hJets_btagged_0',
    'hJets_btagged_1',
    'hJets_btagWP_0',
    'hJets_btagWP_1',
    'hJets_leadingPt',
    'hJets_pt_0',
    'hJets_pt_1',
    'hJets_subleadingPt',
    'jjVPtRatio',
    'lepMetDPhi',
    'minDPhiFromOtherJets',
    'nAddJet_f',
    'nAddJets302p5_puid_',
    'nAddJets252p5_puid',
    'nAddJets_2lep',
    'nJets25_dR06',
    'nJets30_0lep',
    'nJets30_2lep',
    'nLooseBtagsDR0p8',
    'nLooseBtagsDR1p0',
    'otherJetsBestBtag',
    'otherJetsHighestP',

    # kin-fit only variables
    'H_mass_fit_fallback',
    'H_pt_fit_fallback',
    'jjVPtRatio_fit_fallback',
    'HVdPhi_fit_fallback',
    'n_recoil_jets_fit',
    'H_mass_sigma_fit',
)


def make_sys_event_proxies(am, lep_sys_names):
    def mk_single_getter(branch, sys_name):
        # this function gives a local namespace to branch, sys_name
        return lambda e: getattr(e, branch+'_'+sys_name, getattr(e, branch))

    def mk_sys_attr_getters(sys_name):
        return dict(
            (sys_var, mk_single_getter(sys_var, sys_name))
            for sys_var in SYS_VARS
        )

    sys_names = lep_sys_names + list(
        amsys.name
        for amsys in am.systematics
        if str(amsys.name) != 'nominal'
    )
    return list(
        EventProxy(sys_name, mk_sys_attr_getters(sys_name))
        for sys_name in sys_names
    )


def prep_bdt_variables(bdt_info):
    bdt_info = ROOT.BDTInfo(bdt_info)  # copy bdt info instances => fresh TMVA reader

    def mk_raw_vars(bdtvar):
        raw_var = ctypes.c_float()
        if bdtvar.isSpec:
            bdt_info.reader.AddSpectator(bdtvar.varName, raw_var)
        else:
            bdt_info.reader.AddVariable(bdtvar.varName, raw_var)
        return bdtvar.localVarName, raw_var

    bdt_info.raw_vars = tuple(  # just monkey patch the variables for filling later (it's a tuple of tuples)
        mk_raw_vars(bdtvar)
        for bdtvar in bdt_info.bdtVars
    )

    if bdt_info.mvaType == "BDT":
        bdt_info.BookMVA()
    elif bdt_info.mvaType == "DNN" or bdt_info.mvaType == "MultiDNN":
        tokens = str(bdt_info.xmlFile).split(':')
        if len(tokens)==3: #old evaluation 'for DNN, I need bdt_info.xmlFile to be in the form "cfg:scaler_dump:checkPointFile" (separated by colons)'
            import TensorflowEvaluatorSummer2018
            cfg, scaler_dump, checkPointFile = tokens
            bdt_info.tf_evaluator = TensorflowEvaluatorSummer2018.TensorflowEvaluator(
                cfg, scaler_dump, checkPointFile, n_features=len(bdt_info.raw_vars), verbose=VERBOSE)
        elif len(tokens)==1:  #new evaluation only needs the checkpoint file
            import TensorflowEvaluatorRun2Legacy
            checkPointFile = tokens[0]
            bdt_info.tf_evaluator = TensorflowEvaluatorRun2Legacy.TensorflowDNNEvaluator(checkPointFile)
        else:
            print "len(tokens)",len(tokens),"is not 1 or 3.  So this is over"
            assert(0==1)
    else:
        raise RuntimeError("I don't know mvaType %s. Just BDT and DNN." % bdt_info.mvaType)
    return bdt_info


def evaluate_discriminator(ep, bdt_infos):
    for bdt_info in bdt_infos:
        if bdt_info.mvaType == "BDT":
            for varname, raw_var in bdt_info.raw_vars:
                raw_var.value = getattr(ep, varname)
            getattr(ep, bdt_info.bdtname)[0] = bdt_info.reader.EvaluateMVA(bdt_info.methodName)

        elif bdt_info.mvaType == "DNN":
            raw_vars = tuple(getattr(ep, varname) for varname, _ in bdt_info.raw_vars)
            getattr(ep, bdt_info.bdtname)[0] = bdt_info.tf_evaluator.EvaluateDNN(raw_vars)
        elif bdt_info.mvaType == "MultiDNN":
            raw_vars = tuple(getattr(ep, varname) for varname, _ in bdt_info.raw_vars)
            multiDNNOutput = bdt_info.tf_evaluator.EvaluateMultiDNN(raw_vars)
            for iOut in range(len(multiDNNOutput)):
                getattr(ep, bdt_info.bdtname)[0] = multiDNNOutput[iOut]
            getattr(ep, bdt_info.mostProbIndex)[0] = float(numpy.where(multiDNNOutput==max(multiDNNOutput))[0][0])
        else:
            print "What type of MVA is",bdt_info.mvaType


def apply_mva_eval(input_file, output_file, am, is_data, allowed_names, lep_sys_names=tuple()):

    bdt_infos = list(
        prep_bdt_variables(bdt_info)
        for bdt_name, bdt_info in am.bdtInfos
        if bdt_name in allowed_names
    )

    f, t, of = prepare_io(input_file, output_file)
    ot = t.CloneTree(0)

    if is_data:
        eps = [EventProxy()]
    else:
        eps = make_sys_event_proxies(am, lep_sys_names) + [EventProxy()]  # the last one is for nominal

    for ep in eps:
        listOfOutputBranches=[str(bdt_info.bdtname) for bdt_info in bdt_infos]
        for bdt_info in bdt_infos:
            if bdt_info.mvaType == "MultiDNN":
                listOfOutputBranches.append(str(bdt_info.mostProbIndex))
        ep.init_output(ot, set(listOfOutputBranches))

    n_events = 0
    start_time = time.ctime()
    print 'starting mva evaluation loop with %i systematics on %s' % (len(eps) - 1, start_time)
    for e in t:
        n_events += 1

        if n_events % 1000 == 1:
            print 'at event %i' % n_events

        for ep in eps:
            ep.set_event(e)
            #if ep.controlSample != 0 or not ep.twoResolvedJets: continue
            evaluate_discriminator(ep, bdt_infos)
        ot.Fill()

    ot.Write()
    of.Close()
    del f
    print 'started mva evaluation loop at', start_time
    print 'finished mva evaluation loop at', time.ctime()
    print 'total number of events processed:    ', n_events
