from kinfitter import EventProxy, prepare_io, lep_sys_names
import TensorflowEvaluator
import ctypes
import time
import ROOT


VERBOSE = False
SYS_VARS = (
    'H_mass',
    'H_pt',
    'V_pt',
    'nAddJets252p9_puid',
    'HJ1_HJ2_dEta',
    'HJ1_HJ2_dPhi',
    'HJ1_HJ2_dR',
    'HVdPhi',
    'HVdR',
    'MET_pt',
    'Top1_mass_fromLepton_regPT_w4MET',
    'V_mt',
    'hJets_btagged_0',
    'hJets_btagged_1',
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


def make_sys_event_proxies(am):
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
    elif bdt_info.mvaType == "DNN":
        tokens = str(bdt_info.xmlFile).split(':')
        assert len(tokens)==3, 'for DNN, I need bdt_info.xmlFile to be in the form "cfg:scaler_dump:chpt" (separated by colons)'
        cfg, scaler_dump, chpt = tokens
        bdt_info.tf_evaluator = TensorflowEvaluator.TensorflowEvaluator(
            cfg, scaler_dump, chpt, n_features=len(bdt_info.raw_vars), verbose=VERBOSE)
    else:
        raise RuntimeError("I don't know mvaType %s. Just BDT and DNN." % bdt_info.mvaType)
    return bdt_info


def evaluate_discriminator(ep, bdt_infos):
    for bdt_info in bdt_infos:
        if bdt_info.mvaType == "BDT":
            for varname, raw_var in bdt_info.raw_vars:
                raw_var.value = getattr(ep, varname)
            getattr(ep, bdt_info.bdtname)[0] = bdt_info.reader.EvaluateMVA(bdt_info.methodName)

        else:
            raw_vars = tuple(getattr(ep, varname) for varname, _ in bdt_info.raw_vars)
            getattr(ep, bdt_info.bdtname)[0] = bdt_info.tf_evaluator.EvaluateDNN(raw_vars)


def apply_mva_eval(input_file, output_file, am):

    bdt_infos = list(
        prep_bdt_variables(bdt_info)
        for _, bdt_info in am.bdtInfos
    )

    f, t, of = prepare_io(input_file, output_file)
    ot = t.CloneTree(0)

    eps = make_sys_event_proxies(am) + [EventProxy()]  # the last one is for nominal
    for ep in eps:
        ep.init_output(ot, set(str(bdt_info.bdtname) for bdt_info in bdt_infos))
    n_events = 0

    print 'starting mva evaluation loop with %i systematics' % (len(eps) - 1)
    start_time = time.ctime()
    for e in t:
        n_events += 1

        if n_events % 1000 == 1:
            print 'at event %i' % n_events

        for ep in eps:
            ep.set_event(e)
            evaluate_discriminator(ep, bdt_infos)
        ot.Fill()

    ot.Write()
    of.Close()
    del f
    print 'started mva evaluation loop at', start_time
    print 'finished mva evaluation loop at', time.ctime()
    print 'total number of events processed:    ', n_events
