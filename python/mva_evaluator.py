from kinfitter import EventProxy, prepare_io
import ctypes
import time
import ROOT


def prep_bdt_variables(bdt_info):
    bdt_info = ROOT.BDTInfo(bdt_info)  # copy bdt info instances => fresh TMVA reader
    raw_vars = {}
    for bdtvar in bdt_info.bdtVars:
        raw_var = ctypes.c_float()
        raw_vars[bdtvar.localVarName] = raw_var
        if bdtvar.isSpec:
            bdt_info.reader.AddSpectator(bdtvar.varName, raw_var)
        else:
            bdt_info.reader.AddVariable(bdtvar.varName, raw_var)
    bdt_info.raw_vars = raw_vars  # just monkey patch the variables for filling later
    bdt_info.BookMVA()
    return bdt_info


def evaluate_bdt(ep, bdt_infos):
    for bdt_info in bdt_infos:
        for varname, raw_var in bdt_info.raw_vars.iteritems():
            raw_var.value = getattr(ep, varname)

        getattr(ep, bdt_info.bdtname)[0] = bdt_info.reader.EvaluateMVA(bdt_info.methodName)


def apply_mva_eval(input_file, output_file, am):

    bdt_infos = list(
        prep_bdt_variables(bdt_info)
        for _, bdt_info in am.bdtInfos
    )

    f, t, of = prepare_io(input_file, output_file)
    ot = t.CloneTree(0)

    ep = EventProxy(postfix = '_postMVA')  # FIXME remove postfix for the final thing!
    ep.init_output(ot, set(str(bdt_info.bdtname) for bdt_info in bdt_infos))
    n_events = 0

    print 'starting mva evaluation loop'
    start_time = time.ctime()
    for e in t:
        n_events += 1

        if n_events % 1000 == 1:
            print 'at event %i' % n_events

        ep.set_event(e)
        evaluate_bdt(ep, bdt_infos)
        ot.Fill()

    ot.Write()
    of.Close()
    del f
    print 'started mva evaluation loop at', start_time
    print 'finished mva evaluation loop at', time.ctime()
    print 'total number of events processed:    ', n_events
