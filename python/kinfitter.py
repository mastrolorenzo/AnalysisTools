from ROOT.TMath import Sin, Cos, Pi
import numpy
import time
import ROOT


# constants
PI = Pi()
LLVV_PXY_VAR = 8**2  # with no additional jets above pt > 20 GeV and no eta cut
HH4B_RES_SCALE = 0.62
Z_WIDTH = 1.7 * 3
TREE_VARS = {
    # higgs jets
    'n_hj_matched',
    'hj1_pt',        'hj1_eta',        'hj1_phi',        'hj1_mass',
    'hj2_pt',        'hj2_eta',        'hj2_phi',        'hj2_mass',
    'hj12_pt',       'hj12_eta',       'hj12_phi',       'hj12_mass',
    'hj1_gen_pt',    'hj1_gen_eta',    'hj1_gen_phi',    'hj1_gen_mass',
    'hj2_gen_pt',    'hj2_gen_eta',    'hj2_gen_phi',    'hj2_gen_mass',
    'hj12_gen_pt',   'hj12_gen_eta',   'hj12_gen_phi',   'hj12_gen_mass',
    'hj1_reg_pt',    'hj1_reg_eta',    'hj1_reg_phi',    'hj1_reg_mass',
    'hj2_reg_pt',    'hj2_reg_eta',    'hj2_reg_phi',    'hj2_reg_mass',
    'hj12_reg_pt',   'hj12_reg_eta',   'hj12_reg_phi',   'hj12_reg_mass',
    'hj1_fit_pt',    'hj1_fit_eta',    'hj1_fit_phi',    'hj1_fit_mass',
    'hj2_fit_pt',    'hj2_fit_eta',    'hj2_fit_phi',    'hj2_fit_mass',
    'H_pt_fit',      'H_eta_fit',      'H_phi_fit',      'H_mass_fit',
    'hJets_pt_0_fit',
    'hJets_pt_1_fit',
    'H_mass_sigma_fit',
    'n_fsr_jets',
    'n_recoil_jets_fit',

    'HVdPhi_fit',    'HVdR_fit',       'HVdEta_fit',      'jjVPtRatio_fit',

    'hj1_reg_res',
    'hj2_reg_res',
    'hj1_hh4b_res',
    'hj2_hh4b_res',
    'hj1_jme_res',
    'hj2_jme_res',
    'hj1_fit_err',
    'hj2_fit_err',
    'hj12_fit_corr',

    # leptons
    'l1_pt',        'l1_eta',        'l1_phi',        'l1_mass',
    'l2_pt',        'l2_eta',        'l2_phi',        'l2_mass',
    'V_eta',
    'l1_gen_pt',    'l1_gen_eta',    'l1_gen_phi',    'l1_gen_mass',
    'l2_gen_pt',    'l2_gen_eta',    'l2_gen_phi',    'l2_gen_mass',
    'l12_gen_pt',   'l12_gen_eta',   'l12_gen_phi',   'l12_gen_mass',
    'l1_fit_pt',    'l1_fit_eta',    'l1_fit_phi',    'l1_fit_mass',
    'l2_fit_pt',    'l2_fit_eta',    'l2_fit_phi',    'l2_fit_mass',
    'V_pt_fit',     'V_eta_fit',     'V_phi_fit',     'V_mass_fit',

    'l1_res',
    'l2_res',

    # llbb system / recoil jet
    # 'recoil_pt',     'recoil_eta',     'recoil_phi',     'recoil_mass',
    # 'recoil_fit_pt', 'recoil_fit_eta', 'recoil_fit_phi', 'recoil_fit_mass',

    # 'llbb_px',       'llbb_py',
    'llbb_fit_px',   'llbb_fit_py',   'llbb_fit_pt',
    'llbbr_fit_px',  'llbbr_fit_py',  'llbbr_fit_pt',
    'llbb_gen_px',   'llbb_gen_py',
    'll_gen_mass',

    # aux
    'kinfit_fit',
    'kinfit_getNDF',
    'kinfit_getS',
    'kinfit_getF',

    'H_mass_fit_fallback',
    'H_pt_fit_fallback',
    'HVdPhi_fit_fallback',
    'jjVPtRatio_fit_fallback',
    'hJets_pt_0_fit_fallback',
    'hJets_pt_1_fit_fallback',
}

# resolution functions from HH4b analysis
# https://github.com/cvernier/HH4b2016/blob/master/HbbHbb_Component_KinFit.cc
def ErrpT_Signal(pT, eta):
    if abs(eta)<1.4:
        if pT<40: pT=40
        if pT>550: pT=550
        sigmapT = 22.26 - 0.01*pT + 0.00018*pT*pT
    else:
        if pT<40: pT=40
        if pT>350: pT=350
        sigmapT = 17.11 + 0.058 * pT
    return sigmapT*sigmapT

def ErrEta_Signal(pT):
    if pT<40: pT=40
    if pT>500: pT=500
    sigmaEta = 0.033 + (4.1/pT) + (-0.17/pow(pT,0.5))
    return sigmaEta*sigmaEta

def ErrPhi_Signal(pT):
    if pT<40: pT=40
    if pT>500: pT=500
    sigmaPhi = 0.038 + (4.1/pT) + (-0.19/pow(pT,0.5))
    return sigmaPhi*sigmaPhi


# auxillary functions
def deltaPhi(phi1, phi2):
    dphi = (phi1 - phi2 + PI) % (2*PI) - PI
    return dphi

def deltaR(eta1, phi1, eta2, phi2):
    dphi = deltaPhi(phi1, phi2)
    deta = eta1 - eta2
    return (deta*deta + dphi*dphi)**.5

def matrix_add(self, a):
    mtrx = ROOT.TMatrixD(self.GetNrows(), self.GetNcols())
    for i in xrange(self.GetNrows()):
        for j in xrange(self.GetNcols()):
            mtrx[i][j] = self[i][j] + a[i][j]
    return mtrx

ROOT.TMatrixD.__add__ = matrix_add

def mk_lv(pt, eta, phi, m):
    v = ROOT.TLorentzVector()
    v.SetPtEtaPhiM(pt, eta, phi, m)
    return v


# functions for covariances and fit particles from jets ...
def mk_lv_ind(e, ind, regressed=False):
    return mk_lv(
        e.Jet_PtReg[ind] if regressed else e.Jet_Pt[ind],
        e.Jet_eta[ind],
        e.Jet_phi[ind],
        e.Jet_mass[ind],
    )

def mk_jet_cart_cov(v, rho, jme_res_obj):
    jme_par = ROOT.JME.JetParameters().setRho(rho).setJetPt(v.Pt()).setJetEta(v.Eta())
    jme_res = jme_res_obj.getResolution(jme_par) * v.Pt()
    cos_phi = Cos(v.Phi())
    sin_phi = Sin(v.Phi())
    cov = ROOT.TMatrixD(4, 4)
    cov.Zero()
    cov[0][0] = (jme_res*cos_phi)**2
    cov[1][0] = jme_res**2 * cos_phi * sin_phi
    cov[0][1] = cov[1][0]
    cov[1][1] = (jme_res*sin_phi)**2
    return cov

def mk_fit_particle(e, ind, res_scale=1.):
    v = mk_lv_ind(e, ind, True)
    cov = ROOT.TMatrixD(3, 3)
    cov.Zero()
    hh4b_res = ErrpT_Signal(v.Et(), v.Eta())
    reg_res = (e.Jet_bRegRes[ind]*v.Pt())**2
    cov[0][0] = hh4b_res * res_scale**2
    #cov[0][0] = reg_res
    cov[1][1] = ErrEta_Signal(v.Et())
    cov[2][2] = ErrPhi_Signal(v.Et())
    return v, hh4b_res, reg_res, ROOT.TFitParticleEtEtaPhi(v, cov)


# ... and from leptons as well
def mk_mu_lv_ind(e, ind):
    return mk_lv(
        abs(e.Muon_pt_corrected[ind]),
        # e.Muon_pt[ind],
        e.Muon_eta[ind],
        e.Muon_phi[ind],
        e.Muon_mass[ind],
    )

def mk_el_lv_ind(e, ind):
    return mk_lv(
        e.Electron_pt[ind],
        e.Electron_eta[ind],
        e.Electron_phi[ind],
        e.Electron_mass[ind],
    )

def get_sel_leptons(e):
    if e.isZmm:
        return mk_mu_lv_ind(e, e.lepInd1), mk_mu_lv_ind(e, e.lepInd2)
    else:
        return mk_el_lv_ind(e, e.lepInd1), mk_el_lv_ind(e, e.lepInd2)

def mk_fit_particle_lepton(e, v, ind):
    cov = ROOT.TMatrixD(3, 3)
    cov.Zero()
    cov[1][1] = 0.0001 # eta: too small
    cov[2][2] = 0.0001 # phi: too small

    if e.isZmm:
        cov[0][0] = (e.Muon_ptErr[ind])**2
    else:
        cov[0][0] = (e.Electron_energyErr[ind])**2

    return v, cov[0][0], ROOT.TFitParticleEtEtaPhi(v, cov)


class EventProxy(object):
    '''
    Arguments:
    - postfix is added as to the name of all output branches
    - sys_branch (optional) is the name of the branch that should be replaced
    - sys_attr_getter (optional) is a function that takes the event as an argument and returns
      the systematic replacement.
    '''
    def __init__(self, postfix='', sys_branch='', sys_attr_getter=None):
        assert bool(sys_branch) == bool(callable(sys_attr_getter)
            ), 'I need both, sys_branch and a callable sys_attr_getter, or none of them.'
        self.output_postfix = postfix
        self.sys_branch = sys_branch
        self.sys_attr_getter = sys_attr_getter
        self.e = None
        self.cache = None
        self.tree_vars = set()

    def init_output(self, output_tree, tree_vars):
        postfix = '_'+self.output_postfix if self.output_postfix else ''
        self.tree_vars |= tree_vars
        for var in self.tree_vars:
            setattr(self, var, numpy.zeros(1,float))
            output_tree.Branch(var+postfix, getattr(self, var), var+postfix+'/D')
        self.set_vals_to_zero()

    def set_event(self, e):
        self.e = e
        self.cache = None

    def apply_fallback(self):
        self.H_mass_fit_fallback[0] = self.e.H_mass
        self.H_pt_fit_fallback[0] = self.e.H_pt
        self.HVdPhi_fit_fallback[0] = self.e.HVdPhi
        self.jjVPtRatio_fit_fallback[0] = self.e.jjVPtRatio
        if self.e.twoResolvedJets:
            self.hJets_pt_0_fit_fallback[0] = self.e.Jet_PtReg[self.e.hJetInd1]
            self.hJets_pt_1_fit_fallback[0] = self.e.Jet_PtReg[self.e.hJetInd2]

        self.H_mass_sigma_fit[0] = -1
        self.n_recoil_jets_fit[0] = -1

    def set_vals_to_zero(self):
        for var in self.tree_vars:
            getattr(self, var)[0] = 0.

    def fill_pt_eta_phi_mass(self, token, vec):
        getattr(self, token+'_pt'  )[0] = vec.Pt()
        getattr(self, token+'_eta' )[0] = vec.Eta()
        getattr(self, token+'_phi' )[0] = vec.Phi()
        getattr(self, token+'_mass')[0] = vec.M()

    def __getattr__(self, *args):
        def make_variation():
            self.cache = self.sys_attr_getter(self.e)
            return self.cache

        if self.sys_attr_getter and args[0] == self.sys_branch:
            return self.cache or make_variation()

        return getattr(self.e, *args)


def apply_fit_to_event(ep, jme_res_obj):

    # normal higgs jets
    hJidxs = (ep.hJetInd1, ep.hJetInd2)
    v1 = mk_lv_ind(ep, hJidxs[0])
    v2 = mk_lv_ind(ep, hJidxs[1])
    v12 = v1 + v2

    ep.fill_pt_eta_phi_mass('hj1', v1)
    ep.fill_pt_eta_phi_mass('hj2', v2)
    ep.fill_pt_eta_phi_mass('hj12', v12)

    # regressed jets and fit particles
    v1, v1_hh4b_var, v1_reg_var, p1 = mk_fit_particle(ep, hJidxs[0], HH4B_RES_SCALE)
    v2, v2_hh4b_var, v2_reg_var, p2 = mk_fit_particle(ep, hJidxs[1], HH4B_RES_SCALE)

    # fsr
    jet_idxs_for_fsr = (                     # find indices
        i
        for i in xrange(ep.nJet)
        if (
            i not in hJidxs  # exclude higgs jets
            # and ep.Jet_jetId[i] > 0
            and ep.Jet_lepFilter[i] > 0
            and ep.Jet_puId[i] > 0
            and ep.Jet_Pt[i] > 20
            and abs(ep.Jet_eta[i]) < 3.0
        )
    )
    fsr_jet_idxs_vs = (                      # build Lorentz vectors
        (i, mk_lv_ind(ep, i))
        for i in jet_idxs_for_fsr
    )
    fsr_jet_idxs_vs = list(                  # filter for the ones close to hj's
        (i, v)
        for i, v in fsr_jet_idxs_vs
        if min(v.DeltaR(v1), v.DeltaR(v2)) < 0.8
    )
    for _, v in fsr_jet_idxs_vs:             # add to higgs jets
        if v.DeltaR(v1) < v.DeltaR(v2):
            v1 += v
        else:
            v2 += v
    fsr_jet_idxs = list(i for i, _ in fsr_jet_idxs_vs)
    ep.n_fsr_jets[0] = len(fsr_jet_idxs)

    v12 = v1 + v2
    ep.fill_pt_eta_phi_mass('hj1_reg', v1)
    ep.fill_pt_eta_phi_mass('hj2_reg', v2)
    ep.fill_pt_eta_phi_mass('hj12_reg', v12)

    # higgs jet matching
    ep.n_hj_matched[0] = 0
    if deltaR(ep.hj1_reg_eta[0], ep.hj1_reg_phi[0], ep.GenBJ1_eta, ep.GenBJ1_phi) < 0.2:
        v1_gen = mk_lv(ep.GenBJ1_pt,ep.GenBJ1_eta,ep.GenBJ1_phi,ep.GenBJ1_mass)
        ep.fill_pt_eta_phi_mass('hj1_gen', v1_gen)
        ep.n_hj_matched[0] += 1
    elif deltaR(ep.hj1_reg_eta[0], ep.hj1_reg_phi[0], ep.GenBJ2_eta, ep.GenBJ2_phi) < 0.2:
        v1_gen = mk_lv(ep.GenBJ2_pt,ep.GenBJ2_eta,ep.GenBJ2_phi,ep.GenBJ2_mass)
        ep.fill_pt_eta_phi_mass('hj1_gen', v1_gen)
        ep.n_hj_matched[0] += 1

    if deltaR(ep.hj2_reg_eta[0], ep.hj2_reg_phi[0], ep.GenBJ1_eta, ep.GenBJ1_phi) < 0.2:
        v2_gen = mk_lv(ep.GenBJ1_pt,ep.GenBJ1_eta,ep.GenBJ1_phi,ep.GenBJ1_mass)
        ep.fill_pt_eta_phi_mass('hj2_gen', v2_gen)
        ep.n_hj_matched[0] += 1
    elif deltaR(ep.hj2_reg_eta[0], ep.hj2_reg_phi[0], ep.GenBJ2_eta, ep.GenBJ2_phi) < 0.2:
        v2_gen = mk_lv(ep.GenBJ2_pt,ep.GenBJ2_eta,ep.GenBJ2_phi,ep.GenBJ2_mass)
        ep.fill_pt_eta_phi_mass('hj2_gen', v2_gen)
        ep.n_hj_matched[0] += 1

    if ep.n_hj_matched[0] == 2:
        v12_gen = v1_gen + v2_gen
        ep.fill_pt_eta_phi_mass('hj12_gen', v12_gen)

    # fill resolutions
    jme_p1 = ROOT.JME.JetParameters().setRho(ep.fixedGridRhoFastjetAll).setJetPt(ep.hj1_reg_pt[0]).setJetEta(ep.hj1_reg_eta[0])
    jme_p2 = ROOT.JME.JetParameters().setRho(ep.fixedGridRhoFastjetAll).setJetPt(ep.hj2_reg_pt[0]).setJetEta(ep.hj2_reg_eta[0])
    ep.hj1_jme_res[0] = jme_res_obj.getResolution(jme_p1) * ep.hj1_reg_pt[0]
    ep.hj2_jme_res[0] = jme_res_obj.getResolution(jme_p2) * ep.hj2_reg_pt[0]
    ep.hj1_hh4b_res[0] = v1_hh4b_var**.5
    ep.hj2_hh4b_res[0] = v2_hh4b_var**.5
    ep.hj1_reg_res[0] = v1_reg_var**.5
    ep.hj2_reg_res[0] = v2_reg_var**.5

    # base recoil covariance
    cov_recoil = ROOT.TMatrixD(4, 4)
    cov_recoil.Zero()
    cov_recoil[0][0] = LLVV_PXY_VAR
    cov_recoil[1][1] = LLVV_PXY_VAR
    cov_recoil[2][2] = 1
    cov_recoil[3][3] = 1e12

    # recoil vector and covariance additions
    recoil_jet_idxs = (
        i
        for i in xrange(ep.nJet)
        if (
            i not in hJidxs             # exclude higgs jets
            and i not in fsr_jet_idxs   # exclude fsr jets (if fsr is used)
            and ep.Jet_puId[i] > 1
            and ep.Jet_jetId[i] > 1
            and ep.Jet_lepFilter[i] > 0
            and ep.Jet_Pt[i] > 20
            # and abs(ep.Jet_eta) < 2.5
        )
    )
    recoil_jet_vs = list(
        mk_lv_ind(ep, i)
        for i in recoil_jet_idxs
    )
    ep.n_recoil_jets_fit[0] = len(recoil_jet_vs)
    v_recoil = sum(recoil_jet_vs, ROOT.TLorentzVector())
    cov_recoil = sum(
        (mk_jet_cart_cov(v, ep.fixedGridRhoFastjetAll, jme_res_obj) for v in recoil_jet_vs),
        cov_recoil
    )
    p_recoil = ROOT.TFitParticleECart(v_recoil, cov_recoil)

    vl1, vl2 = get_sel_leptons(ep)
    vl12 = vl1 + vl2
    _, vl1_pt_var, pl1 = mk_fit_particle_lepton(ep, vl1, ep.lepInd1)
    _, vl2_pt_var, pl2 = mk_fit_particle_lepton(ep, vl2, ep.lepInd2)

    ep.fill_pt_eta_phi_mass('l1', vl1)
    ep.fill_pt_eta_phi_mass('l2', vl2)
    ep.V_eta[0] = vl12.Eta()

    ep.l1_res[0] = vl1_pt_var**.5
    ep.l2_res[0] = vl2_pt_var**.5

    # constraints
    cons_MZ = ROOT.TFitConstraintMGaus('cons_MZ', '', None, None, 91., Z_WIDTH)
    cons_MZ.addParticle1(pl1)
    cons_MZ.addParticle1(pl2)

    cons_x = ROOT.TFitConstraintEp('cons_x', '', ROOT.TFitConstraintEp.pX)
    cons_x.addParticles(p1, p2, pl1, pl2, p_recoil)

    cons_y = ROOT.TFitConstraintEp('cons_y', '', ROOT.TFitConstraintEp.pY)
    cons_y.addParticles(p1, p2, pl1, pl2, p_recoil)

    # setup fitter
    fitter = ROOT.TKinFitter()
    fitter.addMeasParticles(p1, p2, pl1, pl2, p_recoil)
    fitter.addConstraint(cons_MZ)
    fitter.addConstraint(cons_x)
    fitter.addConstraint(cons_y)

    fitter.setMaxNbIter(30)
    fitter.setMaxDeltaS(1e-2)
    fitter.setMaxF(1e-1)
    fitter.setVerbosity(3)

    # run the fit
    ep.kinfit_fit[0] = fitter.fit() + 1  # 0: not run; 1: fit converged; 2: fit didn't converge
    ep.kinfit_getNDF[0] = fitter.getNDF()
    ep.kinfit_getS[0] = fitter.getS()
    ep.kinfit_getF[0] = fitter.getF()

    # higgs jet output vectors
    v1 = fitter.get4Vec(0)
    v2 = fitter.get4Vec(1)
    v12 = v1 + v2

    ep.fill_pt_eta_phi_mass('hj1_fit', v1)
    ep.fill_pt_eta_phi_mass('hj2_fit', v2)
    ep.hJets_pt_0_fit[0] = v1.Pt()
    ep.hJets_pt_1_fit[0] = v2.Pt()
    ep.H_pt_fit[0] = v12.Pt()
    ep.H_eta_fit[0] = v12.Eta()
    ep.H_phi_fit[0] = v12.Phi()
    ep.H_mass_fit[0] = v12.M()

    # higgs mass uncertainty
    cov_fit = fitter.getCovMatrixFit()
    ep.hj1_fit_err[0] = cov_fit(0,0)**.5
    ep.hj2_fit_err[0] = cov_fit(3,3)**.5
    ep.hj12_fit_corr[0] = cov_fit(0,3) / ep.hj1_fit_err[0] / ep.hj2_fit_err[0]
    dmH_by_dpt1 = ( (v12.X())*Cos(v1.Phi()) + (v12.Y())*Sin(v1.Phi()) )/ep.H_mass_fit[0]
    dmH_by_dpt2 = ( (v12.X())*Cos(v2.Phi()) + (v12.Y())*Sin(v2.Phi()) )/ep.H_mass_fit[0]
    ep.H_mass_sigma_fit[0] = ( dmH_by_dpt1**2          * cov_fit(0,0)
                            + dmH_by_dpt2**2          * cov_fit(3,3)
                            + dmH_by_dpt1*dmH_by_dpt2 * cov_fit(0,3) )**.5

    # lepton output vectors
    vl1 = fitter.get4Vec(2)
    vl2 = fitter.get4Vec(3)
    vl12 = vl1 + vl2

    ep.fill_pt_eta_phi_mass('l1_fit', vl1)
    ep.fill_pt_eta_phi_mass('l2_fit', vl2)
    ep.V_pt_fit[0] = vl12.Pt()
    ep.V_eta_fit[0] = vl12.Eta()
    ep.V_phi_fit[0] = vl12.Phi()
    ep.V_mass_fit[0] = vl12.M()

    # llbb system
    vllbb = v12 + vl12
    ep.llbb_fit_px[0] = vllbb.X()
    ep.llbb_fit_py[0] = vllbb.Y()
    ep.llbb_fit_pt[0] = vllbb.Pt()
    vr = fitter.get4Vec(4)
    vllbbr = vllbb + vr
    ep.llbbr_fit_px[0] = vllbbr.X()
    ep.llbbr_fit_py[0] = vllbbr.Y()
    ep.llbbr_fit_pt[0] = vllbbr.Pt()

    # HV variables
    ep.HVdPhi_fit[0] = abs(deltaPhi(v12.Phi(), vl12.Phi()))
    ep.HVdR_fit[0]   = deltaR(v12.Eta(), v12.Phi(), vl12.Eta(), vl12.Phi())
    ep.HVdEta_fit[0] = abs(v12.Eta() - vl12.Eta())
    ep.jjVPtRatio_fit[0] = v12.Pt() / ep.V_pt

    if (ep.kinfit_fit[0] != 1
        or ep.H_mass_fit[0] < 0
    ):
        ep.apply_fallback()
    else:
        ep.H_mass_fit_fallback[0] = ep.H_mass_fit[0]
        ep.H_pt_fit_fallback[0] = ep.H_pt_fit[0]
        ep.HVdPhi_fit_fallback[0] = ep.HVdPhi_fit[0]
        ep.jjVPtRatio_fit_fallback[0] = ep.jjVPtRatio_fit[0]
        ep.hJets_pt_0_fit_fallback[0] = v1.Pt()
        ep.hJets_pt_1_fit_fallback[0] = v2.Pt()

    # generator
    if (ep.GetBranch('nGenPart')
        and ep.GenLepIndex1 > -1
        and ep.GenLepIndex2 > -1
    ):
        # recalculate, since H jets are not always matched
        gen_H = mk_lv(ep.GenBJJ_pt, ep.GenBJJ_eta, ep.GenBJJ_phi, ep.GenBJJ_mass,)
        gen_l1 = mk_lv(
            ep.GenPart_pt[ep.GenLepIndex1],
            ep.GenPart_eta[ep.GenLepIndex1],
            ep.GenPart_phi[ep.GenLepIndex1],
            ep.GenPart_mass[ep.GenLepIndex1],
        )
        gen_l2 = mk_lv(
            ep.GenPart_pt[ep.GenLepIndex2],
            ep.GenPart_eta[ep.GenLepIndex2],
            ep.GenPart_phi[ep.GenLepIndex2],
            ep.GenPart_mass[ep.GenLepIndex2],
        )
        if (deltaR(gen_l1.Eta(), gen_l1.Phi(), vl1.Eta(), vl1.Phi())
          > deltaR(gen_l1.Eta(), gen_l1.Phi(), vl2.Eta(), vl2.Phi())):
            gen_l1, gen_l2 = gen_l2, gen_l1
        gen_Z = gen_l1 + gen_l2
        ep.fill_pt_eta_phi_mass('l1_gen', gen_l1)
        ep.fill_pt_eta_phi_mass('l2_gen', gen_l2)
        ep.fill_pt_eta_phi_mass('l12_gen', gen_Z)

        gen_ZH = gen_H + gen_Z
        ep.llbb_gen_px[0] = gen_ZH.X()
        ep.llbb_gen_py[0] = gen_ZH.Y()
        ep.ll_gen_mass[0] = gen_Z.M()

    return fitter  # needed for development


def prepare_io(input_file, output_file):
    if isinstance(input_file, ROOT.TChain):
        t = input_file
        of = ROOT.TFile(output_file, 'RECREATE')
    else:
        f = ROOT.TFile(input_file)
        assert not f.IsZombie(), 'input file is zombie, not running kinfit: ' + input_file
        f_other_objects = list((key.GetName(), f.Get(key.GetName()))
                         for key in f.GetListOfKeys()
                         if key.GetName() != 'Events')
        t = f.Get('Events')

        # mk new file and copy all other objects
        of = ROOT.TFile(output_file, 'RECREATE')
        for name, obj in f_other_objects:
            if isinstance(obj, ROOT.TTree):
                obj.CloneTree(-1, 'fast')
            else:
                obj.Clone(name)
            obj.Write(name)

    return f, t, of


def apply_kinfit(input_file, output_file, event_proxies=None):

    # jet resolutions from jme group
    # https://github.com/cms-jet/JRDatabase
    # https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution

    # TODO add 2017 resolutions
    ak4pfchs_ptres = ROOT.JME.JetResolution('aux/Spring16_25nsV6_DATA_PtResolution_AK4PFchs.txt')

    # prepare new tree and new variables
    f, t, of = prepare_io(input_file, output_file)
    ot = t.CloneTree(0)
    eps = event_proxies or [EventProxy()]
    for ep in eps:
        ep.init_output(ot, TREE_VARS)
    n_events = 0
    n_fitted = 0

    print 'starting kinfit loop'
    start_time = time.ctime()
    for e in t:
        n_events += 1

        if n_events % 1000 == 1:
            print 'at event %i' % n_events

        for ep in eps:
            ep.set_event(e)
        if (
            e.twoResolvedJets
            and (e.isZmm or e.isZee)
        ):
            for ep in eps:
                apply_fit_to_event(ep, ak4pfchs_ptres)
            ot.Fill()
            for ep in eps:
                ep.set_vals_to_zero()  # zero values after filling the tree
            n_fitted += 1
        else:
            for ep in eps:
                ep.apply_fallback()
            # don't run kinfit, just write
            ot.Fill()

    ot.Write()
    of.Close()
    del f
    print 'started kinfit loop at', start_time
    print 'finished kinfit loop at', time.ctime()
    print 'total number of events processed:    ', n_events
    print 'number of events with kinfit applied:', n_fitted


ep_mu_up    = EventProxy('mu_up',   'Muon_pt_corrected', lambda e: list(x*1.02 for x in e.Muon_pt_corrected))
ep_mu_down  = EventProxy('mu_down', 'Muon_pt_corrected', lambda e: list(x*0.98 for x in e.Muon_pt_corrected))
ep_el_up    = EventProxy('el_up',   'Electron_pt',       lambda e: list(x*1.02 for x in e.Electron_pt))
ep_el_down  = EventProxy('el_down', 'Electron_pt',       lambda e: list(x*0.98 for x in e.Electron_pt))
