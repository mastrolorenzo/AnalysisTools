from math import pi
import itertools
import numpy
import time
import ROOT
import os


# constants
llbb_pxy_var = 8**2  # with no additional jets above pt > 20 GeV and no eta cut
hh4b_res_scale = 0.62
z_width = 1.7 * 3


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


# jet resolutions from jme group
# https://github.com/cms-jet/JRDatabase
# https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
ak4pfchs_ptres = None


# auxillary functions
def deltaPhi(phi1, phi2):
    dphi = (phi1 - phi2 + pi) % (2*pi) - pi
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
        e.Jet_bReg[ind] if regressed else e.Jet_Pt[ind],
        e.Jet_eta[ind],
        e.Jet_phi[ind],
        e.Jet_mass[ind],
    )

def mk_jet_cart_cov(v, rho):
    jme_par = ROOT.JME.JetParameters().setRho(rho).setJetPt(v.Pt()).setJetEta(v.Eta())
    jme_res = ak4pfchs_ptres.getResolution(jme_par) * v.Pt()
    cos_phi = ROOT.TMath.Cos(v.Phi())
    sin_phi = ROOT.TMath.Sin(v.Phi())
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
    cov[0][0] = hh4b_res * res_scale**2
    cov[1][1] = ErrEta_Signal(v.Et())
    cov[2][2] = ErrPhi_Signal(v.Et())
    return v, hh4b_res, ROOT.TFitParticleEtEtaPhi(v, cov)


# ... and from leptons as well
def mk_mu_lv_ind(e, ind):
    return mk_lv(
        e.Muon_pt[ind],
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


class VarContainer(object):
    def __init__(self, output_tree):
        self.tree_vars = [

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

            'HVdPhi_fit',    'HVdR_fit',       'HVdEta_fit',

            'hj1_reg_res',
            'hj2_reg_res',
            'hj1_hh4b_res',
            'hj2_hh4b_res',
            'hj1_jme_res',
            'hj2_jme_res',

            # leptons
            'l1_pt',        'l1_eta',        'l1_phi',        'l1_mass',
            'l2_pt',        'l2_eta',        'l2_phi',        'l2_mass',
            'V_eta',
            'l1_gen_pt',    'l1_gen_eta',    'l1_gen_phi',    'l1_gen_mass',
            'l2_gen_pt',    'l2_gen_eta',    'l2_gen_phi',    'l2_gen_mass',
            'l12_gen_pt',   'l12_gen_eta',   'l12_gen_phi',   'l12_gen_mass',
            'l1_fit_pt',    'l1_fit_eta',    'l1_fit_phi',    'l1_fit_mass',
            'l2_fit_pt',    'l2_fit_eta',    'l2_fit_phi',    'l2_fit_mass',
            'V_pt_fit',      'V_eta_fit',      'V_phi_fit',      'V_mass_fit',

            'l1_res',
            'l2_res',

            # llbb system / recoil jet
            # 'recoil_pt',     'recoil_eta',     'recoil_phi',     'recoil_mass',
            # 'recoil_fit_pt', 'recoil_fit_eta', 'recoil_fit_phi', 'recoil_fit_mass',

            # 'llbb_px',       'llbb_py',
            'llbb_fit_px',   'llbb_fit_py',
            # 'llbbr_fit_px',  'llbbr_fit_py',
            'llbb_gen_px',   'llbb_gen_py',
            'll_gen_mass',

            # aux
            'kinfit_fit',
            'kinfit_getNDF',
            'kinfit_getS',
            'kinfit_getF',

            'H_mass_fit_fallback',
        ]

        for var in self.tree_vars:
            setattr(self, var, numpy.zeros(1,float))
            output_tree.Branch(var, getattr(self, var), var+'/D')

    def set_vals_to_zero(self):
        for var in self.tree_vars:
            getattr(self, var)[0] = 0.

    def fill_pt_eta_phi_mass(self, token, vec):
        getattr(self, token+'_pt'  )[0] = vec.Pt()
        getattr(self, token+'_eta' )[0] = vec.Eta()
        getattr(self, token+'_phi' )[0] = vec.Phi()
        getattr(self, token+'_mass')[0] = vec.M()


def apply_fit_to_event(e, c):

    # normal jets
    v1 = mk_lv_ind(e, e.hJetInd1) #e.hJCMVAV2idx[0])
    v2 = mk_lv_ind(e, e.hJetInd2) #e.hJCMVAV2idx[1])
    v12 = v1 + v2

    c.fill_pt_eta_phi_mass('hj1', v1)
    c.fill_pt_eta_phi_mass('hj2', v2)
    c.fill_pt_eta_phi_mass('hj12', v12)


    # regressed jets and fit particles
    v1, v1_pt_var, p1 = mk_fit_particle(e, e.hJetInd1, hh4b_res_scale) #e.hJCMVAV2idx[0])
    v2, v2_pt_var, p2 = mk_fit_particle(e, e.hJetInd2, hh4b_res_scale) #e.hJCMVAV2idx[1])
    v12 = v1 + v2
    c.fill_pt_eta_phi_mass('hj1_reg', v1)
    c.fill_pt_eta_phi_mass('hj2_reg', v2)
    c.fill_pt_eta_phi_mass('hj12_reg', v12)

    c.n_hj_matched[0] = 0
    if deltaR(c.hj1_reg_eta[0], c.hj1_reg_phi[0], e.GenBJ1_eta, e.GenBJ1_phi) < 0.2:
        v1_gen = mk_lv(e.GenBJ1_pt,e.GenBJ1_eta,e.GenBJ1_phi,e.GenBJ1_mass)
        c.fill_pt_eta_phi_mass('hj1_gen', v1_gen)
        c.n_hj_matched[0] += 1
    elif deltaR(c.hj1_reg_eta[0], c.hj1_reg_phi[0], e.GenBJ2_eta, e.GenBJ2_phi) < 0.2:
        v1_gen = mk_lv(e.GenBJ2_pt,e.GenBJ2_eta,e.GenBJ2_phi,e.GenBJ2_mass)
        c.fill_pt_eta_phi_mass('hj1_gen', v1_gen)
        c.n_hj_matched[0] += 1

    if deltaR(c.hj2_reg_eta[0], c.hj2_reg_phi[0], e.GenBJ1_eta, e.GenBJ1_phi) < 0.2:
        v2_gen = mk_lv(e.GenBJ1_pt,e.GenBJ1_eta,e.GenBJ1_phi,e.GenBJ1_mass)
        c.fill_pt_eta_phi_mass('hj2_gen', v2_gen)
        c.n_hj_matched[0] += 1
    elif deltaR(c.hj2_reg_eta[0], c.hj2_reg_phi[0], e.GenBJ2_eta, e.GenBJ2_phi) < 0.2:
        v2_gen = mk_lv(e.GenBJ2_pt,e.GenBJ2_eta,e.GenBJ2_phi,e.GenBJ2_mass)
        c.fill_pt_eta_phi_mass('hj2_gen', v2_gen)
        c.n_hj_matched[0] += 1

    if c.n_hj_matched[0] == 2:
        v12_gen = v1_gen + v2_gen
        c.fill_pt_eta_phi_mass('hj12_gen', v12_gen)

    jme_p1 = ROOT.JME.JetParameters().setRho(e.fixedGridRhoFastjetAll).setJetPt(c.hj1_reg_pt[0]).setJetEta(c.hj1_reg_eta[0])
    jme_p2 = ROOT.JME.JetParameters().setRho(e.fixedGridRhoFastjetAll).setJetPt(c.hj2_reg_pt[0]).setJetEta(c.hj2_reg_eta[0])
    c.hj1_jme_res[0] = ak4pfchs_ptres.getResolution(jme_p1) * c.hj1_reg_pt[0]
    c.hj2_jme_res[0] = ak4pfchs_ptres.getResolution(jme_p2) * c.hj2_reg_pt[0]

    c.hj1_hh4b_res[0] = v1_pt_var**.5
    c.hj2_hh4b_res[0] = v2_pt_var**.5

    # TODO hj1_reg_res, hj2_reg_res

    # base recoil covariance
    cov_recoil = ROOT.TMatrixD(4, 4)
    cov_recoil.Zero()
    cov_recoil[0][0] = llbb_pxy_var
    cov_recoil[1][1] = llbb_pxy_var
    cov_recoil[2][2] = 1
    cov_recoil[3][3] = 1e12

    # recoil vector and covariance additions
    recoil_jet_idxs = (
        i
        for i in xrange(e.nJet)
        if (
            i not in (e.hJetInd1, e.hJetInd2)
            and e.Jet_puId[i] > 1
            and e.Jet_jetId[i] > 1
            and e.Jet_lepFilter[i] > 0
            and e.Jet_Pt[i] > 20
            # and abs(e.Jet_eta) < 2.5
        )
    )
    recoil_jet_vs = list(
        mk_lv_ind(e, i)
        for i in recoil_jet_idxs
    )
    v_recoil = sum(recoil_jet_vs, ROOT.TLorentzVector())
    cov_recoil = sum(
        (mk_jet_cart_cov(v, e.fixedGridRhoFastjetAll) for v in recoil_jet_vs),
        cov_recoil
    )
    p_recoil = ROOT.TFitParticleECart(v_recoil, cov_recoil)

    vl1, vl2 = get_sel_leptons(e)
    vl12 = vl1 + vl2
    _, vl1_pt_var, pl1 = mk_fit_particle_lepton(e, vl1, e.lepInd1)
    _, vl2_pt_var, pl2 = mk_fit_particle_lepton(e, vl2, e.lepInd2)

    c.fill_pt_eta_phi_mass('l1', vl1)
    c.fill_pt_eta_phi_mass('l2', vl2)
    c.V_eta[0] = vl12.Eta()

    c.l1_res[0] = vl1_pt_var**.5
    c.l2_res[0] = vl2_pt_var**.5

    # constraints
    cons_MZ = ROOT.TFitConstraintMGaus('cons_MZ', '', None, None, 91., z_width)
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
    c.kinfit_fit[0] = fitter.fit() + 1  # 0: not run; 1: fit converged; 2: fit didn't converge
    c.kinfit_getNDF[0] = fitter.getNDF()
    c.kinfit_getS[0] = fitter.getS()
    c.kinfit_getF[0] = fitter.getF()

    # higgs jet output vectors
    v1 = fitter.get4Vec(0)
    v2 = fitter.get4Vec(1)
    v12 = v1 + v2

    c.fill_pt_eta_phi_mass('hj1_fit', v1)
    c.fill_pt_eta_phi_mass('hj2_fit', v2)
    c.H_pt_fit[0] = v12.Pt()
    c.H_eta_fit[0] = v12.Eta()
    c.H_phi_fit[0] = v12.Phi()
    c.H_mass_fit[0] = v12.M()

    # lepton output vectors
    vl1 = fitter.get4Vec(2)
    vl2 = fitter.get4Vec(3)
    vl12 = vl1 + vl2

    c.fill_pt_eta_phi_mass('l1_fit', vl1)
    c.fill_pt_eta_phi_mass('l2_fit', vl2)
    c.V_pt_fit[0] = vl12.Pt()
    c.V_eta_fit[0] = vl12.Eta()
    c.V_phi_fit[0] = vl12.Phi()
    c.V_mass_fit[0] = vl12.M()

    # llbb system
    vllbb = v12 + vl12
    c.llbb_fit_px[0] = vllbb.X()
    c.llbb_fit_py[0] = vllbb.Y()
    # vr = fitter.get4Vec(4)
    # vllbbr = vllbb + vr
    # c.llbbr_fit_px[0] = vllbbr.X()
    # c.llbbr_fit_py[0] = vllbbr.Y()

    # HV variables
    c.HVdPhi_fit[0] = deltaPhi(v12.Phi(), vl12.Phi())
    c.HVdR_fit[0]   = deltaR(v12.Eta(), v12.Phi(), vl12.Eta(), vl12.Phi())
    c.HVdEta_fit[0] = abs(v12.Eta() - vl12.Eta())

    if (c.kinfit_fit[0] != 1
        or c.H_mass_fit[0] < 0
    ):
        c.H_mass_fit_fallback[0] = c.hj12_reg_mass[0]
    else:
        c.H_mass_fit_fallback[0] = c.H_mass_fit[0]

    # generator
    if (e.GetBranch('nGenPart')
        and e.GenLepIndex1 > -1
        and e.GenLepIndex2 > -1
    ):
        # recalculate, since H jets are not always matched
        gen_H = mk_lv(e.GenBJJ_pt, e.GenBJJ_eta, e.GenBJJ_phi, e.GenBJJ_mass,)
        gen_l1 = mk_lv(
            e.GenPart_pt[e.GenLepIndex1],
            e.GenPart_eta[e.GenLepIndex1],
            e.GenPart_phi[e.GenLepIndex1],
            e.GenPart_mass[e.GenLepIndex1],
        )
        gen_l2 = mk_lv(
            e.GenPart_pt[e.GenLepIndex2],
            e.GenPart_eta[e.GenLepIndex2],
            e.GenPart_phi[e.GenLepIndex2],
            e.GenPart_mass[e.GenLepIndex2],
        )
        if (deltaR(gen_l1.Eta(), gen_l1.Phi(), vl1.Eta(), vl1.Phi())
          > deltaR(gen_l1.Eta(), gen_l1.Phi(), vl2.Eta(), vl2.Phi())):
            gen_l1, gen_l2 = gen_l2, gen_l1
        gen_Z = gen_l1 + gen_l2
        c.fill_pt_eta_phi_mass('l1_gen', gen_l1)
        c.fill_pt_eta_phi_mass('l2_gen', gen_l2)
        c.fill_pt_eta_phi_mass('l12_gen', gen_Z)

        gen_ZH = gen_H + gen_Z
        c.llbb_gen_px[0] = gen_ZH.X()
        c.llbb_gen_py[0] = gen_ZH.Y()
        c.ll_gen_mass[0] = gen_Z.M()


def apply_kinfit(input_file, output_file, is_data, data_year):
    global ak4pfchs_ptres

    if is_data:
        ak4pfchs_ptres = ROOT.JME.JetResolution('aux/Spring16_25nsV6_DATA_PtResolution_AK4PFchs.txt')
    else:
        ak4pfchs_ptres = ROOT.JME.JetResolution('aux/Spring16_25nsV6_MC_PtResolution_AK4PFchs.txt')

    # TODO add 2017 resolutions
    # if data_year == 2016:
    # elif data_year == 2017:
    # else:
    #     raise RuntimeError('data_year is not 2016 or 2017.')

    if isinstance(input_file, ROOT.TChain):
        t = input_file
        of = ROOT.TFile(output_file, 'RECREATE')
    else:
        f = ROOT.TFile(input_file)
        if f.IsZombie():
            raise RuntimeError('input file is zombie, not running kinfit: ' + input_file)
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

    # prepare new tree and new variables
    ot = t.CloneTree(0)
    c = VarContainer(ot)
    n_events = 0
    n_fitted = 0

    print 'starting kinfit loop'
    start_time = time.ctime()
    for e in t:
        n_events += 1

        if n_events % 1000 == 1:
            print 'at event %i' % n_events

        if (
            e.V_pt >= 150
            and e.twoResolvedJets
            and (e.isZmm or e.isZee)
        ):
            apply_fit_to_event(e, c)
            ot.Fill()
            c.set_vals_to_zero()  # zero values after filling the tree
            n_fitted += 1
        else:
            # don't run kinfit, just write
            ot.Fill()

    ot.Write()
    of.Close()
    print 'started kinfit loop at', start_time
    print 'finished kinfit loop at', time.ctime()
    print 'total number of events processed:    ', n_events
    print 'number of events with kinfit applied:', n_fitted
