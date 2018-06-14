import collections
import uuid
import sys

import ROOT
import numpy

import cms_figure


class Plotter(object):

    GROUPS = {
        'ZH': ['ZH125_ZLL_powheg'],
        'ggZH': ['ggZH125_ZLL_powheg'],
        'ZJets': [
            'DYToLL_madgraph',
            'DYToLL_HT100to200_madgraph',
            'DYToLL_HT200to400_madgraph',
            'DYToLL_HT400to600_madgraph',
            'DYToLL_HT600to800_madgraph',
            'DYToLL_HT800to1200_madgraph',
            'DYToLL_HT1200to2500_madgraph',
            'DYToLL_HT2500toInf_madgraph',
            'DYToLL_M4to50_HT70to100_madgraph',
            'DYToLL_M4to50_HT100to200_madgraph',
            'DYToLL_M4to50_HT200to400_madgraph',
            'DYToLL_M4to50_HT400to600_madgraph',
            'DYToLL_M4to50_HT600toInf_madgraph',
        ],
        'WJets': [
            'WJets_madgraph',
            'WJets-HT100To200',
            'WJets-HT200To400',
            'WJets-HT400To600',
            'WJets-HT600To800',
            'WJets-HT800To1200',
            'WJets-HT1200To2500',
        ],
        'TT': [
            'TT_AllHadronic',
            'TT_DiLep',
            'TT_SingleLep',
        ],
        'VV': [
            'WW',
            'WZ',
            'ZZ',
        ],
        'ST': [
            'ST_s-c_4f_lep_PSw',
            'ST_t-c_top_4f_inc',
            'ST_tW_antitop_5f_inc',
            'ST_tW_top_5f_inc_PSw',
        ],
        'QCD': [
            'QCD_HT300to500',
            'QCD_HT500to700',
            'QCD_HT700to1000',
            'QCD_HT1000to1500',
            'QCD_HT1500to2000',
            'QCD_HT2000toInf',
        ],
    }

    def __init__(self, input_pattern):
        for group, input_tokens in self.GROUPS.iteritems():
            setattr(self, group, ROOT.TChain('train'))
            chain = getattr(self, group)
            for input_token in input_tokens:
                chain.Add(input_pattern % input_token)

    def get_histogram(self, sample, bin_edges, selection=''):
        bin_edges = numpy.array(bin_edges, dtype=numpy.float32)
        n_bins = bin_edges.shape[0] - 1
        name = uuid.uuid4().hex
        h = ROOT.TH1F(name, sample, n_bins, bin_edges)
        getattr(self, sample).Project(name, 'CMS_vhbb_BDTG_Zll_LowPT_13TeV', selection)
        h.SetDirectory(0)
        return h

    def get_stack(self, bin_edges):
        histograms = collections.OrderedDict()

        # Because these input trees contain branches following the rebinner's convention, a cheap
        # way to plot signal region events is to select events where the bin index is above -1.
        histograms['QCD'] = self.get_histogram('QCD', bin_edges, 'weight * (bin_index_Zll_lowPt>-1)')
        histograms['Single top'] = self.get_histogram('ST', bin_edges, 'weight * (bin_index_Zll_lowPt>-1)')
        histograms['VV'] = self.get_histogram('VV', bin_edges, 'weight * (bin_index_Zll_lowPt>-1)')
        histograms['TT'] = self.get_histogram('TT', bin_edges, 'weight * (bin_index_Zll_lowPt>-1)')
        histograms['W+Jets'] = self.get_histogram('WJets', bin_edges, 'weight * (bin_index_Zll_lowPt>-1)')
        histograms['Z+Jets'] = self.get_histogram('ZJets', bin_edges, 'weight * (bin_index_Zll_lowPt>-1)')
        histograms['ggZH(b#bar{b})'] = self.get_histogram('ggZH', bin_edges, 'weight * (bin_index_Zll_lowPt>-1)')
        histograms['ZH(b#bar{b})'] = self.get_histogram('ZH', bin_edges, 'weight * (bin_index_Zll_lowPt>-1)')

        histograms['QCD'].SetFillColor(613)
        histograms['Single top'].SetFillColor(70)
        histograms['VV'].SetFillColor(17)
        histograms['TT'].SetFillColor(4)
        histograms['W+Jets'].SetFillColor(81)
        histograms['Z+Jets'].SetFillColor(5)
        histograms['ggZH(b#bar{b})'].SetFillColor(625)
        histograms['ZH(b#bar{b})'].SetFillColor(2)

        stack = ROOT.THStack('stack', '')
        for histogram in histograms.values():
            histogram.SetLineColor(ROOT.kBlack)
            stack.Add(histogram)

        legend1 = ROOT.TLegend(0.51, 0.64, 0.73, 0.92)
        legend1.SetFillColor(0)
        legend1.SetLineColor(0)
        legend1.SetShadowColor(0)
        legend1.SetTextFont(62)
        legend1.SetTextSize(0.03)
        legend1.SetBorderSize(1)

        legend2 = ROOT.TLegend(0.74, 0.64, 0.96, 0.92)
        legend2.SetFillColor(0)
        legend2.SetLineColor(0)
        legend2.SetShadowColor(0)
        legend2.SetTextFont(62)
        legend2.SetTextSize(0.03)
        legend2.SetBorderSize(1)

        for i, name in enumerate(reversed(histograms)):
            if i < 4:
                legend1.AddEntry(histograms[name], name, 'f')
            else:
                legend2.AddEntry(histograms[name], name, 'f')

        ROOT.SetOwnership(stack, False)

        return stack, legend1, legend2


def get_stack_info(stack):
    n_signal = []
    n_background = []
    s_over_b = []

    n_bins = stack.GetHists()[0].GetNbinsX()
    for i in xrange(1, n_bins + 1):
        s, b = 0, 0
        for hist in stack.GetHists():
            if hist.GetTitle() in {'ZH', 'ggZH'}:
                s += hist.GetBinContent(i)
            else:
                b += hist.GetBinContent(i)
        n_signal.append(s)
        n_background.append(b)
        s_over_b.append(0 if b == 0 else s / b)

    return n_signal, n_background, s_over_b


def main():
    ROOT.gROOT.SetBatch(True)

    plotter = Plotter('../SkimTreeBuilder/samples_2017/%s.root')

    best_bin_edges = {
        'asimov': [-1, -0.23956692218780518, 0.23340271413326263, 0.5274667739868164, 0.6607716977596283, 0.7466760575771332, 0.827688992023468, 0.8841321468353271, 1],
        'poisson': [-1, -0.23956692218780518, 0.23340271413326263, 0.5274667739868164, 0.6607716977596283, 0.7466760575771332, 0.827688992023468, 0.8841321468353271, 1],
    }

    best_significance = {
        'asimov': 1.397680535286613,
        'poisson': 1.4166875291979029,
    }

    text = ROOT.TLatex()
    text.SetTextAlign(12)
    text.SetTextFont(42)

    for metric, bin_edges in best_bin_edges.iteritems():
        stack, legend1, legend2 = plotter.get_stack(bin_edges)
        with cms_figure.TDRStyle() as style:
            style.SetTitleOffset(1.2, 'X')
            style.SetTitleOffset(1.5, 'Y')
            style.SetTitleSize(0.04, 'XY')
            style.SetLabelSize(0.04, 'XY')
            canvas = ROOT.TCanvas(uuid.uuid4().hex, '')
            canvas.SetLogy(True)
            stack.SetMinimum(10**-1)
            stack.SetMaximum(10**8)
            stack.Draw('hist')
            stack.GetXaxis().SetTitle('BDT output')
            stack.GetYaxis().SetTitle('Entries')
            legend1.Draw()
            legend2.Draw()
            cms_figure.draw_labels('41.53 fb^{-1} (13 TeV)', extra_text='Simulation Preliminary')
            text.SetTextSize(0.025)
            text.DrawLatexNDC(0.2, 0.79, '{0} Significance Score'.format(metric.title()))
            text.DrawLatexNDC(0.2, 0.74, '{0:.2f}'.format(best_significance[metric]))
            canvas.SaveAs('bdt_2lep_lowPt_rebinned_{0}.png'.format(metric))

        n_signal, n_background, s_over_b = get_stack_info(stack)

        h_s_over_b = stack.GetHists().Last().Clone('h_s_over_b')
        h_s_over_b.SetFillColor(0)
        h_s_over_b.SetLineColor(ROOT.kRed)
        h_s_over_b.SetLineWidth(2)
        for i, value in enumerate(s_over_b, start=1):
            h_s_over_b.SetBinContent(i, value)

        with cms_figure.TDRStyle() as style:
            style.SetTitleOffset(1.2, 'X')
            style.SetTitleOffset(1.5, 'Y')
            style.SetTitleSize(0.04, 'XY')
            style.SetLabelSize(0.04, 'XY')
            canvas = ROOT.TCanvas(uuid.uuid4().hex, '')
            h_s_over_b.SetMinimum(0.)
            h_s_over_b.SetMaximum(h_s_over_b.GetMaximum() * 2)
            h_s_over_b.Draw('hist')
            h_s_over_b.GetXaxis().SetTitle('BDT output')
            h_s_over_b.GetYaxis().SetTitle('S / B')
            cms_figure.draw_labels('41.53 fb^{-1} (13 TeV)', extra_text='Simulation Preliminary')
            text.DrawLatexNDC(0.2, 0.79, 'Signal Yields')
            text.SetTextSize(0.025)
            for i, value in enumerate(n_signal, start=1):
                text.DrawLatexNDC(0.2, 0.79 - i * 0.04, 'Bin #{0!s}: {1:.2f}'.format(i, value))
            canvas.SaveAs('s_over_b_2lep_lowPt_rebinned_{0}.png'.format(metric))


if __name__ == '__main__':

    status = main()
    sys.exit(status)

