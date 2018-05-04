import os
import uuid

import ROOT
import cms_figure


__all__ = [
    'restyle_varial',
]


ROOT.TGaxis.SetMaxDigits(3)


def _restyle_upper_pad(pad):
    # Resize the pad.
    pad.SetPad(0.0, 0.301, 1.0, 1.0)
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.05)
    pad.SetBottomMargin(0.018)
    pad.SetTopMargin(0.08)

    # Then modify...
    primitives = pad.GetListOfPrimitives()
    # the stacked plot,
    stack = primitives[1]
    hists = []
    titles = []
    colors = []
    for hist in stack.GetHists():
        hists.append(hist)
        titles.append(hist.GetTitle())
        colors.append(hist.GetFillColor())
        hist.SetLineColor(1)
    y_max_stack = stack.GetMaximum()
    x_axis, y_axis = stack.GetXaxis(), stack.GetYaxis()
    x_axis.SetNdivisions(510)
    x_axis.SetLabelSize(0)
    x_axis.SetTitleSize(0)
    y_axis.CenterTitle(False)
    y_axis.SetNdivisions(510)
    y_axis.SetLabelSize(0.05)
    y_axis.SetLabelOffset(0.007)
    bin_width = x_axis.GetBinWidth(1)
    y_axis.SetTitle('Events / {0:.1f}'.format(bin_width))
    y_axis.SetTitleSize(0.05)
    y_axis.SetTitleOffset(1.4)
    # the shaded MC uncertainty (converted into a TGraphAsymmErrors
    # because SetErrorX(0) wipes out the error bands globally),
    mc_unc = primitives[2]
    mc_unc_new = ROOT.TGraphAsymmErrors(mc_unc)
    for i in xrange(mc_unc_new.GetN()):
        mc_unc_new.SetPointEXlow(i, bin_width / 2.)
        mc_unc_new.SetPointEXhigh(i, bin_width / 2.)
    mc_unc_new.SetFillColor(923)
    mc_unc.Reset()
    primitives.Add(mc_unc_new, '2SAME')
    # and the data points.
    data = primitives[3]
    y_max_data = data.GetMaximum()
    data_new = data.Clone()
    data_new.SetMarkerSize(0.9)
    data.Reset()
    primitives.Add(data_new, 'E1SAME')
    # Adjust the y-axis for better viewing.
    if y_max_data > y_max_stack:
        stack.SetMaximum(1.7 * y_max_data)
    else:
        stack.SetMaximum(1.7 * y_max_stack)

    # Finally, add two legends...
    legend1 = ROOT.TLegend(0.52, 0.56, 0.77, 0.88, '', 'nbNDC')
    legend1.SetFillColor(0)
    legend1.SetTextSize(0.035)
    legend2 = ROOT.TLegend(0.69, 0.56, 0.94, 0.88, '', 'nbNDC')
    legend2.SetFillColor(0)
    legend2.SetTextSize(0.035)
    # and fill them with entries.
    legend1.AddEntry(data_new, 'Data', 'ep')
    stack_legend_info = zip(reversed(hists), reversed(titles), reversed(colors))
    max_legend_entries = len(hists) + 2
    for i, (hist, title, color) in enumerate(stack_legend_info, start=2):
        if i <= max_legend_entries / 2:
            entry = legend1.AddEntry(hist, title, 'f')
        else:
            entry = legend2.AddEntry(hist, title, 'f')
        entry.SetFillColor(color)
    legend2.AddEntry(mc_unc, 'MC unc. (stat.)', 'f')
    primitives.Add(legend1)
    primitives.Add(legend2)


def _restyle_lower_pad(pad):
    # Resize and modify the pad.
    pad.SetPad(0.0, 0.0, 1.0, 0.299)
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.05)
    pad.SetBottomMargin(0.35)
    pad.SetGridy(0)

    # Then modify...
    primitives = pad.GetListOfPrimitives()
    # the shaded ratio uncertainty (converted into a TGraphAsymmErrors
    # because SetErrorX(0) wipes out the error bands globally),
    ratio_unc = primitives[1]
    bin_width = ratio_unc.GetXaxis().GetBinWidth(1)
    ratio_unc_new = ROOT.TGraphAsymmErrors(ratio_unc)
    for i in xrange(ratio_unc_new.GetN()):
        ratio_unc_new.SetPointEXlow(i, bin_width / 2.)
        ratio_unc_new.SetPointEXhigh(i, bin_width / 2.)
    ratio_unc_new.SetFillColor(923)
    ratio_unc.Reset()
    primitives.Add(ratio_unc_new, '2SAME')
    ratio_unc.SetMaximum(0.999)
    ratio_unc.SetMinimum(-0.999)
    x_axis, y_axis = ratio_unc.GetXaxis(), ratio_unc.GetYaxis()
    x_axis.SetNdivisions(510)
    x_axis.SetLabelOffset(0.007)
    x_axis.SetLabelSize(0.11)
    x_axis.SetTickLength()
    x_axis.SetTitleOffset(1.2)
    x_axis.SetTitleSize(0.11)
    y_axis.SetNdivisions(505)
    y_axis.SetLabelOffset(0.007)
    y_axis.SetLabelSize(0.11)
    y_axis.SetTitleOffset(0.55)
    y_axis.SetTitleSize(0.11)
    # and the ratio points.
    ratio = primitives[2]
    ratio_new = ratio.Clone()
    ratio_new.SetMarkerSize(0.9)
    ratio.Reset()
    primitives.Add(ratio_new, 'E1SAME')

    # Finally, add the legend.
    legend = ROOT.TLegend(0.32, 0.86, 0.93, 0.97, '', 'nbNDC')
    legend.SetTextSize(0.075)
    legend.SetNColumns(2)
    legend.AddEntry(ratio_unc, 'MC unc. (stat.)', 'f')
    primitives.Add(legend)


def restyle_varial(canvas, width=800, height=800, lumi_text='', cms_position='left', extra_text='', dst=None, exts=None):
    """Restyle an AnalysisTools figure produced by Varial for publication.

    The current style conforms to the Publications Committees standards and
    the aesthetic choices made for the figures published in HIG-16-044.

    Parameters
    ----------
    canvas : ROOT.TCanvas
        A TCanvas produced by Varial, typically stored in a file called
        "_varial_rootobjects.root.rt".
    width : int, optional
        The width of the restyled figure in pixels. The default is 800.
    height : int, optional
        The height of the restyled figure in pixels. The default is 800.
    lumi_text : string, optional
        The luminosity label text. Data taking periods must be separated by
        the "+" symbol, e.g. "19.7 fb^{-1} (8 TeV) + 4.9 fb^{-1} (7 TeV)".
        The default is no luminosity label.
    cms_position : string, optional
        The CMS label position on the restyled figure:
        * left : The top left corner inside the frame (default)
        * center : The top center inside the frame
        * right : The top right corner inside the frame
        * outside : The top left corner outside the frame
    extra_text : string, optional
        The sublabel text positioned below the CMS label inside of the frame
        or to the right of the CMS label outside of the frame. Common examples
        are "Preliminary" or "Unpublished". The default is no sublabel.
    dst : path, optional
        The output directory where the restyled figure is saved. The default
        is the current working directory. This parameter is ignored when
        `exts` is None.
    exts : string or iterable of strings, optional
        The file extensions which specify the output file format, e.g. '.pdf',
        '.png', '.root', or '.C'. The default is draw the figure directly in
    """
    if dst is None:
        dst = os.getcwd()
    name = canvas.GetName().split('_')[0]
    primitives = canvas.GetListOfPrimitives()
    upper_pad = primitives[0]
    lower_pad = primitives[1]
    with cms_figure.TDRStyle() as style:
        style.SetErrorX(0)
        # Set up the restyled canvas.
        new_name = 'restyled_canvas_{0!s}'.format(uuid.uuid4().hex)
        new_canvas = ROOT.TCanvas(new_name, '', width, height)
        new_canvas.SetLeftMargin(0.12)
        new_canvas.SetRightMargin(0.04)
        new_canvas.SetBottomMargin(0.12)
        new_canvas.SetTopMargin(0.08)
        # Restyle and add the upper pad.
        _restyle_upper_pad(upper_pad)
        upper_pad.Draw()
        upper_pad.cd()
        cms_figure.draw_labels(lumi_text, cms_position, extra_text)
        upper_pad.Modified()
        upper_pad.Update()
        upper_pad.RedrawAxis()
        # Restyle and add the lower pad.
        new_canvas.cd()
        _restyle_lower_pad(lower_pad)
        lower_pad.Draw()
        lower_pad.Modified()
        lower_pad.Update()
        lower_pad.RedrawAxis()
        # Draw or save the restyled canvas.
        if exts is None:
            new_canvas.Draw()
        else:
            if isinstance(exts, basestring):
                exts = [exts]
            for ext in exts:
                new_canvas.SaveAs(os.path.join(dst, name + ext))

