#!/usr/bin/env python
import os
import sys

import ROOT
import click

from utils import restyle_varial, safe_makedirs


# Varial input file name.
VARIAL_FILE = '_varial_rootobjects.root.rt'

# Width and height of the restyled figure in pixels.
WIDTH = 800
HEIGHT = 800

# Basic annotation settings.
CMS_POSITION = 'left'
LUMI_TEXT = '35.9 fb^{-1} (13 TeV)'
EXTRA_TEXT = 'Preliminary 2016'

# Output file type settings.
EXTS = ['.C', '.png', '.pdf']


@click.command()
@click.argument('src')
@click.argument('dst')
@click.option('--logy', is_flag=True, help='Set the y-axis scale of the restyled plots to logarithmic.')
def main(src, dst, logy):
    """Restyle AnalysisTools figures produced using Varial.

    The SRC argument should be the path to a Varial output directory,
    while DST is the path for the new output directory containing the
    restyled plots.
    """
    ROOT.gROOT.SetBatch(True)
    for dirpath, dirnames, filenames in os.walk(os.path.join(src, 'Plots')):
        if VARIAL_FILE in filenames:
            outdir = os.path.join(dst, os.path.basename(dirpath))
            safe_makedirs(outdir)
            f = ROOT.TFile.Open(os.path.join(dirpath, VARIAL_FILE))
            dirs = [key.ReadObj() for key in f.GetListOfKeys()]
            canvases = [d.GetListOfKeys()[0].ReadObj() for d in dirs]
            for canvas in canvases:
                restyle_varial(canvas, logy, WIDTH, HEIGHT, LUMI_TEXT, CMS_POSITION, EXTRA_TEXT, outdir, EXTS)
                canvas.IsA().Destructor(canvas)
            f.Close()


if __name__ == '__main__':

    main()

