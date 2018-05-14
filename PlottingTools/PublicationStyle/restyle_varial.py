#!/usr/bin/env python
import os
import sys

import ROOT

from restyle import restyle_varial, utils


ROOT.gROOT.SetBatch(True)

# Varial input file name.
VARIAL_FILE = '_varial_rootobjects.root.rt'

# Width and height of the restyled figure in pixels.
WIDTH = 800
HEIGHT = 800

# Basic annotation settings.
CMS_POSITION = 'left'
LUMI_TEXT = '35.9 fb^{-1} (13 TeV)'
EXTRA_TEXT = 'Preliminary'

# Output file type settings.
EXTS = ['.pdf']


def main():
    """Restyle AnalysisTools figures produced using Varial.

    The positional command line arguments are the Varial output
    directory and the output directory for the restyled plots.
    """
    src = os.path.abspath(sys.argv[1])
    dst = os.path.abspath(sys.argv[2])
    for dirpath, dirnames, filenames in os.walk(os.path.join(src, 'Plots')):
        if VARIAL_FILE in filenames:
            outdir = os.path.join(dst, os.path.basename(dirpath))
            utils.safe_makedirs(outdir)
            f = ROOT.TFile.Open(os.path.join(dirpath, VARIAL_FILE))
            dirs = [key.ReadObj() for key in f.GetListOfKeys()]
            canvases = [d.GetListOfKeys()[0].ReadObj() for d in dirs]
            for canvas in canvases:
                restyle_varial(canvas, WIDTH, HEIGHT, LUMI_TEXT, CMS_POSITION, EXTRA_TEXT, outdir, EXTS)
            f.Close()


if __name__ == '__main__':

    status = main()
    sys.exit(status)

