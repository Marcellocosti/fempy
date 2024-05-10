import argparse
import os
import yaml
from rich import print

from ROOT import TFile, TCanvas, gInterpreter, TF1, TLegend, TLine, TH1, TGraph, TGraphErrors, TGraphAsymmErrors, TH1D, gStyle, TLatex

import fempy
from fempy import logger as log
from fempy.utils.format import TranslateToLatex
from fempy.utils.io import Load
from fempy.utils import style

gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/FitFunctions.cxx"')

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfg')
args = parser.parse_args()

# Load configuration file
with open(args.cfg, "r") as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError:
        log.critical('Yaml file not loaded')

style.SetStyle()

for plot in cfg:
    plot = plot["plot"]

    # Define the canvas
    nPanelsX, nPanelsY = fempy.utils.GetNPanels(1)
    cPlot = TCanvas("cPlot", "cPlot", 600*nPanelsX, 600*nPanelsY)
    cPlot.cd()

    inFilePars = TFile(plot['histoparfile'])
    hFilePars = Load(inFilePars, plot['histoparpath'])

    # Load the objects to draw
    funcs = []
    for func in plot['comps']:
        funcs.append(TF1(func['name'], func['name'], plot['drawrange'][0], plot['drawrange'][1], 7))
        funcs[-1].Draw('same')

    cPlot.Modified()
    cPlot.Update()

    # save canvas
    for ext in plot["opt"]["ext"]:
        cPlot.SaveAs(f'{os.path.splitext(plot["output"])[0]}.{ext}')
