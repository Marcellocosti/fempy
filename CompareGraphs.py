import argparse
import os
import yaml
from rich import print

from ROOT import TFile, TCanvas, TLegend, TLine, TH1, TGraph, TGraphErrors, TGraphAsymmErrors, TH1D

import fempy
from fempy import logger as log
from fempy.utils.format import TranslateToLatex
from fempy.utils.io import Load
from fempy.utils import style

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

    panels = {'default': 1}
    if plot['ratio']['enable']:
        panels['ratio'] = len(panels)+1
    if plot['relunc']['enable']:
        panels['relunc'] = len(panels)+1

    # Load the objects to draw
    inObjs = []
    legends = []
    drawOpts = []
    for inputCfg in plot["input"]:
        inFile = TFile(inputCfg['file'])

        inObj = Load(inFile, inputCfg['name'])

        if isinstance(inObj, TH1):
            inObj.SetDirectory(0)
            inObj.Rebin(inputCfg['rebin'])

            if inputCfg['normalize']:
                inObj.Scale(1./inObj.Integral())
            if inputCfg['normalizecf']:
                inObj.Scale(inputCfg['normalizecf'])

        inObj.SetLineColor(style.GetColor(inputCfg['color']))
        inObj.SetMarkerColor(style.GetColor(inputCfg['color']))
        inObj.SetLineWidth(inputCfg.get('thickness', 1))
        drawOpts.append(inputCfg.get('drawopt', 'p' if isinstance(inObj, TH1) else 'pe'))
        inObj.SetMarkerStyle(inputCfg['markerstyle'])
        inObjs.append(inObj)
        legends.append(inputCfg['legend'])

    # Define the canvas
    nPanelsX, nPanelsY = fempy.utils.GetNPanels(len(panels))
    cPlot = TCanvas("cPlot", "cPlot", 600*nPanelsX, 600*nPanelsY)
    cPlot.Divide(nPanelsX, nPanelsY)
    pad = cPlot.cd(1)
    pad.SetLogx(plot["opt"]["logx"])
    pad.SetLogy(plot["opt"]["logy"])

    fx1 = plot['opt']['rangex'][0]
    fy1 = plot['opt']['rangey'][0]
    fx2 = plot['opt']['rangex'][1]
    fy2 = plot['opt']['rangey'][1]
    pad.DrawFrame(fx1, fy1, fx2, fy2, TranslateToLatex(plot['opt']['title']))

    legx1 = plot['opt']['leg']['posx'][0]
    legy1 = plot['opt']['leg']['posy'][0]
    legx2 = plot['opt']['leg']['posx'][1]
    legy2 = plot['opt']['leg']['posy'][1]
    leg = TLegend(legx1, legy1, legx2, legy2)

    for iObj, (inObj, legend) in enumerate(zip(inObjs, legends)):
        if isinstance(inObj, TGraph):
            inObj.Draw('same p')
        elif isinstance(inObj, TH1):
            inObj.Draw("same pe")

        # Compute statistics for hist in the displayed range
        if isinstance(inObj, TH1):
            firstBin = inObj.FindBin(plot['opt']['rangex'][0]*1.0001)
            lastBin = inObj.FindBin(plot['opt']['rangex'][1]*0.9999)
            inObj.GetXaxis().SetRange(firstBin, lastBin)
            print(f'{legend}: mean = {inObj.GetMean()} sigma = {inObj.GetStdDev()}')
            if plot['opt']['leg']['mean']:
                legend += f';  #mu={inObj.GetMean():.3f}'
            if plot['opt']['leg']['sigma']:
                legend += f';  #sigma={inObj.GetStdDev():.3f}'
        leg.AddEntry(inObj, legend, 'lp')
    for line in plot['opt']['lines']:
        x1 = plot['opt']['rangex'][0] if(line['coordinates'][0] == 'min') else line['coordinates'][0]
        y1 = plot['opt']['rangey'][0] if(line['coordinates'][1] == 'min') else line['coordinates'][1]
        x2 = plot['opt']['rangex'][1] if(line['coordinates'][2] == 'max') else line['coordinates'][2]
        y2 = plot['opt']['rangey'][1] if(line['coordinates'][3] == 'max') else line['coordinates'][3]
        inputline = TLine(x1, y1, x2, y2)
        inputline.SetLineColor(style.GetColor(line['color']))
        inputline.SetLineWidth(line['thickness'])
        inputline.Draw("same")
        leg.AddEntry(inputline, TranslateToLatex(line['legendtag']),"l")
        
    leg.SetHeader(TranslateToLatex(plot['opt']['leg']['header']), 'C')
    leg.Draw()

    # Compute ratio wrt the first obj
    if plot['ratio']['enable']:
        pad = cPlot.cd(panels['ratio'])
        pad.SetLogx(plot['ratio']['logx'])
        pad.SetLogy(plot['ratio']['logy'])
        x1 = plot['opt']['rangex'][0]
        y1 = plot['ratio']['rangey'][0]
        x2 = plot['opt']['rangex'][1]
        y2 = plot['ratio']['rangey'][1]
        frame = pad.DrawFrame(x1, y1, x2, y2, TranslateToLatex(plot['opt']['title']))
        frame.GetYaxis().SetTitle('Ratio')
        hDen = inObjs[0].Clone()
        hDen.Rebin(plot['ratio']['rebin'])
        hDen.Sumw2()

        if isinstance(inObj, TH1):
            for inObj in inObjs[1:]:
                hRatio = inObj.Clone()
                hRatio.Rebin(plot['ratio']['rebin'])
                hRatio.Divide(hDen)
                hRatio.Draw('same pe')
        else:
            log.error('Ratio for type %s is not implemented. Skipping this object', type(inObj))
            continue

        line = TLine(plot['opt']['rangex'][0], 1, plot['opt']['rangex'][1], 1)
        line.SetLineColor(13)
        line.SetLineStyle(7)
        line.Draw('same pe')

    # Compute the relative uncertainties
    if plot['relunc']['enable']:
        pad = cPlot.cd(panels['relunc'])
        pad.SetGridx(plot['relunc']['gridx'])
        pad.SetGridy(plot['relunc']['gridy'])
        pad.SetLogx(plot["relunc"]["logx"])
        pad.SetLogy(plot["relunc"]["logy"])
        x1 = plot['opt']['rangex'][0]
        y1 = plot['relunc']['rangey'][0]
        x2 = plot['opt']['rangex'][1]
        y2 = plot['relunc']['rangey'][1]
        pad.SetLeftMargin(0.16)
        frame = pad.DrawFrame(x1, y1, x2, y2, TranslateToLatex(plot['opt']['title']))
        frame.GetYaxis().SetTitle('Relative uncertainty (%)')

        for iObj, inObj in enumerate(inObjs):
            if isinstance(inObj, TH1):
                nBins = inObj.GetNbinsX()
                hRelUnc = TH1D(f'hRelUnc_{iObj}', '', nBins, inObj.GetXaxis().GetXmin(), inObj.GetXaxis().GetXmax())

                for iBin in range(hRelUnc.GetNbinsX() + 1):
                    if inObj.GetBinContent(iBin) > 0:
                        hRelUnc.SetBinContent(iBin, 100 * inObj.GetBinError(iBin)/inObj.GetBinContent(iBin))
                        hRelUnc.SetBinError(iBin, 0)
            elif isinstance(inObj, TGraphErrors):
                hRelUnc = TGraphErrors(1)
                for iPoint in range(inObj.GetN()):
                    x = inObj.GetPointX(iPoint)
                    y = inObj.GetPointY(iPoint)
                    xUnc = inObj.GetErrorX(iPoint)
                    yUnc = inObj.GetErrorY(iPoint)

                    hRelUnc.SetPoint(iPoint, x,  100 * yUnc/y)
                    hRelUnc.SetPointError(iPoint, xUnc, 0)
            elif isinstance(inObj, TGraphAsymmErrors):
                hRelUnc = TGraphAsymmErrors(1)
                for iPoint in range(inObj.GetN()):
                    x = inObj.GetPointX(iPoint)
                    y = inObj.GetPointY(iPoint)
                    yUnc = inObj.GetErrorY(iPoint)  # Computes the average of the uncertainties
                    xUncUpper = inObj.GetErrorXhigh(iPoint)
                    xUncLower = inObj.GetErrorXlow(iPoint)

                    hRelUnc.SetPoint(iPoint, x,  100 * yUnc/y)
                    hRelUnc.SetPointError(iPoint, xUncLower, xUncUpper, 0, 0)
            else:
                log.error('Relative uncertainties for type %s are not implemented. Skipping this object', type(inObj))
                continue

            hRelUnc.SetLineColor(inObj.GetLineColor())
            hRelUnc.SetMarkerColor(inObj.GetMarkerColor())
            if isinstance(inObj, TH1):
                hRelUnc.DrawCopy('same hist')
            elif isinstance(inObj, TGraph):
                hRelUnc.Draw('same p')

    cPlot.Modified()
    cPlot.Update()

    # save canvas
    for ext in plot["opt"]["ext"]:
        cPlot.SaveAs(f'{os.path.splitext(plot["output"])[0]}.{ext}')
