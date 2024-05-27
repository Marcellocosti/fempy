'''
Script to perform the fit on a correlation function.
The output file is CustomNameFromYaml_suffix.root

Usage:
python3 CorrelationFitter.py cfg.yml

'''

import os
import argparse
import yaml

from ROOT import TFile, TCanvas, gInterpreter, TH1, TH1D, TSpline3

from fempy import logger as log
from fempy.utils.io import Load
from fempy.utils.analysis import ChangeUnits

parser = argparse.ArgumentParser()
parser.add_argument('cfg', default='')
parser.add_argument('--debug', default=False, action='store_true')
args = parser.parse_args()

if args.debug:
    gInterpreter.ProcessLine(f'#define LOG_LEVEL 1')
    log.setLevel(1)

gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/CorrelationFitter.hxx"')
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/DrawFitFuncts.hxx"')
from ROOT import CorrelationFitter, DrawFitFuncts

# Load yaml file
with open(args.cfg, "r") as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        log.critical('Yaml configuration could not be loaded. Is it properly formatted?')

# Load input file with data and mc CF
inFileFit = TFile(cfg['infile'])

# Define the output file
oFileBaseName = cfg['ofilebasename']
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'
oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)

fileLines = []
fitCfgIdxs = []
with open(args.cfg, 'r') as file:
    for lineNumber, line in enumerate(file, start=1):
        if('cfpath' in line):
            fitCfgIdxs.append(lineNumber)
        firstNonBlankChar = next((char for char in line if not char.isspace()), None)
        if(firstNonBlankChar == '#'): continue
        fileLines.append(line.rstrip())  # Strip to remove leading/trailing whitespaces

modelFitters = []
drawFits = []

# for loop over the correlation functions
for iFit, fitcf in enumerate(cfg['fitcfs']):
    
    # change unity of measure of histograms from GeV to MeV
    fitHisto = ChangeUnits(Load(inFileFit, fitcf['cfpath']), 1000)
    if fitcf.get('rebin'):
        fitHisto.Rebin(fitcf['rebin'])

    # fit range
    lowFitRange = fitcf['fitrange'][0]
    uppFitRange = fitcf['fitrange'][1]
    if fitcf.get('rejectrange'):
        lowRejectRange = fitcf['rejectrange'][0]
        uppRejectRange = fitcf['rejectrange'][1]
        modelFitters.append(CorrelationFitter(fitHisto, lowFitRange, uppFitRange, lowRejectRange, uppRejectRange))
    else: 
        modelFitters.append(CorrelationFitter(fitHisto, lowFitRange, uppFitRange))
        
    # drawing class constructor
    if fitcf.get('drawrange'):
        lowDrawRange = fitcf['drawrange'][0]
        uppDrawRange = fitcf['drawrange'][1]
        drawFits.append(DrawFitFuncts(fitHisto, lowDrawRange, uppDrawRange))
    else: 
        drawFits.append(DrawFitFuncts(fitHisto, lowFitRange, uppFitRange))

    # directory of the fit
    oFile.mkdir(fitcf['fitname'])
    oFile.cd(fitcf['fitname'])

    baselineIdx = -1
    linesThickness = fitcf['linethick']
    compsToFile = []
    onBaseline = []
    shifts = []
    multNorm = []
    multGlobNorm = []
    saveSubComps = []
    normsSubComps = []
    normsSubCompsLabels = []
    subCompsMothers = []
    subComps = []
    legLabels = []
    legLabels.append(fitcf['datalabel'])
    legLabels.append(fitcf['fitfunclabel'])
    # for loop over the functions entering in the model
    for iTerm, term in enumerate(fitcf['model']):
        
        legLabels.append(term['legentry'])
        onBaseline.append(term.get('onbaseline', 0))
        shifts.append(term.get('shift', 0))
        multNorm.append(term.get('multnorm', 0))
        multGlobNorm.append(term.get('multglobnorm', 0))

        if term.get('subcomps'):
            saveSubComps.append(iTerm-1)
            normsSubComps.append(term['normssubcomps'])
            normsSubCompsLabels.append(term['normssubcompslabels'])
            subCompsMothers.append(term['func'])
            subComps.append(term['subcomps'])
            for iSubComp in range(len(term['subcomps'])):
                print('Reading subcomps!')
                onBaseline.append(term['sub_onbaseline'][iSubComp])                
                legLabels.append(term['sub_legentry'][iSubComp])
                shifts.append(term.get('sub_shifts', [0.]*len(term['subcomps']))[iSubComp])
                multNorm.append(term.get('sub_multnorm', [1]*len(term['subcomps']))[iSubComp])
                multGlobNorm.append(term.get('sub_multglobnorm', [1]*len(term['subcomps']))[iSubComp])

        if term.get('isbaseline'):
            drawFits[-1].SetBasIdx(iTerm)
            baselineIdx = iTerm
                
        if term.get('template'):
            drawFits[-1].AddFitCompName(term['template'])
            templFile = TFile(term['templfile'])
            splinedTempl = Load(templFile, term['templpath'])
            if(isinstance(splinedTempl, TH1)):
                splinedTempl = ChangeUnits(splinedTempl, 1000)
                if term.get('rebin'):
                    splinedTempl.Rebin(term['rebin'])
            initPars = [(key, term['params'][key][0], term['params'][key][1], 
                         term['params'][key][2]) for key in term['params']]    
            modelFitters[-1].Add(term['template'], splinedTempl, initPars, term['addmode'])
            cSplinedTempl = TCanvas(f'c{term["template"]}', '', 600, 600)
            modelFitters[-1].DrawSpline(cSplinedTempl, splinedTempl)
            drawFits[-1].AddSplineHisto(splinedTempl)
            oFile.cd(fitcf['fitname'])
            cSplinedTempl.Write()
        
        elif term.get('spline'):
            drawFits[-1].AddFitCompName(term['spline'])
            histoFile = TFile(term['histofile'])
            toBeSplinedHisto = ChangeUnits(Load(histoFile, term['histopath']), term['changeunits'])
            drawFits[-1].AddSplineHisto(toBeSplinedHisto)
            initPars = []
            normPar = list(term['params'].keys())[0]
            initPars.append((normPar, term['params'][normPar][0], 
                             term['params'][normPar][1], term['params'][normPar][2]))
            for nKnot, xKnot in enumerate(term['params']['xknots']):
                initPars.append([f'xKnot{nKnot}', xKnot, xKnot, xKnot])
            for nKnot, xKnot in enumerate(term['params']['xknots']):
                nBin = toBeSplinedHisto.FindBin(xKnot)
                yKnot = toBeSplinedHisto.GetBinContent(nBin)
                initPars.append([f'yKnot{nKnot}', yKnot, yKnot - (yKnot/100)*term['yboundperc'], yKnot + (yKnot/100)*term['yboundperc']])
            modelFitters[-1].Add(term['spline'], initPars, term['addmode'])
        
        elif term.get('func'):
            drawFits[-1].AddFitCompName(term['func'])
            if term.get('subcomps'):
                for iSubComp in range(len(term['subcomps'])):
                    drawFits[-1].AddFitCompName(term['subcomps'][iSubComp])

            if term.get('fixparsfromfuncts'):
                histoFuncFiles = []
                histoFuncHistoPars = []
                for histoFuncFile, histoFuncHistoParsPath in zip(term['histofuncfile'], term['histofuncpath']):
                    histoFuncFiles.append(TFile(histoFuncFile))
                    histoFuncHistoPars.append(Load(histoFuncFiles[-1], histoFuncHistoParsPath))

            initPars = []
            for key in term['params']:
                if "fromfunct" == term['params'][key][0]:
                    initPars.append((key, histoFuncHistoPars[term['params'][key][1]].GetBinContent(term['params'][key][2]), 0, -1))
                else:
                    initPars.append((key, term['params'][key][0], term['params'][key][1], 
                                     term['params'][key][2]))

            modelFitters[-1].Add(term['func'], initPars, term['addmode'])
                
    # perform the fit and save the result
    if fitcf.get('globnorm'):
        modelFitters[-1].AddGlobNorm('globnorm', fitcf['globnorm'][0], fitcf['globnorm'][1], fitcf['globnorm'][2])    
        drawFits[-1].SetGlobNorm(True)    
    
    oFile.cd(fitcf['fitname'])
    modelFitters[-1].BuildFitFunction()
    modelFitters[-1].Fit()
    fitHisto.Write()
    fitFunction = modelFitters[-1].GetFitFunction()
    fitFunction.Write()
    hChi2DOF = TH1D('hChi2DOF', 'hChi2DOF', 1, 0, 1)
    hChi2DOF.Fill(0.5, modelFitters[-1].GetChi2Ndf())
    hChi2DOF.Write()
    modelFitters[-1].SaveFreeFixPars().Write()
    modelFitters[-1].SaveFitPars().Write()
    modelFitters[-1].PullDistribution().Write()
    if fitcf.get('isfitcf'):
        modelFitters[-1].SaveScatPars().Write()

    hAllCompsParHisto = modelFitters[-1].SaveFitParsSplitComponents(subCompsMothers, saveSubComps, normsSubComps, subComps, normsSubCompsLabels)
    hAllCompsParHisto.Write()    
    cFit = TCanvas('cFit', '', 600, 600)
    drawFits[-1].SetTotalFitFunc(fitFunction)
    drawFits[-1].SetParHist(hAllCompsParHisto)
    if fitcf.get('drawsumcomps'):
        drawFits[-1].EvaluateToBeDrawnComponents(onBaseline, multNorm, multGlobNorm, shifts,
                                             baselineIdx, fitcf['drawsumcomps'])
    else:        
        drawFits[-1].EvaluateToBeDrawnComponents(onBaseline, multNorm, multGlobNorm, shifts,
                                             baselineIdx)
    
    drawFits[-1].Draw(cFit, legLabels, fitcf['legcoords'], linesThickness)
        
    cFit.Write()

oFile.Close()
print(f'output saved in {oFileName}')