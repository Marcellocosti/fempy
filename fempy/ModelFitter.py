'''
Script to perform the fit on a correlation function.
The output file is CustomNameFromYaml_suffix.root

Usage:
python3 ModelFitter.py cfg.yml

'''

import os
import argparse
import yaml

from ROOT import TFile, TCanvas, gInterpreter, TH1D
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/ModelFitter.hxx"')
from ROOT import ModelFitter

from fempy import logger as log
from fempy.utils.io import Load
from fempy.utils.analysis import ChangeUnits

parser = argparse.ArgumentParser()
parser.add_argument('cfg', default='')
parser.add_argument('--debug', default=False, action='store_true')
args = parser.parse_args()

if args.debug:
    log.setLevel(1)

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

# for loop over the correlation functions
for iFit, fitcf in enumerate(cfg['fitcfs']):
    
    # change unity of measure of histograms from GeV to MeV
    fitHisto = ChangeUnits(Load(inFileFit, fitcf['cfpath']), 1000)

    # fit range
    lowFitRange = fitcf['fitrange'][0]
    uppFitRange = fitcf['fitrange'][1]
    if('rejectrange' in fitcf):
        lowRejectRange = fitcf['rejectrange'][0]
        uppRejectRange = fitcf['rejectrange'][1]
        modelFitters.append(ModelFitter(fitHisto, lowFitRange, uppFitRange, lowRejectRange, uppRejectRange))
    else: 
        modelFitters.append(ModelFitter(fitHisto, lowFitRange, uppFitRange))
        
    # directory of the fit
    oFile.mkdir(fitcf['fitname'])
    oFile.cd(fitcf['fitname'])

    baselineIdx = -1
    linesThickness = fitcf['linethick']
    compsToFile = []
    onBaseline = []
    colors = []
    legLabels = []
    legLabels.append(fitcf['datalabel'])
    legLabels.append(fitcf['fitfunclabel'])
    # for loop over the functions entering in the model
    for iTerm, term in enumerate(fitcf['model']):
        
        if('isbaseline' in term):
            if(term['isbaseline']):
                baselineIdx = iTerm
        if('savetofile' in term):
            compsToFile.append(iTerm)
        
        legLabels.append(term['legentry'])
        colors.append(term['linecolor'])
        onBaseline.append(term['onbaseline'])

        if('template' in term):
            histoFile = TFile(term['histofile'])
            splinedHisto = ChangeUnits(Load(histoFile, term['histopath']), 1000)
            if('rebin' in term):
                splinedHisto.Rebin(term['rebin'])
            initPars = [(key, term['params'][key][0], term['params'][key][1], 
                         term['params'][key][2]) for key in term['params']]
            modelFitters[-1].Add(term['template'], splinedHisto, initPars, term['addmode'])
            cSplinedHisto = TCanvas(f'c{term["template"]}', '', 600, 600)
            modelFitters[-1].DrawSpline(cSplinedHisto, splinedHisto)
            oFile.cd(fitcf['fitname'])
            cSplinedHisto.Write()
        
        elif('spline' in term):
            histoFile = TFile(term['histofile'])
            toBeSplinedHisto = ChangeUnits(Load(histoFile, term['histopath']), term['changeunits'])
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
        
        elif('func' in term):
            if('fixparsfromfuncts' in term):
                histoFuncFile = TFile(term['histofuncfile'])
                histoFuncParams = Load(histoFuncFile, term['histofuncpath'])

            initPars = []
            for key in term['params']:
                if("fromfunct" == term['params'][key][0]):
                    initPars.append((key, histoFuncParams.GetBinContent(term['params'][key][1]), 0, -1))
                else:
                    initPars.append((key, term['params'][key][0], term['params'][key][1], 
                                     term['params'][key][2]))

            print(initPars)
            modelFitters[-1].Add(term['func'], initPars, term['addmode'])
                
    # perform the fit and save the result
    oFile.cd(fitcf['fitname'])
    modelFitters[-1].BuildFitFunction()
    oFile.cd(fitcf['fitname'])
    modelFitters[-1].Fit()
    cFit = TCanvas('cFit', '', 600, 600)
    print(colors)
    if('drawsumcomps' in fitcf):
        modelFitters[-1].Draw(cFit, legLabels, fitcf['legcoords'], onBaseline,
                              linesThickness, baselineIdx, fitcf['drawsumcomps'])
    else:
        modelFitters[-1].Draw(cFit, legLabels, fitcf['legcoords'], onBaseline,
                              linesThickness, baselineIdx)
        
    cFit.Write()
    fitHisto.Write()
    fitFunction = modelFitters[-1].GetFitFunction()
    fitFunction.Write()
    for iCompToFile, compToFile in enumerate(compsToFile):
        if('spline' not in fitcf['model'][iTerm] and 'template' not in fitcf['model'][iTerm]):
            modelFitters[-1].GetComponent(compToFile, baselineIdx).Write(term['func'])
        modelFitters[-1].GetComponentPars(compToFile).Write('h' + fitcf['model'][compToFile]['func'][0].upper() + 
                                                            fitcf['model'][compToFile]['func'][1:])
    for iPar in range(fitFunction.GetNpar()):
        cfg[f'Fit n°{iFit}, par {iPar}'] = fitFunction.GetParName(iPar) + ", " + str(fitFunction.GetParameter(iPar))
    modelFitters[-1].Debug()
    #pdfFileName = fitcf['fitname'] + cfg["suffix"] + ".pdf"
    #pdfFilePath = os.path.join(cfg['odir'], pdfFileName) 
    #cFit.SaveAs(pdfFilePath)
    
oFileNameCfg = os.path.join(cfg['odir'], oFileBaseName + '_cfg.txt')            
print(cfg)
with open(oFileNameCfg, 'w') as outfile:
    yaml.dump(cfg, outfile, default_flow_style=False)
    
oFile.Close()
print(f'Config saved in {oFileNameCfg}')
print(f'output saved in {oFileName}')
