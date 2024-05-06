'''
Script to perform the fit on a correlation function.
The output file is CustomNameFromYaml_suffix.root

Usage:
python3 CorrelationFitter.py cfg.yml

'''

import os
import argparse
import yaml

from ROOT import TFile, TCanvas, gInterpreter, TH1D
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/CorrelationFitter.hxx"')
from ROOT import CorrelationFitter

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
    if('rebin' in fitcf):
        fitHisto.Rebin(fitcf['rebin'])

    # fit range
    lowFitRange = fitcf['fitrange'][0]
    uppFitRange = fitcf['fitrange'][1]
    if('rejectrange' in fitcf):
        lowRejectRange = fitcf['rejectrange'][0]
        uppRejectRange = fitcf['rejectrange'][1]
        modelFitters.append(CorrelationFitter(fitHisto, lowFitRange, uppFitRange, lowRejectRange, uppRejectRange))
    else: 
        modelFitters.append(CorrelationFitter(fitHisto, lowFitRange, uppFitRange))
        
    # directory of the fit
    oFile.mkdir(fitcf['fitname'])
    oFile.cd(fitcf['fitname'])

    baselineIdx = -1
    linesThickness = fitcf['linethick']
    compsToFile = []
    onBaseline = []
    shifts = []
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
        onBaseline.append(term['onbaseline'])
        if('shift' in term):
            shifts.append(term['shift'])
        else: 
            shifts.append(0.)
        
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
                histoFuncFiles = []
                histoFuncHistoPars = []
                for histoFuncFile, histoFuncHistoParsPath in zip(term['histofuncfile'], term['histofuncpath']):
                    histoFuncFiles.append(TFile(histoFuncFile))
                    histoFuncHistoPars.append(Load(histoFuncFiles[-1], histoFuncHistoParsPath))

            initPars = []
            for key in term['params']:
                if("fromfunct" == term['params'][key][0]):
                    initPars.append((key, histoFuncHistoPars[term['params'][key][1]].GetBinContent(term['params'][key][2]), 0, -1))
                else:
                    initPars.append((key, term['params'][key][0], term['params'][key][1], 
                                     term['params'][key][2]))

            modelFitters[-1].Add(term['func'], initPars, term['addmode'])
                
    # perform the fit and save the result
    if('globnorm' in fitcf):
        modelFitters[-1].AddGlobNorm('globnorm', fitcf['globnorm'][0], fitcf['globnorm'][1], fitcf['globnorm'][2])    
    oFile.cd(fitcf['fitname'])
    modelFitters[-1].BuildFitFunction()
    oFile.cd(fitcf['fitname'])
    modelFitters[-1].Fit()
    hChi2DOF = TH1D('hChi2DOF', 'hChi2DOF', 1, 0, 1)
    hChi2DOF.Fill(0.5, modelFitters[-1].GetChi2Ndf())
    print('Chi2 / DOF: ' + str(modelFitters[-1].GetChi2Ndf()))
    print('\n\n')
    cFit = TCanvas('cFit', '', 600, 600)
    if('drawsumcomps' in fitcf):
        modelFitters[-1].Draw(cFit, legLabels, fitcf['legcoords'], onBaseline,
                              linesThickness, baselineIdx, fitcf['drawsumcomps'])
    else:
        modelFitters[-1].Draw(cFit, legLabels, fitcf['legcoords'], linesThickness,
                              onBaseline, shifts, baselineIdx)
    if('debug' in fitcf):
        modelFitters[-1].Debug()
        
    cFit.Write()
    fitHisto.Write()
    fitFunction = modelFitters[-1].GetFitFunction()
    modelFitters[-1].SaveFitPars().Write()
    modelFitters[-1].SaveScatPars().Write()
    modelFitters[-1].SaveFreeFixPars().Write()
    fitFunction.Write()
    hChi2DOF.Write()
    for iCompToFile, compToFile in enumerate(compsToFile):
        if('spline' not in fitcf['model'][iTerm] and 'template' not in fitcf['model'][iTerm]):
            modelFitters[-1].GetComponent(compToFile, baselineIdx, fitcf['model'][iTerm]['onbaseline']).Write(term['func'])
        modelFitters[-1].GetComponentPars(compToFile).Write('h' + fitcf['model'][compToFile]['func'][0].upper() + 
                                                            fitcf['model'][compToFile]['func'][1:])
    for iPar in range(fitFunction.GetNpar()):
        cfg[f'Fit n°{iFit}, par {iPar}'] = fitFunction.GetParName(iPar) + ", " + str(fitFunction.GetParameter(iPar))
    #pdfFileName = fitcf['fitname'] + cfg["suffix"] + ".pdf"
    #pdfFilePath = os.path.join(cfg['odir'], pdfFileName) 
    #cFit.SaveAs(pdfFilePath)
    
#os.makedirs(cfg['odir'], exist_ok=True)
#oFileNameCfg = os.path.join(cfg['odir'], oFileBaseName + '_cfg.txt')            
#with open(oFileNameCfg, 'w') as file:
#    for line in fileLines:
#        file.write(line + '\n')
#with open(oFileNameCfg, 'w') as outfile:
#    yaml.dump(cfg, outfile, default_flow_style=False)

oFile.Close()
#print(f'Config saved in {oFileNameCfg}')
print(f'output saved in {oFileName}')