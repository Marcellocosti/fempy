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

oFileNameCfg = os.path.join(cfg['odir'], oFileBaseName + '_cfg.txt')
with open(oFileNameCfg, 'w') as file:
    for line in fileLines:
        file.write(line)
        file.write('\n')

modelFitters = []

# for loop over the correlation functions
for fitcf in cfg['fitcfs']:
    
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

    compsToFile = []
    # for loop over the functions entering in the model
    for iTerm, term in enumerate(fitcf['model']):
        
        if('isbaseline' in term):
            if(term['isbaseline']):
                modelFitters[-1].SetBaselineIdx(iTerm)
        if('savetofile' in term):
            compsToFile.append(iTerm)
        
        # fit function parameters initialization
        #initPars = []
    
        if('template' in term):
            print(term['template'])
            histoFile = TFile(term['histofile'])
            splinedHisto = ChangeUnits(Load(histoFile, term['histopath']), 1000)
            if('rebin' in term):
                splinedHisto.Rebin(term['rebin'])
            initPars = [(key, term['params'][key][0], term['params'][key][1], 
                         term['params'][key][2]) for key in term['params']]
            modelFitters[-1].AddSplineHisto(term['template'], splinedHisto, initPars, term['addmode'], term['onbaseline'])
            cSplinedHisto = TCanvas(f'c{term["template"]}', '', 600, 600)
            modelFitters[-1].DrawSpline(cSplinedHisto, splinedHisto)
            oFile.cd(fitcf['fitname'])
            cSplinedHisto.Write()
        
        elif('func' in term):
            print(term['func'])
            #if('spline3' in term['func']):
            #    for nKnot, xKnot in enumerate(term['xknots']):
            #        initPars.append([f'xKnot{nKnot}', xKnot, xKnot, xKnot])
            #    for nKnot, xKnot in enumerate(term['xknots']):
            #        nBin = prefitHisto.FindBin(xKnot)
            #        yKnot = prefitHisto.GetBinContent(nBin)
            #        initPars.append([f'yKnot{nKnot}', yKnot, yKnot - (yKnot/100)*30, yKnot + (yKnot/100)*30])
            #else:
            
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
            modelFitters[-1].Add(term['func'], initPars, term['addmode'], term['onbaseline'])
        
        else:
            pass
                
    # perform the fit and save the result
    oFile.cd(fitcf['fitname'])
    modelFitters[-1].BuildFitFunction()
    oFile.cd(fitcf['fitname'])
    modelFitters[-1].Fit()
    cFit = TCanvas('cFit', '', 600, 600)
    if('drawsumcomps' in fitcf):
        modelFitters[-1].Draw(cFit, fitcf['drawsumcomps'])
    else:
        modelFitters[-1].Draw(cFit)
        
    modelFitters[-1].DrawLegend(cFit, fitcf['legcoords'][0], fitcf['legcoords'][1], fitcf['legcoords'][2], 
                             fitcf['legcoords'][3], fitcf['legentries'])
    cFit.Write()
    fitHisto.Write()
    fitFunction = modelFitters[-1].GetFitFunction()
    fitFunction.Write()
    for compToFile in compsToFile:
        if('spline' not in fitcf['model'][iTerm]['func']):
            modelFitters[-1].GetComponent(compToFile).Write(term['func'])
        modelFitters[-1].GetComponentPars(compToFile).Write('h' + fitcf['model'][compToFile]['func'][0].upper() + 
                                                            fitcf['model'][compToFile]['func'][1:])
    
    with open(oFileNameCfg, 'a') as file:
        file.write('-----------------------------------')
        file.write('\n')
        file.write('Parameters obtained from the fit')    
        file.write('\n')
        for iPar in range(fitFunction.GetNpar()):
            file.write(fitFunction.GetParName(iPar) + ": " + str(fitFunction.GetParameter(iPar)))
            file.write('\n')
    modelFitters[-1].Debug()
    #pdfFileName = fitcf['fitname'] + cfg["suffix"] + ".pdf"
    #pdfFilePath = os.path.join(cfg['odir'], pdfFileName) 
    #cFit.SaveAs(pdfFilePath)
            

oFile.Close()
print(f'Config saved in {oFileNameCfg}')
print(f'output saved in {oFileName}')
#