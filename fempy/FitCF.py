'''
Script to perform the fit on a correlation function.
The output file is FitCF_suffix.root

Usage:
python3 fitCF.py cfg.yml

'''

import os
import argparse
import yaml

from ROOT import TFile, TCanvas, gInterpreter, TF1, TDatabasePDG
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/CorrelationFitter.hxx"')
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/Fitter.hxx"')
from ROOT import CorrelationFitter, BreitWigner, Fitter

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
inFileData = TFile(cfg['cfinfile'])
inFileMC = TFile(cfg['mcinfile'])

# Define the output file
oFileBaseName = 'FitCF'
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

oFileNameCfg = '/home/mdicostanzo/an/LPi/fitstry/' + oFileBaseName + '_cfg.txt' 
with open(oFileNameCfg, 'w') as file:
    for line in fileLines:
        file.write(line)
        file.write('\n')

cfFitters = []
preFitters = []
prePreFitters = []

# for loop over the correlation functions
for nFit, fitcf in enumerate(cfg['fitcfs']):

    # change unity of measure of histograms from GeV to MeV
    dataCF = ChangeUnits(Load(inFileData, fitcf['cfpath']), 1000)
    mcCF = ChangeUnits(Load(inFileMC, fitcf['cfpath']), 1000)

    # fit range
    lowFitRange = fitcf['fitrange'][0]
    uppFitRange = fitcf['fitrange'][1]
    # directory of the fit
    oFile.mkdir(fitcf['fitname'])
    oFile.cd(fitcf['fitname'])

    cfFitters.append(CorrelationFitter(dataCF, mcCF, lowFitRange, uppFitRange))

    ancIdx = 0
    # for loop over the functions entering in the model
    for funcIdx, func in enumerate(fitcf['model']):
        print(func['funcname'])
        if('isbaseline' in func):
            if(func['isbaseline']):
                cfFitters[-1].SetBaselineIdx(funcIdx)
                #cfFitters[-1].SetBaselineIdx(6)
        if('isancestor' in func and ancIdx==0):
            if(func['isancestor']):
                cfFitters[-1].SetAncestorIdx(ancIdx)
        
        # fit function parameters initialization
        initPars = []
        
        #print('CIAO')
        if('splinehisto' in func['funcname']):
            #print('CIAO1')
            histoFile = TFile(func['histofile'])
            splinedHisto = ChangeUnits(Load(histoFile, func['histopath']), 1000)
            if('rebin' in func):
                splinedHisto.Rebin(func['rebin'])
            initPars = [(func['norm'][0], func['norm'][1], func['norm'][2], func['norm'][3])]
            cfFitters[-1].AddSplineHisto(func['funcname'], splinedHisto, initPars, func['addmode'], func['onbaseline'])
            cSplinedHisto = TCanvas(f'c{func["funcname"]}', '', 600, 600)
            cfFitters[-1].DrawSpline(cSplinedHisto, splinedHisto)
            oFile.cd(fitcf['fitname'])
            cSplinedHisto.Write()
            continue
        
        # check if the function is to be prefitted
        if('prefitcomp' not in func and func['prefitfile'] is not None):
                
                #print('CIAO2')
                prefitFile = TFile(func['prefitfile'])
                prefitHisto = ChangeUnits(Load(prefitFile, func['prefitpath']), 1000)

                # prefit function parameters initialization
                preInitPars = []
                lowPrefitRange = func['prefitrange'][0]
                uppPrefitRange = func['prefitrange'][1]

                preFitters.append(Fitter(prefitHisto, lowPrefitRange, uppPrefitRange))
                cPrefit = TCanvas(f'cPrefit_{func["funcname"]}', '', 600, 600)
                # the splines need a different implementation
                if('spline3' in func['funcname']):
                    for nKnot, xKnot in enumerate(func['xknots']):
                        preInitPars.append([f'xKnot{nKnot}', xKnot, xKnot, xKnot])
                    for nKnot, xKnot in enumerate(func['xknots']):
                        nBin = prefitHisto.FindBin(xKnot)
                        yKnot = prefitHisto.GetBinContent(nBin)
                        preInitPars.append([f'yKnot{nKnot}', yKnot, yKnot - (yKnot/100)*func['bounds'], 
                                            yKnot + (yKnot/100)*func['bounds']])
                else:
                    preInitPars = [(func[f'p{iPar}'][0], func[f'p{iPar}'][1], func[f'p{iPar}'][2], 
                                    func[f'p{iPar}'][3]) for iPar in range(func['npars'])]
                preFitters[-1].Add(func['funcname'], preInitPars)

                # include other functions of the prefitting model
                for prefitFunc in func['prefitmodel']:
                    prefitInitPars = []
                    if('spline3' in prefitFunc['funcname']):
                        for nKnot, xKnot in enumerate(prefitFunc['xknots']):
                            prefitInitPars.append([f'xKnot{nKnot}', xKnot, xKnot, xKnot])
                        for nKnot, xKnot in enumerate(prefitFunc['xknots']):
                            nBin = prefitHisto.FindBin(xKnot)
                            yKnot = prefitHisto.GetBinContent(nBin)
                            prefitInitPars.append([f'yKnot{nKnot}', yKnot, yKnot - (yKnot/100)*prefitFunc['bounds'], 
                                               yKnot + (yKnot/100)*prefitFunc['bounds']])
                    else: 
                        prefitInitPars = [(prefitFunc[f'p{iPar}'][0], prefitFunc[f'p{iPar}'][1], prefitFunc[f'p{iPar}'][2], 
                                           prefitFunc[f'p{iPar}'][3]) for iPar in range(prefitFunc['npars'])]

                    preFitters[-1].Add(prefitFunc['funcname'], prefitInitPars)

                print('PREFITTING')
                preFitters[-1].Fit()
                preFitters[-1].Draw(cPrefit)
                oFile.cd(fitcf['fitname'])
                preFitters[-1].GetFunction().Write()
                cPrefit.Write()

                prefitRes = preFitters[-1].GetFunction()

                # save prefit results for fit parameter initialization
                if(func['fixprefit']):
                    lowBound = 1
                    uppBound = 1
                else:
                    lowBound = 0.8
                    uppBound = 1.2

                if('spline3' in func['funcname']):
                    nKnots = int(int(func['npars'])/2)
                    for iPar in range(nKnots):
                        initPars.append([f'xKnot{iPar}', prefitRes.GetParameter(iPar), 
                                         prefitRes.GetParameter(iPar), prefitRes.GetParameter(iPar)])            
                    for iPar in range(nKnots):
                        initPars.append([f'yKnot{iPar}', prefitRes.GetParameter(iPar + nKnots), 
                                         prefitRes.GetParameter(iPar + nKnots) * lowBound, prefitRes.GetParameter(iPar + nKnots) * uppBound])

                else:
                    for iPar in range(func['npars']):
                        if(func[f'p{iPar}'][2] > func[f'p{iPar}'][3]):
                            initPars.append([func[f'p{iPar}'][0], func[f'p{iPar}'][1], func[f'p{iPar}'][2], func[f'p{iPar}'][3]])
                        else:
                            if(prefitRes.GetParameter(iPar) >= 0):
                                initPars.append([func[f'p{iPar}'][0], prefitRes.GetParameter(iPar), 
                                     prefitRes.GetParameter(iPar) * lowBound, prefitRes.GetParameter(iPar) * uppBound])
                            else:
                                initPars.append([func[f'p{iPar}'][0], prefitRes.GetParameter(iPar), 
                                                 prefitRes.GetParameter(iPar) * uppBound, prefitRes.GetParameter(iPar) * lowBound])

        # no prefit case
        else:
            #print('CIAO3')
            if('spline3' in func['funcname']):
                for nKnot, xKnot in enumerate(func['xknots']):
                    initPars.append([f'xKnot{nKnot}', xKnot, xKnot, xKnot])
                for nKnot, xKnot in enumerate(func['xknots']):
                    nBin = prefitHisto.FindBin(xKnot)
                    yKnot = prefitHisto.GetBinContent(nBin)
                    initPars.append([f'yKnot{nKnot}', yKnot, yKnot - (yKnot/100)*30, yKnot + (yKnot/100)*30])
            else:
                if('splinehisto' in func['funcname']):
                    initPars = [(['splinecoeff', 1, 0, -1])]
                else:
                    initPars = [(func[f'p{iPar}'][0], func[f'p{iPar}'][1], func[f'p{iPar}'][2], 
                                 func[f'p{iPar}'][3]) for iPar in range(func['npars'])]
            
            #if('prefitcomp' in func):
            #    if(func['prefitcomp']):
            #        prefitFile = TFile(func['prefitfile'])
            #        prefitHisto = ChangeUnits(Load(prefitFile, func['prefitpath']), 1000)
            #
            #        # prefit function parameters initialization
            #        lowPrefitRange = func['prefitrange'][0]
            #        uppPrefitRange = func['prefitrange'][1]
            #        
            #        cfFitters[-1].PrefitComponent(func['compfuncname'], prefitHisto, initPars, lowPrefitRange, 
            #                                      uppPrefitRange, func['startnewpar'])

        #print('Init pars')
        #print(initPars)

        if('lambdapar' in func):
            lambdaParam = [("lambdapar_" + func['funcname'], func['lambdapar'], 0, -1)]
            initPars = lambdaParam + initPars
            print('ADDING FUNCTION')
            print(initPars)
            cfFitters[-1].Add(func['funcname'], initPars, func['addmode'], func['onbaseline'])
        if('lambdagen' in func):
            lambdaGen = [("lambda_gen_" + func['funcname'], func['lambdagen'], 0, -1)]
            initPars = lambdaGen + initPars
            print('ADDING FUNCTION')
            print(initPars)            
            cfFitters[-1].Add(func['funcname'], initPars, func['addmode'], func['onbaseline'])
        if('norm' in func):
            normParam = [(func['norm'][0], func['norm'][1], func['norm'][2], func['norm'][3])]
            initPars = normParam + initPars
            print('ADDING FUNCTION')
            print(initPars)
            if('splinehisto' in func['funcname']):
                cfFitters[-1].Add('pol0', initPars, func['addmode'],  func['onbaseline'])
            else:    
                cfFitters[-1].Add(func['funcname'], initPars, func['addmode'], func['onbaseline'])
                
    # perform the fit and save the result
    oFile.cd(fitcf['fitname'])
    #cfFitters[-1].PrefitMC()
    #cPrefit = TCanvas('cPrefit', '', 600, 600)
    #cfFitters[-1].DrawPrefit(cPrefit)
    #cPrefit.Write()
    cfFitters[-1].BuildFitFunction()
    for funcIdx, func in enumerate(fitcf['model']):
        #print('CHECK PREFIT COMP')
        if('prefitcomp' in func):
            if(func['prefitcomp']):
                #print('PREFITTING COMP')
                prefitFile = TFile(func['prefitfile'])
                prefitHisto = ChangeUnits(Load(prefitFile, func['prefitpath']), 1000)
                lowPrefitRange = func['prefitrange'][0]
                uppPrefitRange = func['prefitrange'][1]
                lowRejectRange = func['rejectrange'][0]
                uppRejectRange = func['rejectrange'][1]
                startPar = func['startnewpar']
                nParsComp = func['nparscomp']
                compFuncName = func['compfuncname']
                
                cCompPrefit = TCanvas('cCompPrefit_' + func['funcname'], '', 600, 600)
                #print('READY FOR THE COMPONENT PREFIT')
                oFile.cd(fitcf['fitname'])
                #print('CHANGED TO FILE DIRECTORY')
                cfFitters[-1].PrefitComponent(cCompPrefit, prefitHisto, compFuncName, startPar, nParsComp, 
                                              lowPrefitRange, uppPrefitRange, lowRejectRange, uppRejectRange)
                #print('COMPONENT PREFITTED')
                cCompPrefit.Write()
                #print('WRITTEN TO FILE')
    
    print('READY FOR THE FIT')
    oFile.cd(fitcf['fitname'])
    cfFitters[-1].Fit().Write()
    cFit = TCanvas('cFit', '', 600, 600)
    if('drawsumcomps' in fitcf):
        cfFitters[-1].Draw(cFit, fitcf['drawsumcomps'])
    else:
        cfFitters[-1].Draw(cFit)
    cfFitters[-1].DrawLegend(cFit, fitcf['legcoords'][0], fitcf['legcoords'][1], fitcf['legcoords'][2], 
                             fitcf['legcoords'][3], fitcf['legentries'])
    cFit.Write()
    dataCF.Write()
    fitFunction = cfFitters[-1].GetFunction()
    fitFunction.Write()
    with open(oFileNameCfg, 'a') as file:
        file.write('-----------------------------------')
        file.write('\n')
        file.write('Parameters obtained from the fit')    
        file.write('\n')
        for iPar in range(fitFunction.GetNpar()):
            file.write(fitFunction.GetParName(iPar) + ": " + str(fitFunction.GetParameter(iPar)))
            file.write('\n')
    
    print('Evaluation of the fit function: ' + str(cfFitters[-1].GetFunction().Eval(100)) )

    #pdfFileName = fitcf['fitname'] + cfg["suffix"] + ".pdf"
    #pdfFilePath = os.path.join(cfg['odir'], pdfFileName) 
    #cFit.SaveAs(pdfFilePath)
    #cfFitters[-1].Debug()

oFile.Close()
print(f'Config saved in {oFileNameCfg}')
print(f'output saved in {oFileName}')