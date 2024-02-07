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
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/CorrelationFitterNew.hxx"')
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/Fitter.hxx"')
from ROOT import CorrelationFitterNew, BreitWigner, Fitter

from fempy import logger as log
from fempy.utils.io import Load
from utils.analysis import ChangeUnits

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

oFileNameCfg = '/home/mdicostanzo/an/LPi/fits/' + oFileBaseName + '_cfg.txt' 
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

    cfFitters.append(CorrelationFitterNew(dataCF, mcCF, lowFitRange, uppFitRange))

    # for loop over the functions entering in the model
    for func in fitcf['model']:

        # fit function parameters initialization
        initPars = []

        if('splinehisto' in func['funcname']):
            histoFile = TFile(func['histofile'])
            splinedHisto = ChangeUnits(Load(histoFile, func['histopath']), 1000)
            if('rebin' in func):
                splinedHisto.Rebin(func['rebin'])
            initPars = [(func['norm'][0], func['norm'][1], func['norm'][2], func['norm'][3])]
            cfFitters[-1].AddSplineHisto(func['funcname'], splinedHisto, initPars, func['addmode'])
            cSplinedHisto = TCanvas(f'cSplinedHisto_{func["funcname"]}', '', 600, 600)
            cfFitters[-1].DrawSpline(cSplinedHisto, splinedHisto)
            oFile.cd(fitcf['fitname'])
            cSplinedHisto.Write()
            continue

        # check if the function is to be prefitted
        if(func['prefitfile'] is not None):

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
                        print(xKnot)
                        prefitInitPars.append([f'xKnot{nKnot}', xKnot, xKnot, xKnot])
                    for nKnot, xKnot in enumerate(prefitFunc['xknots']):
                        nBin = prefitHisto.FindBin(xKnot)
                        yKnot = prefitHisto.GetBinContent(nBin)
                        print(yKnot)
                        prefitInitPars.append([f'yKnot{nKnot}', yKnot, yKnot - (yKnot/100)*prefitFunc['bounds'], 
                                           yKnot + (yKnot/100)*prefitFunc['bounds']])
                else: 
                    prefitInitPars = [(prefitFunc[f'p{iPar}'][0], prefitFunc[f'p{iPar}'][1], prefitFunc[f'p{iPar}'][2], 
                                       prefitFunc[f'p{iPar}'][3]) for iPar in range(prefitFunc['npars'])]
    
                preFitters[-1].Add(prefitFunc['funcname'], prefitInitPars)
                print('Prefit model function added')
    
            preFitters[-1].Fit()
            preFitters[-1].Draw(cPrefit)
            oFile.cd(fitcf['fitname'])
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
                    print(func['funcname'] + ' N pars')
                    print(range(func['npars']))
                    print(func['funcname'] + ' pars')
                    print(initPars)

        if('lambdapar' in func):
            lambdaParam = [("lambdapar_" + func['funcname'], func['lambdapar'], 0, -1)]
            initPars = lambdaParam + initPars
            cfFitters[-1].Add(func['funcname'], initPars, func['addmode'])
        if('lambdagen' in func):
            lambdaGen = [("lambda_gen_" + func['funcname'], func['lambdagen'], 0, -1)]
            initPars = lambdaGen + initPars
            cfFitters[-1].Add(func['funcname'], initPars, func['addmode'])
            antiLambdaGen = [("anti_lambda_gen_coeff", 1, 0, -1)]
            antiLambdaGen.append(("anti_lambda_gen_" + func['funcname'], 1-func['lambdagen'], 0, -1))
            cfFitters[-1].Add('pol0', antiLambdaGen, 'sum')
        if('norm' in func):
            normParam = [(func['norm'][0], func['norm'][1], func['norm'][2], func['norm'][3])]
            print(normParam)
            initPars = normParam + initPars
            print(initPars)
            print('\n\n\n')
            if(func['funcname'] == 'splinehisto'):
                cfFitters[-1].Add('pol0', initPars, func['addmode'])
            else:    
                cfFitters[-1].Add(func['funcname'], initPars, func['addmode'])


    # perform the fit and save the result
    cfFitters[-1].Fit()
    cFit = TCanvas('cFit', '', 600, 600)
    cfFitters[-1].Draw(cFit)
    cFit.SaveAs(f'./Try.pdf')
    oFile.cd(fitcf['fitname'])
    cFit.Write()
    dataCF.Write()
    cfFitters[-1].GetFunction().Write()

oFile.Close()
print(f'output saved in {oFileName}')


                # the splines need a different implementation
                #if('prefitfile' in prefitFunc):
                #    prePrefitFile = TFile(prefitFunc['prefitfile'])
                #    prePrefitHisto = ChangeUnits(Load(prePrefitFile, prefitFunc['prefitpath']), 1000)
                #
                #    # prefit function parameters initialization
                #    prePreInitPars = []
                #    lowPrePrefitRange = prefitFunc['prefitrange'][0]
                #    uppPrePrefitRange = prefitFunc['prefitrange'][1]
                #    prePreFitters.append(Fitter(prePrefitHisto, lowPrePrefitRange, uppPrePrefitRange))
                #    #prePreFitters[-1].Add(func['funcname'], preInitPars)
                #
                #    if('spline3' in prefitFunc['funcname']):
                #        for nKnot, xKnot in enumerate(prefitFunc['xknots']):
                #            print(xKnot)
                #            prePreInitPars.append([f'xKnot{nKnot}', xKnot, xKnot, xKnot])
                #        for nKnot, xKnot in enumerate(prefitFunc['xknots']):
                #            nBin = prefitHisto.FindBin(xKnot)
                #            yKnot = prefitHisto.GetBinContent(nBin)
                #            print(yKnot)
                #            prePreInitPars.append([f'yKnot{nKnot}', yKnot, yKnot - (yKnot/100)*prefitFunc['bounds'], 
                #                                   yKnot + (yKnot/100)*prefitFunc['bounds']])
                #    else: 
                #        print('PREFITTINGGGG')
                #        prePreInitPars = [(prefitFunc[f'p{iPar}'][0], prefitFunc[f'p{iPar}'][1], prefitFunc[f'p{iPar}'][2], 
                #                           prefitFunc[f'p{iPar}'][3]) for iPar in range(prefitFunc['npars'])]   
                #    print(prefitFunc['funcname'])
                #    print(prePreInitPars)
                #    prePreFitters[-1].Add(prefitFunc['funcname'], prePreInitPars)
                #    prePreFitters[-1].Fit()
                #    cPrePrefit = TCanvas(f'cPrePrefit_{func["funcname"]}', '', 600, 600)
                #    prePreFitters[-1].Draw(cPrePrefit)
                #    print('DRAWN FUNCTION')
                #    oFile.cd(fitcf['fitname'])
                #    cPrePrefit.Write()
                #    prePrefitRes = prePreFitters[-1].GetFunction()
                #    nKnots = int(int(func['npars'])/2)
                #    if('spline3' in prefitFunc['funcname']):
                #        for iPar in range(nKnots):
                #            prefitInitPars.append([f'xKnot{iPar}', prePrefitRes.GetParameter(iPar), 
                #                             prePrefitRes.GetParameter(iPar), prePrefitRes.GetParameter(iPar)])            
                #        for iPar in range(nKnots):
                #            prefitInitPars.append([f'yKnot{iPar}', prePrefitRes.GetParameter(iPar + nKnots), 
                #                             prePrefitRes.GetParameter(iPar + nKnots) - (prePrefitRes.GetParameter(iPar + nKnots)/100)*prefitFunc['bounds'], 
                #                             prePrefitRes.GetParameter(iPar + nKnots) + (prePrefitRes.GetParameter(iPar + nKnots)/100)*prefitFunc['bounds']])
                #        print('PREFIT SPLINE')
                #    else: 
                #        prefitInitPars = [(prefitFunc[f'p{iPar}'][0], prefitFunc[f'p{iPar}'][1], prefitFunc[f'p{iPar}'][2], 
                #                           prefitFunc[f'p{iPar}'][3]) for iPar in range(prefitFunc['npars'])]    
                #    preFitters[-1].Add(prefitFunc['funcname'], prefitInitPars)
                #else:

                    #for iPar in range(func['npars']):
                    #    if(prePrefitRes.GetParameter(iPar) >= 0):
                    #        prefitInitPars.append([prefitFunc[f'p{iPar}'][0], prePrefitRes.GetParameter(iPar), 
                    #             prePrefitRes.GetParameter(iPar), prePrefitRes.GetParameter(iPar)])

                    #    else:
                    #        prefitInitPars.append([prefitFunc[f'p{iPar}'][0], prePrefitRes.GetParameter(iPar), 
                    #                         prePrefitRes.GetParameter(iPar), prePrefitRes.GetParameter(iPar)])      