'''
Script to compute the purities of correlated particles.
The output file is RawPurities_suffix.root

Usage:
python3 ComputePurities.py cfg.yml

'''
import os
import sys
import argparse
import yaml
import numpy as np

from ROOT import TFile, TH1F, gInterpreter, TH1D, TCanvas, gPad
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}/fempy/MassFitter.hxx"')
from ROOT import MassFitter
from fempy import logger as log
from fempy.utils.io import Load
from fempy.utils.analysis import WeightedAverage
from fempy.utils import DivideCanvas

parser = argparse.ArgumentParser()
parser.add_argument('cfg', help='Configuration file')
parser.add_argument('--debug', default=False, action='store_true', help='Run in debug mode (verbose)')
args = parser.parse_args()

if args.debug:
    log.setLevel(1)

# Load yaml file
with open(args.cfg, "r") as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        log.critical('Yaml configuration could not be loaded. Is it properly formatted?')

# Define the output file
oFileBaseName = 'Purities'
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'
oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)

runSuffix = cfg['runsuffix']

# purities of particles with MC method
if cfg['infilemc'] != '' and cfg['infilemc'] is not None:
    
    inFileMC = TFile(cfg['infilemc'])
    
    for partNum, MCpart in enumerate(cfg['MCpurityparts']):
        
        folderName = f'p{partNum}'
        
        # load MC histos relevant for purity
        hIdentPartPt = Load(inFileMC, MCpart['idpartfile'])
        hCorrPartPt = Load(inFileMC, MCpart['corridpartfile'])
        
        # compute purity for each pT bin
        hCorrPartPt.Divide(hIdentPartPt)
        hPurity = hCorrPartPt.Clone("hPurity")
        hPurity.SetTitle(';#it{p}_{T} (GeV/#it{c});Purity')
        hPurity.SetName(f'hPurity_{folderName}')
        
        # save purity with its error for each bin
        for iBin in range(hPurity.GetNbinsX()):
            nTotPart = hIdentPartPt.GetBinContent(iBin + 1)
            if(nTotPart != 0):
                binPurity = hPurity.GetBinContent(iBin + 1) 
                purityBinErr = np.sqrt( (binPurity * (1 - binPurity)) / nTotPart )
                hPurity.SetBinError(iBin + 1, purityBinErr)
            else:
                hPurity.SetBinError(iBin + 1, 0)
                
        oFile.cd()
        hPurity.Write()
        
        for partNum, pairedPartFile in enumerate(MCpart['avgptpairs']):
        
            partNum = partNum + len(cfg['IMpurityparts']) + 1
            partPairFolder = f'{folderName}{partNum}'
            print(partPairFolder)
            oFile.mkdir(partPairFolder)
            oFile.cd(partPairFolder)

            # load histo to compute averaged pT            
            kStarPt = Load(inFileMC, pairedPartFile)
            
            # compute yield of pions for each pT bin in k* range [0, 200] MeV/c2
            hPt = kStarPt.ProjectionY(f'femto{folderName}', 1, kStarPt.GetXaxis().FindBin(0.2*0.9999))
            hPt.Scale(1. / hPt.Integral())
            hPt.SetTitle(';#it{p}_{T} (GeV/#it{c});Counts')
            hPt.SetName(f'hPt_0_200')
            hPt.Write()
            
            # compute weights and their errors
            hWeights = hPurity.Clone("hWeights")
            hWeights.SetTitle(';#it{p}_{T} (GeV/#it{c});Counts')
            hWeights.SetName('hWeights')
            if(hPurity.GetBinWidth(1) == hPt.GetBinWidth(1)):
                for iBin in range(hWeights.GetNbinsX()):
                    hWeights.SetBinContent(iBin + 1, hPt.GetBinContent(iBin + 1))
                    hWeights.SetBinError(iBin + 1, hPt.GetBinError(iBin + 1))
            else:
                print(f"\033[33mWarning\033[0m: Histos have different bin widths, weights computed by interpolation!")
                
                # In the case of different bin widths, the weights for pTs corresponding to 
                # hPurity bin centers are obtained from interpolation of the two closest bins
                # in hPt histogram
                
                for iBin in range(hWeights.GetNbinsX()):
                    
                    purityPt = hPurity.GetBinCenter(iBin + 1)
                    hWeights.SetBinContent(iBin + 1, hPt.Interpolate(purityPt))
                    
                    purityBin = hPt.FindBin(purityPt)
                    purityBinCenter = hPt.GetBinCenter(purityBin)
                    
                    if(purityPt <= purityBinCenter):
                        x0 = hPt.GetBinCenter(purityBin-1)
                        x1 = hPt.GetBinCenter(purityBin)
                        y0 = hPt.GetBinContent(purityBin-1)
                        y1 = hPt.GetBinContent(purityBin)
                        sy0 = hPt.GetBinError(purityBin-1)
                        sy1 = hPt.GetBinError(purityBin)
                        xFactor = (purityPt - x0) / (x1 - x0) 
                        sWeight = np.sqrt( (1 - xFactor)**2 * sy0**2 + (xFactor * sy1)**2 )
                        hWeights.SetBinError(iBin + 1, sWeight)
                    
                    if(purityPt > purityBinCenter):
                        x0 = hPt.GetBinCenter(purityBin)
                        x1 = hPt.GetBinCenter(purityBin + 1)
                        y0 = hPt.GetBinContent(purityBin)
                        y1 = hPt.GetBinContent(purityBin + 1)
                        sy0 = hPt.GetBinError(purityBin)
                        sy1 = hPt.GetBinError(purityBin + 1)
                        xFactor = (purityPt - x0) / (x1 - x0) 
                        sWeight = np.sqrt( (1 - xFactor)**2 * sy0**2 + (xFactor * sy1)**2 )
                        hWeights.SetBinError(iBin + 1, sWeight)
                        
            # save histo with average purity value for each k* bin
            hSESample = Load(inFileMC, MCpart['clonedSEhisto'])
            hAvgPurity = hSESample.Clone('hAvgPurity')
            hAvgPurity.Reset('ICESM')
            hAvgPurity.SetTitle(';#it{k}^{*} (GeV/#it{c});Purity')
            avgPurity, avgPurityErr = WeightedAverage(hPurity, hWeights)
            for iBin in range(hAvgPurity.GetNbinsX()):
                hAvgPurity.SetBinContent(iBin, avgPurity) 
                hAvgPurity.SetBinError(iBin, avgPurityErr) 
            hAvgPurity.Write()
            hWeights.Write()

# purities of particles with IM fit method
if cfg['infiledata'] != '' and cfg['infiledata'] is not None:
    
    inFileData = TFile(cfg['infiledata'])
    IMPartList = cfg['IMpurityparts']
    
    for partNum, IMPart in enumerate(IMPartList):
        
        if(IMPart['nsigma']):
            nSigma = IMPart['nsigma']
        if(IMPart['intrange']):
            lowIntEdge = IMPart['PDGmass'] - IMPart['intrange']
            uppIntEdge = IMPart['PDGmass'] + IMPart['intrange']
        ptMins = IMPart['ptMins']
        ptMaxs = IMPart['ptMaxs']
        sgnFuncName = IMPart['sgnfuncname']
        bkgFuncName = IMPart['bkgfuncname']
        fitRange = IMPart['fitrange']
        splineYLimits = IMPart['splineylimits']
        fitSettingList = IMPart['fitsettings']
        fitDrawOpts = IMPart['drawopt']
        purities = []
        puritiesUncs = []
        hAvgPurity = TH1D('hAvgPurity', ';#it{p}_{T} (GeV/#it{c});Width (GeV/#it{c})', len(ptMins), 
                                 np.array(ptMins + [ptMaxs[-1]], 'd'))

        for iPartNum, iPartInvMassFile in enumerate(cfg['IMpurityparts'][partNum]['invmassfilepaths']):
        
            folderName = iPartNum + len(cfg['MCpurityparts'])
            oFile.mkdir(f'p{folderName}')

            # load TH2D with invariant mass spectra
            hInvMassPt = Load(inFileData, iPartInvMassFile)

            # create useful histos
            hInvMassWidthVsPt = TH1D('hInvMassWidthVsPt', ';#it{p}_{T} (GeV/#it{c});Width (GeV/#it{c})', len(ptMins), 
                                     np.array(ptMins + [ptMaxs[-1]], 'd'))
            hInvMassMeanVsPt = TH1D('hInvMassMeanVsPt', ';#it{p}_{T} (GeV/#it{c});Mean (GeV/#it{c})', len(ptMins), 
                                    np.array(ptMins + [ptMaxs[-1]], 'd'))
            hPurityVsPtAllUnc = TH1D('hPurityVsPtAllUnc', ';#it{p}_{T} (GeV/#it{c});Purity', len(ptMins), 
                               np.array(ptMins + [ptMaxs[-1]], 'd'))
            hPurityVsPtSystUnc = TH1D('hPurityVsPtSystUnc', ';#it{p}_{T} (GeV/#it{c});Purity', len(ptMins), 
                               np.array(ptMins + [ptMaxs[-1]], 'd'))
            hPurityVsPtStatUnc = TH1D('hPurityVsPtStatUnc', ';#it{p}_{T} (GeV/#it{c});Purity', len(ptMins), 
                               np.array(ptMins + [ptMaxs[-1]], 'd'))
            hChi2VsPt = TH1D('hChi2VsPt', ';#it{p}_{T} (GeV/#it{c});#chi2/ndf', len(ptMins), 
                             np.array(ptMins + [ptMaxs[-1]], 'd'))
            hSgnWindowChi2VsPt = TH1D('hSgnWindowChi2VsPt', ';#it{p}_{T} (GeV/#it{c});#chi2/ndf', len(ptMins), 
                             np.array(ptMins + [ptMaxs[-1]], 'd'))
            hYieldVsPt = TH1D('hYieldVsPt', ';#it{p}_{T} (GeV/#it{c});Yield', len(ptMins), 
                               np.array(ptMins + [ptMaxs[-1]], 'd'))
            hSgnVsPt = TH1D('hSgnVsPt', ';#it{p}_{T} (GeV/#it{c});Counts', len(ptMins), 
                             np.array(ptMins + [ptMaxs[-1]], 'd'))
            hBkgVsPt = TH1D('hBkgVsPt', ';#it{p}_{T} (GeV/#it{c});Counts', len(ptMins), 
                             np.array(ptMins + [ptMaxs[-1]], 'd'))
            hRatioSgnBkgVsPt = TH1D('hRatioSgnBkgVsPt', ';#it{p}_{T} (GeV/#it{c});S/B', len(ptMins), 
                             np.array(ptMins + [ptMaxs[-1]], 'd'))
            hRelUncVsPt = TH1D('hRelUncVsPt', ';#it{p}_{T} (GeV/#it{c});#sigma_P/P', len(ptMins), 
                             np.array(ptMins + [ptMaxs[-1]], 'd'))

            # create list of fit configuration
            fitters = []
            #cInvMass = TCanvas('cInvMass', '', 1200, 800)
            cInvMass = TCanvas('cInvMass', '', 600, 600)
            DivideCanvas(cInvMass, len(ptMins))

            for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
                cInvMass.cd(iPt+1)

                # project the TH2D to obtain the invariant mass spectra to fit 
                lowBin = hInvMassPt.GetXaxis().FindBin(ptMin*1.0001) 
                uppBin = hInvMassPt.GetXaxis().FindBin(ptMax*0.9999) 
                hMassProj = hInvMassPt.ProjectionY(f'hMass_{ptMin}_{ptMax}', lowBin, uppBin)
                hMassProj.SetTitle(f'{ptMin:.1f} < #it{{p}}_{{T}} < {ptMax:.1f} GeV/#it{{c}}')

                # execute fit with settings from config file
                fitters.append(MassFitter(hMassProj, sgnFuncName, bkgFuncName, fitRange[0], fitRange[1]))
                fitPars = []
                if('spline' in bkgFuncName):
                    prefitPars = [] 
                for parName, parValue in fitSettingList[iPt].items():
                    if(parName == 'xknots'):
                        for nKnot, xKnot in enumerate(parValue):
                            nBin = hMassProj.FindBin(xKnot)
                            xCoord = hMassProj.GetBinCenter(nBin)
                            fitPars.append([f'xKnot{nKnot}', xCoord, xCoord, xCoord])
                            prefitPars.append([f'xKnot{nKnot}', xCoord, xCoord, xCoord])
                        for nKnot, xKnot in enumerate(parValue):
                            nBin = hMassProj.FindBin(xKnot)
                            yCoord = hMassProj.GetBinContent(nBin)
                            fitPars.append([f'yKnot{nKnot}', yCoord, yCoord - (yCoord/100)*splineYLimits, yCoord + (yCoord/100)*splineYLimits])
                            prefitPars.append([f'yKnot{nKnot}', yCoord, yCoord - (yCoord/100)*splineYLimits, yCoord + (yCoord/100)*splineYLimits])
                    else:
                        fitPars.append([parValue[0], parValue[1], parValue[2], parValue[3]]) 
                fitters[-1].Add(fitPars)
                if('spline' in bkgFuncName):
                    fitters[-1].AddPrefit(prefitPars)  
                    fitters[-1].Prefit()            
                fitters[-1].Fit()
                if(IMPart['nsigma']):
                    fitters[-1].SetIntegrationEdges(nSigma)
                if(IMPart['intrange']):
                    fitters[-1].SetIntegrationEdges(lowIntEdge, uppIntEdge)
                fitters[-1].Draw(gPad, fitDrawOpts, "data_minus_bkg")

                ### fill histos with fit information
                hInvMassWidthVsPt.SetBinContent(iPt+1, fitters[-1].GetWidth())
                hInvMassWidthVsPt.SetBinError(iPt+1, fitters[-1].GetWidthUnc())
                hInvMassMeanVsPt.SetBinContent(iPt+1, fitters[-1].GetMean())
                hInvMassMeanVsPt.SetBinError(iPt+1, fitters[-1].GetMeanUnc())
                hChi2VsPt.SetBinContent(iPt+1, fitters[-1].GetChi2Ndf())
                hSgnWindowChi2VsPt.SetBinContent(iPt+1, fitters[-1].GetSgnWindowChi2Ndf())

                sgnInt = fitters[-1].GetSignal("data_minus_bkg")
                sgnIntUnc = fitters[-1].GetSignalUnc("data_minus_bkg")
                bkgInt = fitters[-1].GetBackground()
                bkgIntUnc = fitters[-1].GetBackgroundUnc()
                hSgnVsPt.SetBinContent(iPt+1, sgnInt)
                hSgnVsPt.SetBinError(iPt+1, sgnIntUnc)
                hBkgVsPt.SetBinContent(iPt+1, bkgInt)
                hBkgVsPt.SetBinError(iPt+1, bkgIntUnc)
                hRatioSgnBkgVsPt.SetBinContent(iPt+1, sgnInt / bkgInt)
                hRatioSgnBkgVsPt.SetBinError(iPt+1, np.sqrt( (sgnIntUnc / bkgInt)**2 + 
                                                             ( (sgnInt * bkgIntUnc) / bkgIntUnc**2 )**2 ))
                
                ### compute purity and error
                purity = ( fitters[-1].GetPurity("data_minus_bkg") + 
                          fitters[-1].GetPurity("data_minus_bkg") ) / 2
                purities.append(purity)
                purityUnc = ( abs( fitters[-1].GetPurity("data_minus_bkg") - 
                          fitters[-1].GetPurityPrefit("data_minus_bkg") ) ) / 2
                puritiesUncs.append(purityUnc)
                hPurityVsPtSystUnc.SetBinContent(iPt+1, purity)
                hPurityVsPtSystUnc.SetBinError(iPt+1, fitters[-1].GetPuritySystUnc("data_minus_bkg"))
                hPurityVsPtStatUnc.SetBinContent(iPt+1, purity)
                hPurityVsPtStatUnc.SetBinError(iPt+1, fitters[-1].GetPurityStatUnc("data_minus_bkg"))
                hRelUncVsPt.SetBinContent(iPt+1, fitters[-1].GetPurityAllUnc("data_minus_bkg")/purity)
                hRelUncVsPt.SetBinError(iPt+1, 0.00000000000000)

                for iBin in range(hPurityVsPtAllUnc.GetNbinsX()):
                    hPurityVsPtAllUnc.SetBinContent(iPt+1, purity)
                    hPurityVsPtAllUnc.SetBinError(iPt+1, fitters[-1].GetPurityAllUnc("data_minus_bkg"))
                
                fitters[-1].SetIntegrationEdges(5)
                sgnYield = fitters[-1].GetSignal("data_minus_bkg")
                sgnYieldUnc = fitters[-1].GetSignalUnc("data_minus_bkg")
                bkgYield = fitters[-1].GetBackground()
                bkgYieldUnc = fitters[-1].GetBackgroundUnc()
                hYieldVsPt.SetBinContent(iPt+1, sgnYield + bkgYield)
                hYieldVsPt.SetBinError(iPt+1, np.sqrt( sgnYieldUnc**2 + bkgYieldUnc**2) )
                
                oFile.cd(f'p{folderName}')
                hMassProj.Write()
            
            # load histo to compute averaged pT
            for iPair in range(2):            
                kStarPt = Load(inFileData, IMPart['avgptpairs'][iPair + 2*iPartNum])
            
                # compute yield of pions for each pT bin in k* range [0, 200] MeV/c2
                hPt = kStarPt.ProjectionY(f'femto{folderName}', 1, kStarPt.GetXaxis().FindBin(0.2*0.9999))
                hPt.Scale(1. / hPt.Integral())
                hPt.SetTitle(';#it{p}_{T} (GeV/#it{c});Counts')
                hPt.SetName(f'hPt_0_200')
                hPt.Write()

                # compute weights and their errors
                hWeights = hPurityVsPtAllUnc.Clone("hWeights")
                hWeights.SetTitle(';#it{p}_{T} (GeV/#it{c});Counts')
                hWeights.SetName('hWeights')
                if(hPurityVsPtAllUnc.GetBinWidth(1) == hPt.GetBinWidth(1)):
                    for iBin in range(hWeights.GetNbinsX()):
                        hWeights.SetBinContent(iBin + 1, hPt.GetBinContent(iBin + 1))
                        hWeights.SetBinError(iBin + 1, hPt.GetBinError(iBin + 1))
                
                elif((hPurityVsPtAllUnc.GetBinWidth(1) % hPt.GetBinWidth(1)) < 0.09):
                    print(f"\033[33mWarning\033[0m: Histos have different bin widths but share larger binwidth, weights histogram rebinned by {hPurityVsPtAllUnc.GetBinWidth(1) % hPt.GetBinWidth(1)}!")
                    print(hWeights.GetNbinsX())
                    hPt.Rebin(int(hPurityVsPtAllUnc.GetBinWidth(1) / hPt.GetBinWidth(1)))
                    for iBin in range(hWeights.GetNbinsX()):
                        hWeights.SetBinContent(iBin + 1, hPt.GetBinContent(iBin + 1))
                        hWeights.SetBinError(iBin + 1, hPt.GetBinError(iBin + 1))
                else:
                    print(f"\033[33mWarning\033[0m: Histos have different bin widths and don't share any binwidths, weights computed by interpolation!")

                    # In the case of different bin widths, the weights for pTs corresponding to 
                    # hPurityVsPtAllUnc bin centers are obtained from interpolation of the two closest bins
                    # in hPt histogram

                    for iBin in range(hWeights.GetNbinsX()):

                        purityPt = hPurityVsPtAllUnc.GetBinCenter(iBin + 1)
                        hWeights.SetBinContent(iBin + 1, hPt.Interpolate(purityPt))

                        purityBin = hPt.FindBin(purityPt)
                        purityBinCenter = hPt.GetBinCenter(purityBin)

                        if(purityPt <= purityBinCenter):
                            x0 = hPt.GetBinCenter(purityBin-1)
                            x1 = hPt.GetBinCenter(purityBin)
                            y0 = hPt.GetBinContent(purityBin-1)
                            y1 = hPt.GetBinContent(purityBin)
                            sy0 = hPt.GetBinError(purityBin-1)
                            sy1 = hPt.GetBinError(purityBin)
                            xFactor = (purityPt - x0) / (x1 - x0) 
                            sWeight = np.sqrt( (1 - xFactor)**2 * sy0**2 + (xFactor * sy1)**2 )
                            hWeights.SetBinError(iBin + 1, sWeight)

                        if(purityPt > purityBinCenter):
                            x0 = hPt.GetBinCenter(purityBin)
                            x1 = hPt.GetBinCenter(purityBin + 1)
                            y0 = hPt.GetBinContent(purityBin)
                            y1 = hPt.GetBinContent(purityBin + 1)
                            sy0 = hPt.GetBinError(purityBin)
                            sy1 = hPt.GetBinError(purityBin + 1)
                            xFactor = (purityPt - x0) / (x1 - x0) 
                            sWeight = np.sqrt( (1 - xFactor)**2 * sy0**2 + (xFactor * sy1)**2 )
                            hWeights.SetBinError(iBin + 1, sWeight)

                # save histo with average purity value for each k* bin
                hSESample = Load(inFileData, IMPart['clonedSEhistos'][iPair + iPartNum])
                hAvgPurity = hSESample.Clone('hAvgPurity')
                hAvgPurity.Reset('ICESM')
                hAvgPurity.SetTitle(';#it{k}^{*} (GeV/#it{c});Purity')
                avgPurity, avgPurityErr = WeightedAverage(hPurityVsPtAllUnc, hWeights)
                for iBin in range(hAvgPurity.GetNbinsX()):
                    hAvgPurity.SetBinContent(iBin, avgPurity) 
                    hAvgPurity.SetBinError(iBin, avgPurityErr) 
                folderNamePair = str(iPartNum + len(cfg['MCpurityparts'])) + str(iPair)
                oFile.mkdir(f'p{folderNamePair}')
                oFile.cd(f'p{folderNamePair}')
                hAvgPurity.Write()
                hWeights.Write()

            oFile.cd(f'p{folderName}')
            hInvMassWidthVsPt.Write()
            hInvMassMeanVsPt.Write()
            hPurityVsPtAllUnc.Write()
            hPurityVsPtSystUnc.Write()
            hPurityVsPtStatUnc.Write()
            hRelUncVsPt.Write()
            hChi2VsPt.Write()
            hSgnWindowChi2VsPt.Write()
            hYieldVsPt.Write()
            hSgnVsPt.Write()
            hBkgVsPt.Write()
            hRatioSgnBkgVsPt.Write()
            cInvMass.Write()
        
oFile.Close()
print(f'output saved in {oFileName}')

