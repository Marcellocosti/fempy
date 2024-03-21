'''
Compute the invariant mass of the pairs from the relative momentum k*.
'''

import argparse
import yaml
    
from ROOT import TFile, TH1D

import os

from ROOT import TFile, TCanvas, gInterpreter, TF1, TDatabasePDG, TSpline3
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/CorrelationFitter.hxx"')
gInterpreter.ProcessLine(f'#include "{os.environ.get("FEMPY")}fempy/Fitter.hxx"')
from ROOT import Fitter

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

with open(args.cfg, 'r') as file:
    cfg = yaml.safe_load(file)

# Define the output file
oFileBaseName = 'InvMass'
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'
oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')
oFile = TFile(oFileName, 'recreate')

# Load Data AnalysisResults
inFileData = TFile(cfg['data']['file'])
hSEData = Load(inFileData, cfg['data']['se_path'])
hMEData = Load(inFileData, cfg['data']['me_path'])
hSEData = ChangeUnits(hSEData, 1000)
hMEData = ChangeUnits(hMEData, 1000)

# Load MC AnalysisResults
inFileMC = TFile(cfg['mc']['file'])
hSEMC = Load(inFileMC, cfg['mc']['se_path'])
hMEMC = Load(inFileMC, cfg['mc']['me_path'])
hSEMC = ChangeUnits(hSEMC, 1000)
hMEMC = ChangeUnits(hMEMC, 1000)

oFile.cd()

# Normalize mixed event to same event
firstBin = hSEData.GetXaxis().FindBin(cfg['norm_range'][0] * 1.001)
lastBin = hSEData.GetXaxis().FindBin(cfg['norm_range'][1] * 0.999)
hMEData.Scale(hSEData.Integral(firstBin, lastBin) / hMEData.Integral(firstBin, lastBin))
hMEMC.Scale(hSEMC.Integral(firstBin, lastBin) / hMEMC.Integral(firstBin, lastBin))

# Bin error as sqrt of bin content
hSEData.Sumw2()
hMEData.Sumw2()
hSEMC.Sumw2()
hMEMC.Sumw2()

# Subtraction on data
hSubSEME_Data = hSEData - hMEData
hSubSEME_Data.Scale(1/hSubSEME_Data.Integral())
hSubSEME_Data.Sumw2()
hSubSEME_Data.Write('SE-ME_Data')

# Subtraction on MC
hSubSEME_MC = hSEMC - hMEMC
hSubSEME_MC.Scale(1/hSubSEME_MC.Integral())
hSubSEME_MC.Sumw2()
hSubSEME_MC.Write('SE-ME_MC')

# Write SE and ME histograms to file
#hSEMC.Scale(hSEData.Integral(firstBin, lastBin) / hSEMC.Integral(firstBin, lastBin))
#hMEMC.Scale(hSEData.Integral(firstBin, lastBin) / hMEMC.Integral(firstBin, lastBin))
hSEData.Write('hSEData')
hMEData.Write('hMEData')
hSEMC.Write('hSEMC')
hMEMC.Write('hMEMC')

# Fit the subtraction of SE and ME in MC, excluding the peak
MCFit = Fitter(hSubSEME_MC, 0, 400)
initParsMCFitPol4 = [("p0", 0, -1, 1), ("p1", 0, -1, 1), ("p2", 0, -1, 1), ("p3", 0, -1, 1), ("p4", 0, -1, 1)]
MCFit.Add("pol4", initParsMCFitPol4)
MCFit.Fit(150, 300)
cFit = TCanvas('cFit_SE-ME_MC', '', 600, 600)
MCFit.Draw(cFit)
cFit.Write()

# Subtraction between data and MC fit
MCFitFunc = MCFit.GetFunction() 
#hInvMass = hSubSEME_Data.Clone('ICESM')
#for iBin in range(hInvMass.GetNbinsX()):
hInvMass = TH1D("Data-MCFit", "Data-MCFit", 100, 0, 400)
for iBin in range(100):
    hInvMass.SetBinContent(iBin + 1, hSubSEME_Data.GetBinContent(iBin + 1) - 
                                     MCFitFunc.Eval(hSubSEME_Data.GetBinCenter(iBin + 1)))
hInvMass.Write('Data-MCFit')

# Fit the subtraction describing the peak with a Breit Wigner
MCFitWithBW = Fitter(hSubSEME_MC, 0, 400)
#initParsMCFitPol4BW = [("p0", 0, -1, 1), ("p1", 0, -1, 1), ("p2", 0, -1, 1), ("p3", 0, -1, 1), ("p4", 0, -1, 1)]
initParsMCFitPol4BW = [("p0", 5.78801e-05, 1, -1), ("p1", -2.87197e-05, 1, -1), 
                       ("p2", 1.06536e-06, 1, -1), ("p3", -4.45476e-09, 1, -1), 
                       ("p4", 5.09726e-12, 1, -1)]
MCFitWithBW.Add("pol4", initParsMCFitPol4BW)
initParsMCFitWithBW = [("yield", 0.35, 0.3, 0.5), ("mean", 204, 202, 206), ("width", 18, 14, 20)]
MCFitWithBW.Add("gaus", initParsMCFitWithBW)
MCFitWithBW.Fit()
cFit = TCanvas('cFit_SE-ME_MC_with_bw', '', 600, 600)
MCFitWithBW.Draw(cFit)
cFit.Write()

# Subtraction between data and MC fit including the Breit Wigner
MCFitWithBWFunc = MCFitWithBW.GetFunction() 
hInvMassBW = TH1D("Data-MCFit_bw", "Data-MCFit_bw", 100, 0, 400)
for iBin in range(100):
    hInvMassBW.SetBinContent(iBin + 1, hSubSEME_Data.GetBinContent(iBin + 1) - 
                                     MCFitWithBWFunc.Eval(hSubSEME_Data.GetBinCenter(iBin + 1)))
hInvMassBW.Write('Data-MCFit_bw')

# Subtraction between data and MC spline
splineMCHisto = TSpline3(hSubSEME_MC)
hInvMassSpline = TH1D("Data-MC_Spline", "Data-MC_Spline", 100, 0, 400)
for iBin in range(100):
    hInvMassSpline.SetBinContent(iBin + 1, hSubSEME_Data.GetBinContent(iBin + 1) - 
                                     splineMCHisto.Eval(hSubSEME_Data.GetBinCenter(iBin + 1)))
splineMCHisto.Write()
hInvMassSpline.Write('Data-MC_Spline')

############ INCLUDING PYTHIA SIMULATIONS ##############
inFileXi1530 = TFile("/home/mdicostanzo/an/LPi/Simulation/outputs/NewSimStandalonePythia_MERGED.root")
hXi1530 = Load(inFileXi1530, "_2113122_SMEARED/hSE_2113122_3324_SMEARED")
thFistPrimYieldXi1530 = 0.027351
# Charge factor
BRrewYieldXi1530 = thFistPrimYieldXi1530 * (1/3)
hXi1530.Scale(BRrewYieldXi1530/hXi1530.Integral())

inFileLambda1520 = TFile("/home/mdicostanzo/an/LPi/DPG/SimLambda1520_KinemXi1530_42.root")
hLambda1520 = Load(inFileLambda1520, "3122211_-3122-211_smeared")
thFistPrimYieldLambda1520 = 0.0328154
BRrewYieldLambda1520 = 0.0328154
hLambda1520.Scale(BRrewYieldLambda1520/hLambda1520.Integral())

inFileSigma1385Indir = TFile("/home/mdicostanzo/an/LPi/Simulation/outputs/SimSigma1385BkgNew.root")
hSigma1385Indir = Load(inFileSigma1385Indir, "hSEPiPlusLambda")
# BR Sigma0 and charge factor
thFistPrimYieldSigma1385 = 0.0578607
BRrewYieldSigma1385 = 0.0578607 * (11.7/100) * (1/3)
hSigma1385Indir.Scale(BRrewYieldSigma1385/hSigma1385Indir.Integral())

print('Relative abundances in the correlation: ' + str(BRrewYieldLambda1520) + ' ' + str(BRrewYieldXi1530) + ' ' + str(BRrewYieldSigma1385))

hReweightedResonances = hXi1530.Clone('ICESM')
for iBin in range(hReweightedResonances.GetNbinsX()):
    hReweightedResonances.SetBinContent(iBin + 1, hSigma1385Indir.GetBinContent(iBin + 1) + hLambda1520.GetBinContent(iBin + 1) + hXi1530.GetBinContent(iBin + 1))
oFile.cd()
hReweightedResonances.Write('Reweighted_Resonances')
hReweightedResonances.Scale(hInvMassBW.Integral()/hReweightedResonances.Integral())
hReweightedResonances.Write('Reweighted_Resonances_Rescaled')

# Computing the weight to add the distribution to the MC
hMCWithResonances = hReweightedResonances.Clone('ICESM')
print(hReweightedResonances.GetNbinsX())
print(hSubSEME_MC.GetNbinsX())
for iBin in range(hReweightedResonances.GetNbinsX()):
    hMCWithResonances.SetBinContent(iBin + 1, 3*hReweightedResonances.GetBinContent(iBin + 1) + hSubSEME_MC.GetBinContent(iBin + 1))
hMCWithResonances.Scale(1/hMCWithResonances.Integral())
oFile.cd()
hMCWithResonances.Write('MC_with_resonances')

splineMCHistoWithRes = TSpline3(hMCWithResonances)
hInvMassWithResSpline = TH1D("Data-MC_SplineRes", "Data-MC_SplineRes", 100, 0, 400)
for iBin in range(100):
    hInvMassWithResSpline.SetBinContent(iBin + 1, hSubSEME_Data.GetBinContent(iBin + 1) - 
                                        splineMCHistoWithRes.Eval(hSubSEME_Data.GetBinCenter(iBin + 1)))
splineMCHistoWithRes.Write()
hInvMassWithResSpline.Write('Data-MC_SplineRes')

oFile.Close()
print(f'Output saved in: {oFileName}')