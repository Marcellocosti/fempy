'''
Script to compute the raw correlation function.
The output file is RawCF_suffix.root

Usage:
python3 ComputeRawCF.py cfg.yml

'''
import os
import argparse
import yaml

from ROOT import TFile, TH2F, TH1D

from fempy import logger as log
from fempy.utils.io import Load, GetKeyNames

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

regions = ['sgn']
combs = ['p02', 'p03', 'p12', 'p13']

# Load input file with same- and mixed-event distributions
inFile = TFile(cfg['infile'])

# Define the output file
oFileBaseName = 'RawCF'
if cfg['suffix'] != '' and cfg['suffix'] is not None:
    oFileBaseName += f'_{cfg["suffix"]}'
oFileName = os.path.join(cfg['odir'], oFileBaseName + '.root')

# Open the output file
try:
    oFile = TFile.Open(oFileName, 'recreate')
except OSError:
    log.critical('The output file %s is not writable', oFileName)

hSE = {}
hME = {}
hMErew = {}
for comb in combs:
    iPart1 = int(comb[1])
    iPart2 = int(comb[2])
    oFile.mkdir(comb)
    oFile.cd(comb)
    
    hSE[comb] = {}
    hME[comb] = {}
    hMErew[comb] = {}
    for region in regions:
        runSuffix = cfg['runsuffix']
        if f'HMResults{runSuffix}' in GetKeyNames(inFile): # Make correlation functions from FemtoDream
            fdcomb = f'Particle{iPart1}_Particle{iPart2}'
            # The histograms are casted to TH1D with TH1::Copy to avoid NotImplementedError when computing hSE/hME
            hSE[comb][region] = TH1D()
            Load(inFile, f'HMResults{runSuffix}/HMResults{runSuffix}/{fdcomb}/SEDist_{fdcomb}').Copy(hSE[comb][region])

            # The histograms are casted to TH1D with TH1::Copy to avoid NotImplementedError when computing hSE/hME
            hME[comb][region] = TH1D()
            Load(inFile, f'HMResults{runSuffix}/HMResults{runSuffix}/{fdcomb}/MEDist_{fdcomb}').Copy(hME[comb][region])

            # No need to cast the these because the projection of TH2D and TH2F is always TH1D
            hSEmultk = Load(inFile, f'HMResults{runSuffix}/HMResults{runSuffix}/{fdcomb}/SEMultDist_{fdcomb}')
            hMEmultk = Load(inFile, f'HMResults{runSuffix}/HMResults{runSuffix}/{fdcomb}/MEMultDist_{fdcomb}')

        elif comb in GetKeyNames(inFile): # Make correlation functions from ALICE3 simulations
            # The histograms are casted to TH1D with TH1::Copy to avoid NotImplementedError when computing hSE/hME
            hSE[comb][region] = TH1D()
            Load(inFile, f'{comb}/hSE').Copy(hSE[comb][region])

            # The histograms are casted to TH1D with TH1::Copy to avoid NotImplementedError when computing hSE/hME
            hME[comb][region] = TH1D()
            Load(inFile, f'{comb}/hME').Copy(hME[comb][region])

            # Mult reweighting not implemented. Keep dummy histogram for now
            nbins = hSE[comb][region].GetNbinsX()
            xMin = hSE[comb][region].GetXaxis().GetXmin()
            xMax = hSE[comb][region].GetXaxis().GetXmax()
            hSEmultk = TH2F('hSEMult', '', nbins, xMin, xMax, 200, 0, 200)
            hMEmultk = TH2F('hMEMult', '', nbins, xMin, xMax, 200, 0, 200)

        nbins = hMEmultk.ProjectionX().GetNbinsX()
        hMEreweightk = TH1D("MErewdistr", "MErewdistr", nbins, hMEmultk.GetXaxis().GetXmin(), hMEmultk.GetXaxis().GetXmax())
        for iBin in range(hMEmultk.ProjectionY().GetNbinsX() + 2): # Loop over underflow, all bins, and overflow
            hSEbinmult = hSEmultk.ProjectionX(f'{comb}SEdistr', iBin, iBin)
            hMEbinmult = hMEmultk.ProjectionX(f'{comb}MEdistr', iBin, iBin)
            if(hMEbinmult.Integral() > 0):
                hMEreweightk.Add(hMEbinmult, hSEbinmult.Integral()/hMEbinmult.Integral())
        hMErew[comb][region] = hMEreweightk

# Sum pair and antipair
for comb in combs:
    hSE['p02_13'] = {}
    hSE['p03_12'] = {}
    hME['p02_13'] = {}
    hME['p03_12'] = {}
    hMErew['p02_13'] = {}
    hMErew['p03_12'] = {}

    for region in regions:
        hSE['p02_13'][region] = hSE['p02'][region] + hSE['p13'][region]
        hME['p02_13'][region] = hME['p02'][region] + hME['p13'][region]
        hMErew['p02_13'][region] = hMErew['p02'][region] + hMErew['p13'][region]
        
        hSE['p03_12'][region] = hSE['p03'][region] + hSE['p12'][region]
        hME['p03_12'][region] = hME['p03'][region] + hME['p12'][region]
        hMErew['p03_12'][region] = hMErew['p03'][region] + hMErew['p12'][region]

# Compute the CF and write to file
for comb in combs + ['p02_13', 'p03_12']:
    for region in regions:
        rebin = round(float(cfg['binwidth']) / (hSE[comb][region].GetBinWidth(1) * 1000))
        hSE[comb][region].Rebin(rebin)
        hME[comb][region].Rebin(rebin)
        hMErew[comb][region].Rebin(rebin)

        if cfg['norm'] is None:  # if not specified, normalize to the yields
            norm = hME[comb][region].GetEntries() / hSE[comb][region].GetEntries()
            normrew = hMErew[comb][region].GetEntries() / hSE[comb][region].GetEntries()
        else:
            lowEdge = hSE[comb][region].FindBin(cfg['norm'][0]*1.0001)
            uppEdge = hSE[comb][region].FindBin(cfg['norm'][1]*0.9999)
            norm = hME[comb][region].Integral(lowEdge, uppEdge) / hSE[comb][region].Integral(lowEdge, uppEdge)
            normrew = hMErew[comb][region].Integral(lowEdge, uppEdge) / hSE[comb][region].Integral(lowEdge, uppEdge)

        hSE[comb][region].Sumw2()
        hME[comb][region].Sumw2() # Just to trigger the same draw option as for hSE
        
        hCF = norm * hSE[comb][region] / hME[comb][region]
        hCFrew = normrew * hSE[comb][region] / hMErew[comb][region]

        oFile.mkdir(f'{comb}/{region}')
        oFile.cd(f'{comb}/{region}')

        hSE[comb][region].SetName('hSE')
        hSE[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')
        hSE[comb][region].Write()

        hME[comb][region].SetName('hME')
        hME[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')
        hME[comb][region].Write()
        
        hMErew[comb][region].SetName('hMErew')
        hMErew[comb][region].SetTitle(';#it{k}* (GeV/#it{c});Counts')

        hCF.SetName('hCF')
        hCF.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
        hCF.Write()
        
        hCFrew.SetName('hCFrew')
        hCFrew.SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
        hCFrew.Write()

oFile.Close()
print(f'output saved in {oFileName}')
