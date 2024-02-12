'''
Utility functions and helpers.
'''

import sys

import yaml

from ROOT import TDatabasePDG, TGraphErrors, TH1F
import math
from math import sqrt
import fempy

class Pair:
    def __init__(self, pair):
        self.name = pair
        if pair == 'DstarPi' or pair == 'DstarK':
            self.max_kstar = 3000
            self.binwidths = [1, 2, 4, 5, 10, 20, 40, 50]
            self.mass_regions = ['sgn', 'sbr']
            self.hpdg = 413
            self.lpdg = 211
            self.norm_range = [1500, 2000]
            self.heavy_mass_label = '#it{M}(K#pi#pi) #minus #it{M}(K#pi) (GeV/#it{c}^{2})'
            self.cfg_sidebands = None

        elif pair == 'DPi':
            self.max_kstar = 3000
            self.binwidths = [1, 2, 4, 5, 10, 20, 40, 50]
            self.mass_regions = ['sgn', 'sbl', 'sbr']
            self.hpdg = 411
            self.lpdg = 211
            self.norm_range = [1000, 1500]
            self.heavy_mass_label = '#it{M}(K#pi#pi) (GeV/#it{c})^{2}'
            self.cfg_sidebands = None
        elif pair == 'DK':
            self.max_kstar = 3000
            self.binwidths = [1, 2, 4, 5, 10, 20, 40, 50]
            self.mass_regions = ['sgn', 'sbl', 'sbr']
            self.hpdg = 411
            self.lpdg = 321
            self.norm_range = [1500, 2000]
            self.heavy_mass_label = '#it{M}(K#pi#pi) (GeV/#it{c})^{2}'
            self.cfg_sidebands = None
        else:
            print("Error: pair not implemented. Exit!")
            sys.exit()

        # with open(cfg_file, "r") as stream:
        #     try:
        #         cfg = yaml.safe_load(stream)
        #     except yaml.YAMLError as exc:
        #         print(exc)
        #         sys.exit()

        #     self.name = 'DstarPi'
        #     print(cfg.get('sidebands'))
        #     self.max_kstar = cfg['max_kstar']
        #     self.binwidths = cfg['binwidths']
        #     self.mass_regions = cfg['mass_regions']
        #     self.hpdg = cfg['hpdg']
        #     self.lpdg = cfg['lpdg']
        #     self.norm_range = cfg['norm_range']
        #     self.heavy_mass_label = cfg['heavy_mass_label']
        #     self.cfg_sidebands = cfg.get('sidebands')


def is_mass_selected(mass, pt, pdg=413, selection='sgn', nsigma_mass=2., nsigma_offset=5., sideband_width=0.2, lower_Dstar_removal=1.992, upper_Dstar_removal=2.028) -> bool:
    '''
    function to perform the pT dependent selection of D+, D* candidates in
    different regions of the invariant mass distribution. Based on
    AliAnalysisTaskCharmingFemto::MassSelection.
    '''

    if selection == 'any':
        return True

    # mass shift observed in all Run2 data samples for all D-meson species
    mass_mean = TDatabasePDG.Instance().GetParticle(pdg).Mass() + 0.0025
    if pdg == 411:
        massWidth = 0.006758 + pt * 0.0005124
    elif pdg == 413:
        mDstarPDG = TDatabasePDG.Instance().GetParticle(413).Mass()
        mD0PDG = TDatabasePDG.Instance().GetParticle(421).Mass()
        mass_mean = mDstarPDG-mD0PDG  # no extra mass shift because it is deltamass
        massWidth = 0.00124673 - pt * 0.000340426 + pt * pt * 4.40729e-05
        if pt > 4 and pt < 5:
            massWidth = 0.00104329 - 0.000113275 * pt
        elif pt >= 5:
            massWidth = 0.000519861 - 8.58874e-06 * pt

    # select D mesons mass window
    if selection == 'sgn':
        lower_selection = mass_mean - nsigma_mass * massWidth
        upper_selection = mass_mean + nsigma_mass * massWidth
    elif selection == 'sbl':
        lower_selection = mass_mean - nsigma_offset * massWidth - sideband_width
        upper_selection = mass_mean - nsigma_offset * massWidth
    elif selection == 'sbr':
        lower_selection = mass_mean + nsigma_offset * massWidth
        upper_selection = mass_mean + nsigma_offset * massWidth + sideband_width

        if pdg == 411 and mass > lower_Dstar_removal and mass < upper_Dstar_removal:
            return False

    if mass > lower_selection and mass < upper_selection:
        return True

    return False

def GetNormFactor(se, me, fromVal, toVal):
    firstBin = se.FindBin(fromVal*1.0001)
    lastBin = se.FindBin(toVal*0.9999)

    return me.Integral(firstBin, lastBin) / se.Integral(firstBin, lastBin)


def ComputeBinBrackets(hist, name='gBrackets'):
    gBrackets = TGraphErrors(1)
    gBrackets.SetName(name)
    for iBin in range(hist.GetNbinsX()):
        gBrackets.SetPoint(iBin, hist.GetBinCenter(iBin+1), hist.GetBinContent(iBin+1))
        gBrackets.SetPointError(iBin, hist.GetBinWidth(iBin+1)/2, 0)
    return gBrackets

def WeightedAverage(hQuantity, hWeight):

    totWeight = hWeight.Integral()
    avgQuantity = 0
    avgQuantityErr = 0
    for iBin in range(hQuantity.GetNbinsX()):
        
        binQuantity = hQuantity.GetBinContent(iBin + 1)
        binWeight = hWeight.GetBinContent(iBin + 1)
        avgQuantity += binQuantity * binWeight
        
    avgQuantity = avgQuantity / totWeight

    for iBin in range(hQuantity.GetNbinsX()):
        
        binQuantity = hQuantity.GetBinContent(iBin + 1)
        binWeight = hWeight.GetBinContent(iBin + 1)

        binQuantityErr = hQuantity.GetBinError(iBin + 1)
        binWeightErr = hWeight.GetBinError(iBin + 1)
        weightErr = ( binQuantity / totWeight - ( avgQuantity / (totWeight**2) ) ) * binWeightErr 
        quantityErr = binWeight * binQuantityErr
        avgQuantityErr += weightErr**2 + quantityErr**2 
    
    avgQuantityErr = math.sqrt(avgQuantityErr)        

    return avgQuantity, avgQuantityErr

def ChangeUnits(hist, multiplier):
    '''
    Only for histogram with constant binwidth!
    '''
    nbins = hist.GetNbinsX()
    lowEdge = hist.GetBinLowEdge(1)
    upEdge = hist.GetBinLowEdge(nbins+1)
    hNew = TH1F(f'{hist.GetName()}_new', '', nbins, lowEdge*multiplier, upEdge*multiplier)
    for i in range(0, nbins+2):
        hNew.SetBinContent(i, hist.GetBinContent(i))
        hNew.SetBinError(i, hist.GetBinError(i))
    return hNew
