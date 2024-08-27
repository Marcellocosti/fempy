#ifndef FEMPY_CORRELATIONFITTER_HXX_
#define FEMPY_CORRELATIONFITTER_HXX_

#include <map>
#include <string>
#include <tuple>
#include <stdexcept>

#include "TF1.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TFitResultPtr.h"
#include "FitFunctions.cxx"

#if LOG_LEVEL_FIT
#define DEBUG(msg) std::cout << msg << std::endl
#else
#define DEBUG(msg)
#endif

class CorrelationFitter {
   public:
    CorrelationFitter(TH1 *fithist, double fitRangeMin, double fitRangeMax, 
                         double rejectMin=0, double rejectMax=-1) {
        this->fFitHist = reinterpret_cast<TH1 *>(fithist->Clone());
        this->fFitRangeMin = fitRangeMin;
        this->fFitRangeMax = fitRangeMax;
        this->fRejectMin = rejectMin;
        this->fRejectMax = rejectMax;
        this->fGlobNorm = false;
        this->fNPars = {0};  // The first parameter has index zero
    }

    void DrawSpline(TVirtualPad *pad, TH1* hist, 
                    std::string name=";k* (MeV/c);Counts") {
        pad->cd();
        
        TH1D* histo = static_cast<TH1D*>(hist);
        double yMaxDraw = histo->GetBinContent(histo->GetMaximumBin())*1.2;
        
        gPad->DrawFrame(fFitRangeMin, 0, fFitRangeMax, yMaxDraw, name.data());

        TSpline3* sp3 = new TSpline3(histo);
        sp3->SetNpx(300);
        sp3->SetLineColor(kRed);
        sp3->Draw("same");

        histo->SetMarkerSize(0.3);
        histo->SetMarkerStyle(20);
        histo->SetMarkerColor(kBlack);
        histo->SetLineColor(kBlack);
        histo->SetLineWidth(3);
        histo->Draw("same pe");
        pad->Update();
    }

    void DrawSpline(TVirtualPad *pad, TGraph* graph, 
                    std::string name=";k* (MeV/c);Counts") {
        
        pad->cd();
        double yMaxDraw = TMath::MaxElement(graph->GetN(), graph->GetY())*1.2; 
        double yMinDraw = TMath::MinElement(graph->GetN(), graph->GetY())*0.8; 
        gPad->DrawFrame(fFitRangeMin, yMinDraw, fFitRangeMax, yMaxDraw, name.data());

        TSpline3* sp3graph = new TSpline3(graph->GetTitle(), graph);
        sp3graph->SetNpx(300);
        sp3graph->SetLineColor(kRed);
        sp3graph->Draw("same");
        
        graph->SetMarkerSize(0.3);
        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(kBlack);
        graph->SetLineColor(kBlack);
        graph->SetLineWidth(3);
        graph->Draw("same pe");
        pad->Update();
    }

    /*
    Add a term to the CF model. Available options:
        - pol0
        - bw (Breit-Wigner)

    In pyroot, the fit parameters (pars) are specified as a list of tuples.

    The syntax is: `(name, initialValue, lowerLim, upperLim)`. To fix a parameter set `lowerLim > upperLim`, e.g
    `fitter.Add('bw',
    [('p0', 0.9, 1, -1)])`

    Example:
    ```fitter.Add('bw', [
        ('yield', 0.03, 0.001, 0.1),
        ('mean', 0.126, 0.125, 0.127),
        ('gamma', 0.004, 0.003, 0.005),
    ])```

    It is your responsibility to give the correct number of parameters. All fit parameters must be initialized.
    */

    void Add(TString name, std::vector<std::tuple<std::string, double, double, double>> pars, std::string addmode) {
        
        DEBUG("Adding " << name);
        if(functions.find(name)!=functions.end()){
            this->fFitFunc.push_back(std::get<0>(functions[name]));
            this->fFitFuncComps.push_back(name);
            // -1 needed because pars includes the norm of the term
            if(pars.size()-1 != std::get<1>(functions[name])) {
                printf("Error: wrong number of parameters for function '%s'!\n", name.Data());
                exit(1);                
            } else {
                this->fNPars.push_back(pars.size());
            }
        } else {
            printf("Error: function '%s' is not implemented. Exit!", name.Data());
            exit(1);
        }

        this->fAddModes.push_back(addmode);
        
        // Save fit settings
        for (const auto &par : pars) {
            this->fFitPars.insert({this->fFitPars.size(), par});
        }
    }

    void Add(TString name, TH1* hist, std::vector<std::tuple<std::string, double, double, double>> pars, std::string addmode) {
        TH1D *splineHisto = static_cast<TH1D*>(hist);
        TSpline3* sp3 = new TSpline3(hist);
        this->Add(name, sp3, pars, addmode);
    }

    void Add(TString name, TGraph* graph, std::vector<std::tuple<std::string, double, double, double>> pars, std::string addmode) {
        TSpline3* sp3 = new TSpline3(graph->GetTitle(), graph);
        this->Add(name, sp3, pars, addmode);
    }

    void Add(TString name, TSpline3* spline, std::vector<std::tuple<std::string, double, double, double>> pars, std::string addmode) {
        this->fFitSplines.push_back(spline);
        this->fFitFuncComps.push_back(name);
        this->fAddModes.push_back(addmode); 
        this->fNPars.push_back(pars.size());
        // Save fit settings
        for (const auto &par : pars) {
            this->fFitPars.insert({this->fFitPars.size(), par});
        }
    }

    void AddGlobNorm(std::string globnorm, double initval, double lowedge, double uppedge) {
        this->fGlobNorm = true;
        this->fFitFunc.push_back(std::get<0>(functions[globnorm]));
        this->fFitFuncComps.push_back(globnorm);
        // -1 needed because pars includes the norm of the term
        
        this->fNPars.push_back(1);
        this->fAddModes.push_back("*");
        
        // Save fit settings
        std::tuple<std::string, double, double, double> init = {globnorm, initval, lowedge, uppedge}; 
        this->fFitPars.insert({this->fFitPars.size(), init});
    }

    TH1D *GetComponentPars(int icomp) {
        if(!this->fFit) {
            throw std::invalid_argument("Fit not performed, component cannot be evaluated!");
        }
        int startPar = accumulate(fNPars.begin(), std::next(fNPars.begin(), icomp+1), 0);
        std::vector<double> compPars;
        TH1D *histoPars = new TH1D("histoPars_" + this->fFitFuncComps[icomp], "histoPars_" + this->fFitFuncComps[icomp], 
                                    this->fNPars[icomp+1], 0, this->fNPars[icomp+1]);
        
        for(int iCompPar=0; iCompPar<this->fNPars[icomp+1]; iCompPar++) {
            histoPars->SetBinContent(iCompPar+1, this->fFit->GetParameter(startPar+iCompPar));
        }
        return histoPars;
    }

    TF1 *GetGenuine() {

        int genuineIdx;
        for(int iFunc=0; iFunc<this->fFitFuncComps.size(); iFunc++) {
            if(this->fFitFuncComps[iFunc].Contains("Lednicky")) {
                genuineIdx = iFunc;
            }
        }
        
        int startGenuinePar;
        for(int iPar=0; iPar<this->fFit->GetNpar(); iPar++) {
            std::string parName = this->fFit->GetParName(iPar);
            if (parName.find("re_a0") != std::string::npos) {
                startGenuinePar = iPar;
            }
        }

        TF1 *fGenuine = new TF1("fGenuine", std::get<0>(functions[this->fFitFuncComps[genuineIdx]]), fFitRangeMin, fFitRangeMax, 
                                std::get<1>(functions[this->fFitFuncComps[genuineIdx]]));
        for(int iGenPar=0; iGenPar<fGenuine->GetNpar(); iGenPar++) {
            fGenuine->FixParameter(iGenPar, this->fFit->GetParameter(iGenPar + startGenuinePar));
        }

        return fGenuine;
    }

    TH1D *SaveFreeFixPars() {
        if(!this->fFit) {
            throw std::invalid_argument("Fit not performed, component cannot be evaluated!");
        }
        
        TH1D *histoFreeFixPars = new TH1D("hFreeFixPars", "hFreeFixPars", this->fFit->GetNpar()+2, 0, this->fFit->GetNpar()+2);
        
        histoFreeFixPars->SetBinContent(1, this->fFitRangeMin);
        histoFreeFixPars->GetXaxis()->SetBinLabel(1, "Low fit edge");
        
        histoFreeFixPars->SetBinContent(2, this->fFitRangeMax);
        histoFreeFixPars->GetXaxis()->SetBinLabel(2, "Upp fit edge");

        double lowParEdge = 0.;
        double uppParEdge = 0.;
        for(int iPar=0; iPar<this->fFit->GetNpar(); iPar++) {
            this->fFit->GetParLimits(iPar, lowParEdge, uppParEdge);
            if(lowParEdge >= uppParEdge) {
                histoFreeFixPars->SetBinContent(iPar+3, -1);
            } else {
                histoFreeFixPars->SetBinContent(iPar+3, 1);
            }
            histoFreeFixPars->GetXaxis()->SetBinLabel(iPar+3, this->fFit->GetParName(iPar));
        }

        histoFreeFixPars->SetStats(0);
        histoFreeFixPars->GetXaxis()->SetLabelSize(100);
        return histoFreeFixPars;
    }

    TH1D *SaveFitPars() {
        if(!this->fFit) {
            throw std::invalid_argument("Fit not performed, component cannot be evaluated!");
        }
        
        TH1D *histoPars = new TH1D("hFitPars", "hFitPars", this->fFit->GetNpar(), 0, this->fFit->GetNpar());
        for(int iPar=0; iPar<this->fFit->GetNpar(); iPar++) {
            histoPars->SetBinContent(iPar+1, this->fFit->GetParameter(iPar));
            histoPars->GetXaxis()->SetBinLabel(iPar+1, this->fFit->GetParName(iPar));
        }
        histoPars->SetStats(0);
        histoPars->GetXaxis()->SetLabelSize(100);
        return histoPars;
    }

    TH1D *SaveFitParsSplitComponents(std::vector<std::string> subcompsmothers, std::vector<int> compstosplit, std::vector<std::vector<int>> compsnormsidx, 
                                     std::vector<std::vector<std::string>> compsnames, std::vector<std::vector<std::string>> normslabels) {
        if(!this->fFit) {
            throw std::invalid_argument("Fit not performed, component cannot be evaluated!");
        }
        
        // Calculate the number of bins needed, considering that each component will have its normalization factor
        int nCompsPars = 0;
        for(int iSplitComp=0; iSplitComp<compsnames.size(); iSplitComp++) {
            DEBUG("Extracting subcomps of funct " << subcompsmothers[iSplitComp]);
            for(int iSubComp=0; iSubComp<compsnames[iSplitComp].size(); iSubComp++) {
                DEBUG(compsnames[iSplitComp][iSubComp]);
                nCompsPars += std::get<1>(functions[compsnames[iSplitComp][iSubComp]]) + 1;
            }
            DEBUG("");
        }
        int nBinNewCompPar = this->fFit->GetNpar();
        DEBUG("Global fit function parameters no. " << nBinNewCompPar << endl; 
              cout << "Total parameters of subcomponents: " << nCompsPars);
        TH1D *histoAllCompsPars = new TH1D("hAllCompsPars", "hAllCompsPars", this->fFit->GetNpar() + nCompsPars, 0, this->fFit->GetNpar() + nCompsPars);
        DEBUG("Number of bins of the histo containing subcomps pars: " << histoAllCompsPars->GetNbinsX());
        if(this->fGlobNorm) {
            DEBUG("Global norm"); 
            // put the global norm at the beginning
            histoAllCompsPars->SetBinContent(1, this->fFit->GetParameter(this->fFit->GetNpar()-1));
            histoAllCompsPars->GetXaxis()->SetBinLabel(1, this->fFit->GetParName(this->fFit->GetNpar()-1));
            for(int iPar=0; iPar<this->fFit->GetNpar()-1; iPar++) {
                histoAllCompsPars->SetBinContent(iPar+2, this->fFit->GetParameter(iPar));
                histoAllCompsPars->GetXaxis()->SetBinLabel(iPar+2, this->fFit->GetParName(iPar));
            }
        } else {
            DEBUG("No global norm"); 
            for(int iPar=0; iPar<this->fFit->GetNpar(); iPar++) {
                DEBUG("Setting ipar " << iPar << " of " << this->fFit->GetNpar()); 
                histoAllCompsPars->SetBinContent(iPar+1, this->fFit->GetParameter(iPar));
                histoAllCompsPars->GetXaxis()->SetBinLabel(iPar+1, this->fFit->GetParName(iPar));
            }
        }

        DEBUG("Parameters set");
        DEBUG("Size of compstosplit: " << compstosplit.size());
        for(int iSplitComp=0; iSplitComp<compstosplit.size(); iSplitComp++) {
            int previousCompsPars = 0; 
            DEBUG("Size of fNpars: " << fNPars.size());
            DEBUG("Compstosplit: " << compstosplit[iSplitComp]);
            DEBUG("fGlobNorm: " << this->fGlobNorm);
            // int startCompPar = accumulate(fNPars.begin(), std::next(fNPars.begin(), 1), 0) + this->fGlobNorm + compstosplit[iSplitComp];
            int startCompPar = accumulate(fNPars.begin(), std::next(fNPars.begin(), compstosplit[iSplitComp]), 0) + this->fGlobNorm + compstosplit[iSplitComp];
            DEBUG("Initializing subcomponents");
            for(int iComp=0; iComp<compsnames[iSplitComp].size(); iComp++) {

                    std::cout << std::showpos;
                    cout.precision(4);
                    std::cout << std::scientific;
                    
                    DEBUG("--------------------------------"); 
                    DEBUG("SUBCOMPONENT " << compsnames[iSplitComp][iComp] << endl;
                          cout << "Start picking the parameters from  bin no. " << startCompPar);

                    // Setting the norm
                    if(compsnormsidx[iSplitComp][iComp] > 0) {
                        DEBUG(std::setw(30) << "[Norm] Filling bin no. " << nBinNewCompPar+1;
                        cout << " with " << histoAllCompsPars->GetBinContent(startCompPar + compsnormsidx[iSplitComp][iComp] - 1);
                        cout << ", copied from bin no. " << startCompPar + compsnormsidx[iSplitComp][iComp] - 1);
                        
                        histoAllCompsPars->SetBinContent(nBinNewCompPar+1, histoAllCompsPars->GetBinContent(startCompPar+compsnormsidx[iSplitComp][iComp]-1));
                    } 
                    if (compsnormsidx[iSplitComp][iComp] < 0) {
                        DEBUG(std::setw(30) << "[Norm] Filling bin no. " << nBinNewCompPar+1;
                        cout << " with " << 1-histoAllCompsPars->GetBinContent(startCompPar + compsnormsidx[iSplitComp][iComp] - 1);
                        cout << ", complementary of bin no. " << startCompPar + compsnormsidx[iSplitComp][iComp] - 1);
                        histoAllCompsPars->SetBinContent(nBinNewCompPar+1, 1 - histoAllCompsPars->GetBinContent(startCompPar+compsnormsidx[iSplitComp][iComp]-1));
                    } 
                    if (compsnormsidx[iSplitComp][iComp] == 0) {
                        DEBUG(std::setw(30) << "[Norm] Filling bin no. " << nBinNewCompPar+1;
                        cout << " with unitary norm");
                        histoAllCompsPars->SetBinContent(nBinNewCompPar+1, 1);
                    } 
                    
                    histoAllCompsPars->GetXaxis()->SetBinLabel(nBinNewCompPar+1, normslabels[iSplitComp][iComp].c_str());
                    
                    int compPar = std::get<1>(functions[compsnames[iSplitComp][iComp]]);
                    for(int iCompPar=0; iCompPar<compPar; iCompPar++) {
                        DEBUG(std::setw(30) << "Filling bin no. " << nBinNewCompPar+iCompPar+2;
                        cout << " with " << histoAllCompsPars->GetBinContent(startCompPar+iCompPar+previousCompsPars+1);
                        cout << ", copied from bin no. " << startCompPar+iCompPar+previousCompsPars+1);
                        histoAllCompsPars->SetBinContent(nBinNewCompPar+iCompPar+2, histoAllCompsPars->GetBinContent(startCompPar+iCompPar+previousCompsPars+1));
                        histoAllCompsPars->GetXaxis()->SetBinLabel(nBinNewCompPar+iCompPar+2, this->fFit->GetParName(startCompPar+iCompPar+previousCompsPars-this->fGlobNorm));
                    }       
                    std::cout << std::noshowpos;  
                    previousCompsPars += compPar;
                    nBinNewCompPar += compPar + 1;
                    DEBUG("Finished component");
                    DEBUG("--------------------------------");    
            }
        }

        histoAllCompsPars->SetStats(0);
        histoAllCompsPars->GetXaxis()->SetLabelSize(100);
        return histoAllCompsPars;
    }

    TGraph *SaveScatPars() {
        if(!this->fFit) {
            throw std::invalid_argument("Fit not performed, component cannot be evaluated!");
        }
        double reScatLength;
        double imScatLength;
        double reScatLengthError;
        double imScatLengthError;
        double effRange;
        for(int iPar=0; iPar<this->fFit->GetNpar(); iPar++) {
            std::string parName = this->fFit->GetParName(iPar);
            if (parName.find("re_a0") != std::string::npos) {
                reScatLength = this->fFit->GetParameter(iPar);
                reScatLengthError = this->fFit->GetParError(iPar);
            } 
            if (parName.find("im_a0") != std::string::npos) {
                imScatLength = this->fFit->GetParameter(iPar);
                imScatLengthError = this->fFit->GetParError(iPar);
            } 
        }
        DEBUG("--------------------------------");

        TGraphErrors *gScatPars = new TGraphErrors(1, &reScatLength, &imScatLength, &reScatLengthError, &imScatLengthError); // 1);
        gScatPars->SetTitle("Scattering parameters;Re_a0;Im_a0");
        gScatPars->SetName("gScatPars");
        // gScatPars->SetPoint(1, reScatLength, imScatLength);
        // gScatPars->SetPointError(1, reScatLengthError, imScatLengthError);
        return gScatPars;
    }
    
    TH1D *PullDistribution() {
        if(!this->fFit) {
            throw std::invalid_argument("Fit not performed, pulls cannot be calculated!");
        }
    
        std::vector<double> pulls;
        for(int iBin=0; iBin<this->fFitHist->GetNbinsX(); iBin++) {    
            if(this->fFitHist->GetBinCenter(iBin+1) >= this->fFitRangeMin &&
               this->fFitHist->GetBinCenter(iBin+1) <= this->fFitRangeMax) {
                    pulls.push_back( (this->fFitHist->GetBinContent(iBin+1) - this->fFit->Eval(this->fFitHist->GetBinCenter(iBin+1))) /         
                                      this->fFitHist->GetBinError(iBin+1));
            }
        }

        TH1D *histoPulls = new TH1D("hPulls", "hPulls", pulls.size(), this->fFitRangeMin, this->fFitRangeMax);
        for(int iBin=0; iBin<this->fFitHist->GetNbinsX(); iBin++) {    
            histoPulls->SetBinContent(iBin+1, pulls[iBin]);         
        }

        return histoPulls;
    }

    void BuildFitFunction() {
        cout << "----------- Building the fit function -----------" << endl;
        // Build the fit function
        this->fFit = new TF1(
            "fFit",
            [&, this](double *x, double *pars) {
                if(x[0] >= this->fRejectMin && x[0] <= this->fRejectMax) {
                    TF1::RejectPoint();
                } 
                double result = 0;
                int nTerms = this->fFitFunc.size() + this->fFitSplines.size();
                int nPar = 0;
                int nSplineComp = 0;
                int nFuncComp = 0;
                for (int iTerm = 0; iTerm < nTerms; iTerm++) {
                    if(fFitFuncComps[iTerm].Contains("splinehisto")) {
                        if(fAddModes[iTerm] == "*") {
                            double partResult = result; 
                            if(this->fNPars[iTerm+1] == 2){
                                result = pars[nPar]*this->fFitSplines[nSplineComp]->Eval(x[0] - pars[nPar+1])*partResult;
                            } else {
                                result = pars[nPar]*this->fFitSplines[nSplineComp]->Eval(x[0])*partResult;
                            }
                        }
                        else {
                            if(this->fNPars[iTerm+1] == 2){
                                result += pars[nPar]*this->fFitSplines[nSplineComp]->Eval(x[0] - pars[nPar+1]);
                            } else {
                                result += pars[nPar]*this->fFitSplines[nSplineComp]->Eval(x[0]);
                            }
                        }
                        nSplineComp++;
                    } else {
                        auto func = this->fFitFunc[nFuncComp];
                        if(fAddModes[iTerm] == "*") {
                            double partResult = result; 
                            result = pars[nPar]*func(x, &pars[nPar+1])*partResult;
                        } else {
                            result += pars[nPar]*func(x, &pars[nPar+1]);
                        }
                        nFuncComp++;
                    }                            
                    nPar += this->fNPars[iTerm+1];   
                }
            return result;}, 
            fFitRangeMin, fFitRangeMax, this->fFitPars.size());

        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            auto pars = this->fFitPars[iPar];
            std::cout << std::showpos;
            cout.precision(4);
            std::cout << std::scientific;
            DEBUG(std::setw(25) << std::get<0>(pars) << "  ----->  ";
                  cout << std::setw(8) << " Init: " << static_cast<double>(std::get<1>(pars)) << ",";
                  cout << std::setw(8) << " Lowlim: " << static_cast<double>(std::get<2>(pars)) << ",";
                  cout << std::setw(8) << " Upplim: " << static_cast<double>(std::get<3>(pars)));
            std::cout << std::noshowpos;    

            this->fFit->SetParName(iPar, std::get<0>(pars).data());
            if (std::get<2>(pars) >= std::get<3>(pars)) {
                this->fFit->FixParameter(iPar, std::get<1>(pars));
            } else {
                this->fFit->SetParameter(iPar, std::get<1>(pars));
                this->fFit->SetParLimits(iPar, std::get<2>(pars), std::get<3>(pars));
            }
        }
    }

    TFitResultPtr Fit() {
        DEBUG("----------- Fit Parameter initialization -----------");
        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            double lowParLimit;
            double uppParLimit;
            this->fFit->GetParLimits(iPar, lowParLimit, uppParLimit);
            std::cout << std::showpos;
            cout.precision(4);
            std::cout << std::scientific;
            DEBUG(std::setw(25) << this->fFit->GetParName(iPar) << "  ----->  ";
                  cout << std::setw(8) << " Init: " << static_cast<double>(this->fFit->GetParameter(iPar)) << ",";
                  cout << std::setw(8) << " Lowlim: " << static_cast<double>(lowParLimit) << ",";
                  cout << std::setw(8) << " Upplim: " << static_cast<double>(uppParLimit));
            std::cout << std::noshowpos;        }

        DEBUG("Bin content: " << fFitHist->GetBinContent(1));
        TFitResultPtr fitResults = fFitHist->Fit(this->fFit, "SMR+0", "");

        #ifdef DEBUG
            Debug();
        #endif

        return fitResults;
    }

    std::vector <TH1D *> FitBootstrap(int ntries) {
        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            double lowParLimit;
            double uppParLimit;
            this->fFit->GetParLimits(iPar, lowParLimit, uppParLimit);
            std::cout << std::showpos;
            cout.precision(4);
            std::cout << std::scientific;
            DEBUG(std::setw(25) << this->fFit->GetParName(iPar) << "  ----->  ";
                  cout << std::setw(8) << " Init: " << static_cast<double>(this->fFit->GetParameter(iPar)) << ",";
                  cout << std::setw(8) << " Lowlim: " << static_cast<double>(lowParLimit) << ",";
                  cout << std::setw(8) << " Upplim: " << static_cast<double>(uppParLimit));
            std::cout << std::noshowpos;        
        }

        TFitResultPtr originalFitResults = fFitHist->Fit(this->fFit, "SMR+0", "");

        // cout << "CIAO1" << endl;
        // Initialize histogram for saving fit parameters
        std::vector <TH1D *> hFitParsBT;
        TH2D *hBootstrapSummary = new TH2D("hBootstrapSummary", "", ntries, 0, ntries, 
                                           this->fFit->GetNpar()+1, 0, this->fFit->GetNpar()+1);
        for(int iPar=0; iPar<this->fFit->GetNpar(); iPar++) {
            double parLowLim, parUppLim;
            this->fFit->GetParLimits(iPar, parLowLim, parUppLim);
            double histoMean = originalFitResults->Parameter(iPar);
            double histoBound = originalFitResults->ParError(iPar);
            hFitParsBT.push_back(new TH1D(this->fFit->GetParName(iPar), this->fFit->GetParName(iPar), 1000, 
                                          histoMean - 5*histoBound, histoMean + 5*histoBound));
        }
        hFitParsBT.push_back(new TH1D("Chi2/DOF", "Chi2/DOF", 1000, 0, 100));

        int nBinsFitCF = this->fFitRangeMax/this->fFitHist->GetBinWidth(1);

        // loop over the tries
        // cout << "CIAO2" << endl;
        for(int iTry=0; iTry<ntries; iTry++) {
            // cout << "CIAO3" << endl;
            // compute CF with gaussian sampling
            TH1D *sampledCF = new TH1D(Form("hTry_%i", iTry), "", nBinsFitCF, 
                                       this->fFitRangeMin, this->fFitRangeMax);
            // cout << "CIAO4" << endl;
            for(int iSampledBin=0; iSampledBin<sampledCF->GetNbinsX(); iSampledBin++) {
                double CFvalue = this->fFitHist->GetBinContent(iSampledBin+1);
                double CFerror = this->fFitHist->GetBinError(iSampledBin+1);
                sampledCF->SetBinContent(iSampledBin+1, gRandom->Gaus(CFvalue, CFerror));
                sampledCF->SetBinError(iSampledBin+1, CFerror);
            }
            // cout << "CIAO5" << endl;
            TFitResultPtr fitResults = sampledCF->Fit(this->fFit, "SMR+0", "");
            // cout << "CIAO6" << endl;
            for(int iPar=0; iPar<this->fFit->GetNpar(); iPar++) {
                // cout << "CIAO7" << endl;
                hFitParsBT[iPar]->Fill(fitResults->Parameter(iPar));
            }
            // cout << "CIAO8" << endl;
            hFitParsBT[this->fFit->GetNpar()]->Fill(fFit->GetChisquare() / fFit->GetNDF());
            // cout << "CIAO9" << endl;
        }

        int index = 0;
        hFitParsBT.erase(std::remove_if(hFitParsBT.begin(), hFitParsBT.end(), 
                                        [&index, originalFitResults](const auto&) {
                                            bool isParFixed = originalFitResults->IsParameterFixed(index);
                                            index++; 
                                            return isParFixed;
                                        }), hFitParsBT.end());
        return hFitParsBT;
    }

    TF1 *GetFitFunction() {
        return this->fFit;
    } 

    double GetChi2Ndf() { return fFit->GetChisquare() / fFit->GetNDF(); }

    double GetChi2NdfManual() { 
        double chi2 = 0.;
        for(int iBin=0; iBin<this->fFitHist->GetNbinsX(); iBin++) {
            if(this->fFitHist->GetBinError(iBin+1) != 0) {
                chi2 += ( (this->fFit->Eval(this->fFitHist->GetBinCenter(iBin+1)) - this->fFitHist->GetBinContent(iBin+1)) * 
                          (this->fFit->Eval(this->fFitHist->GetBinCenter(iBin+1)) - this->fFitHist->GetBinContent(iBin+1)) ) /
                          (this->fFitHist->GetBinError(iBin+1) * this->fFitHist->GetBinError(iBin+1));
            }
        }

        return chi2 / fFit->GetNDF(); 
    }

   private:

    void Debug() {
        DEBUG("");
        DEBUG("");
        DEBUG("");
        int normParNumber = 0;
        for(int iNorm=0; iNorm<fNPars.size()-1; iNorm++) {
            fNorms.push_back(this->fFit->GetParameter(normParNumber));
            normParNumber += fNPars[iNorm+1];
        }
        DEBUG("########## Debugging fit components ##########");
        DEBUG("Total spline terms: " << this->fFitSplines.size());
        DEBUG("Total func terms: " << this->fFitFunc.size());
        DEBUG("--------------------------");
        int nTerms = this->fFitFunc.size() + this->fFitSplines.size();
        DEBUG("Size of norms vector: " << this->fNorms.size());
        std::cout << std::showpos;
        cout.precision(4);
        std::cout << std::scientific;
        for(int iTerm=0; iTerm<nTerms; iTerm++) {
            DEBUG("Term name: " << this->fFitFuncComps[iTerm]);
            DEBUG("Add Mode: " << this->fAddModes[iTerm]);
            DEBUG("No. pars of term: " << std::to_string(this->fNPars[iTerm+1]));
            int startPar = accumulate(fNPars.begin(), std::next(fNPars.begin(), iTerm+1), 0);
            for(int iPar=0; iPar<this->fNPars[iTerm+1]; iPar++) {
                DEBUG(std::setw(30) << this->fFit->GetParName(startPar + iPar) << ": " << std::setw(5) << this->fFit->GetParameter(startPar + iPar));
            }
            DEBUG("--------------------------");
        }
        std::cout << std::noshowpos;
        DEBUG("");
        DEBUG("");
        DEBUG("");
    }

    TH1 *fFitHist = nullptr;
    TF1 *fFit = nullptr;

    std::vector<double (*)(double *x, double *par)> fFitFunc;       // List of function describing each term of the CF model
    std::vector<TString> fFitFuncComps;                             // Function names of fit components
    std::vector<TSpline3 *> fFitSplines;                            // Fit components evaluated after the fitting
    std::vector<int> fNPars;                                        // Keeps track of how many parameters each function has
    std::vector<std::string> fAddModes;                             // Select mode of adding the contributions to the model
    std::vector<double> fNorms;                                     // Vector saving the norm factor of each term
    bool fGlobNorm; 

    double fFitRangeMin;
    double fFitRangeMax;
    double fRejectMin;
    double fRejectMax;
    std::map<int, std::tuple<std::string, double, double, double>> fFitPars;    // List of fit parameters
};

#endif  // FEMPY_CORRELATIONFITTER_HXX_
