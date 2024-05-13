#ifndef FEMPY_CORRELATIONFITTER_HXX_
#define FEMPY_CORRELATIONFITTER_HXX_

#include <map>
#include <string>
#include <tuple>
#include <stdexcept>

#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TFitResultPtr.h"
#include "FitFunctions.cxx"

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
        
        // cout << "Adding " << name << endl;
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
        cout << "First bin content: " << splineHisto->GetBinContent(1) << endl;
        cout << "First bin error: " << splineHisto->GetBinError(1) << endl;
        TSpline3* sp3 = new TSpline3(hist);
        
        this->Add(name, sp3, pars, addmode);
    }

    void Add(TString name, TGraph* graph, std::vector<std::tuple<std::string, double, double, double>> pars, std::string addmode) {
        TSpline3* sp3 = new TSpline3(graph->GetTitle(), graph);
        this->Add(name, sp3, pars, addmode);
    }

    void Add(TString name, TSpline3* spline, std::vector<std::tuple<std::string, double, double, double>> pars, std::string addmode) {
        cout << "Spline eval first bin: " << spline->Eval(2) << endl;
        cout << "Spline eval second bin: " << spline->Eval(6) << endl;
        cout << "Spline eval flat: " << spline->Eval(89) << endl;
        cout << "Spline eval flat: " << spline->Eval(570) << endl;
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

    std::vector<TF1 *> GetComponents(std::vector<int> compsidxs, std::vector<bool> onbaseline, int idxBaseline=-1) {
        if(!this->fFit) {
            throw std::invalid_argument("Fit not performed, component cannot be evaluated!");
        }
        if(compsidxs.size() != onbaseline.size()) {
            throw std::invalid_argument("Number of components does not match onbaseline information!");
        }

        std::vector<TF1 *> comps;
        for(int iToBeSavedComp=0; iToBeSavedComp<compsidxs.size(); iToBeSavedComp++) { 
            if(idxBaseline != -1 && onbaseline[iToBeSavedComp]) {
                // cout << "Evaluating with baseline" << endl;
                // cout << this->fFitFuncComps[compsidxs[iToBeSavedComp]] << endl;
                // cout << "Evaluating with baseline" << endl;
                // cout << this->fFitFuncEval[compsidxs[iToBeSavedComp]]->Eval(2) << endl;
                // cout << this->fNorms[idxBaseline] << endl;
                // cout << this->fFitFuncEval[idxBaseline]->Eval(2) << endl;
                comps.push_back(new TF1(this->fFitFuncComps[compsidxs[iToBeSavedComp]],
                    [&, this, compsidxs, iToBeSavedComp, idxBaseline](double *x, double *pars) {
                       return ( this->fFitFuncEval[compsidxs[iToBeSavedComp]]->Eval(x[0]) - 
                              ( this->fFitFuncEval[idxBaseline]->Eval(x[0])/this->fNorms[idxBaseline] ) )
                              / this->fNorms[compsidxs[iToBeSavedComp]];  
                }, this->fFitRangeMin, this->fFitRangeMax, 0));
                // return compWithoutNormAndBaseline;
            } else {
                // cout << "Evaluating without baseline" << endl;
                // cout << compsidxs[iToBeSavedComp] << " " << this->fFitFuncEval[compsidxs[iToBeSavedComp]]->Eval(2) << this->fNorms[compsidxs[iToBeSavedComp]] << endl;
                comps.push_back(new TF1(this->fFitFuncComps[compsidxs[iToBeSavedComp]],
                    [&, this, compsidxs, iToBeSavedComp](double *x, double *pars) {
                       return this->fFitFuncEval[compsidxs[iToBeSavedComp]]->Eval(x[0]) / this->fNorms[compsidxs[iToBeSavedComp]];  
                }, this->fFitRangeMin, this->fFitRangeMax, 0));
                // return compWithoutNorm;
            }
        }
        return comps;
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
                return this->fFitFuncEval[iFunc];
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
            cout << "iPar" << iPar << ": " << std::get<0>(pars) << " " << std::get<1>(pars) << " " << std::get<2>(pars) << " " << std::get<3>(pars) << endl;        
            this->fFit->SetParName(iPar, std::get<0>(pars).data());
            if (std::get<2>(pars) >= std::get<3>(pars)) {
                this->fFit->FixParameter(iPar, std::get<1>(pars));
            } else {
                this->fFit->SetParameter(iPar, std::get<1>(pars));
                this->fFit->SetParLimits(iPar, std::get<2>(pars), std::get<3>(pars));
            }
        }
    }

    // Perform the fit
    TFitResultPtr Fit() {
        // cout << "----------- Fit Parameter initialization -----------" << endl;
        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            double lowParLimit;
            double uppParLimit;
            this->fFit->GetParLimits(iPar, lowParLimit, uppParLimit);
            // cout << "iPar" << iPar << ": " << this->fFit->GetParName(iPar) << " " << this->fFit->GetParameter(iPar) 
            //  << " " << lowParLimit << " " << uppParLimit << endl;
        }

        // cout << "Bin content: " << fFitHist->GetBinContent(1) << endl;
        TFitResultPtr fitResults = fFitHist->Fit(this->fFit, "SMR+0", "");
        return fitResults;
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

    /*
    Define a canvas before calling this function and pass gPad as TVirtualPad
    */
    void Draw(TVirtualPad *pad, std::vector<TString> legLabels, std::vector<double> legCoords, int linesThickness, 
              std::vector<bool> onBaseline, std::vector<double> shifts, int basIdx=-1, std::vector<TString> addComps = {""},
              double lowRangeUser=0.0, double uppRangeUser=1.05, std::string title=";k* (MeV/c);C(k*)") {

        EvaluateComponents(basIdx, onBaseline, shifts, addComps); 
        pad->cd();
        double yMinDraw = lowRangeUser;
        double yMaxDraw = uppRangeUser + fFitHist->GetMaximum();
        
        TLegend *legend = new TLegend(legCoords[0], legCoords[1], legCoords[2], legCoords[3]);
        legend->AddEntry(this->fFitHist, legLabels[0].Data(), "lp");
        legend->AddEntry(this->fFit, legLabels[1].Data(), "l");

        //gPad->DrawFrame(fFitRangeMin, yMinDraw, fFitRangeMax, yMaxDraw, title.data());
        gPad->DrawFrame(1, yMinDraw, fFitRangeMax, yMaxDraw, title.data());
        //gPad->DrawFrame(fFitRangeMin, yMinDraw, 2000, yMaxDraw, title.data());
        //gPad->DrawFrame(fFitRangeMin, yMinDraw, 2000, yMaxDraw, Form("%s;%s;%s", 
        //                this->fFitHist->GetTitle(), this->fFitHist->GetXaxis()->GetTitle(),
        //                this->fFitHist->GetYaxis()->GetTitle()));
        
        if(basIdx == -1){
            std::cerr << "Warning: Baseline is not fixed!" << std::endl;
        }
        std::vector<TF1 *> gaussians;
        std::vector<Color_t> colors = {kMagenta + 3, kAzure + 2, kGreen, kOrange, kBlue + 2, kCyan, kBlack, kGreen+2};
        cout << "Number of functions to evaluate " << fFitFuncEval.size() << endl; 
        cout << "Number of legend entries " << legLabels.size() << endl; 
        for(int iFuncEval=0; iFuncEval<fFitFuncEval.size(); iFuncEval++) {
            // cout << "Iteration: " << iFuncEval << " " << fFitFuncEval.size() << endl;
            this->fFitFuncEval[iFuncEval]->SetNpx(300);
            this->fFitFuncEval[iFuncEval]->SetLineColor(colors[iFuncEval]); //.data());
            this->fFitFuncEval[iFuncEval]->SetLineWidth(linesThickness);
            this->fFitFuncEval[iFuncEval]->DrawF1(fFitRangeMin+1,fFitRangeMax,"same");
            this->fFitFuncEval[iFuncEval]->DrawF1(1,fFitRangeMax,"same");
            cout << fFitFuncComps[iFuncEval] << endl; 
            if(fFitFuncComps[iFuncEval] == "gaus") {
                cout << "Ciao gaus " << fFitFuncComps[iFuncEval] << endl;
                gaussians.push_back(new TF1(Form("Gaus%i", iFuncEval),
                    [&, this, iFuncEval](double *x, double *pars) {
                    return this->fFitFuncEval[iFuncEval]->Eval(x[0]);
                    }, fFitRangeMin, fFitRangeMax, 0));
                gaussians.back()->SetLineColor(colors[iFuncEval]); //.data());
                gaussians.back()->DrawF1(fFitRangeMin+1,fFitRangeMax,"same");
                gaussians.back()->Draw("same");
            }
            this->fFitFuncEval[iFuncEval]->Draw("same");
            cout << "Evaluating component at 100 MeV: " << this->fFitFuncEval[iFuncEval]->Eval(120) << endl; 
            pad->Update();
            cout << "Legend label: " << legLabels[iFuncEval+2] << endl;
            if(legLabels[iFuncEval+2].Contains("lambda_flat")) continue;
            legend->AddEntry(this->fFitFuncEval[iFuncEval], legLabels[iFuncEval+2].Data(), "l");
        }
        cout << "Evaluate last fit component: " << this->fFitFuncEval[fFitFuncEval.size()-1]->Eval(120) << " "
                                                << this->fFitFuncEval[fFitFuncEval.size()-1]->Eval(130) << " "
                                                << this->fFitFuncEval[fFitFuncEval.size()-1]->Eval(140) << " "
                                                << this->fFitFuncEval[fFitFuncEval.size()-1]->Eval(150) << " "
                                                << this->fFitFuncEval[fFitFuncEval.size()-1]->Eval(160) << " "
                                                << this->fFitFuncEval[fFitFuncEval.size()-1]->Eval(170) << " "
                                                << this->fFitFuncEval[fFitFuncEval.size()-1]->Eval(180) << " " << endl;

        // for(int iGaus=0; iGaus<gaussians.size(); iGaus++) {
        //     gaussians[iGaus]->DrawF1(fFitRangeMin+1,fFitRangeMax,"same");
        // }

        this->fFit->SetNpx(300);
        this->fFit->SetLineColor(kRed);
        this->fFit->SetLineWidth(linesThickness);
        //this->fFit->DrawF1(fFitRangeMin+1,fFitRangeMax,"same");
        this->fFit->DrawF1(1,fFitRangeMax,"same");
        pad->Update();

        // cout << "Drawing" << endl;
        fFitHist->GetYaxis()->SetRangeUser(yMinDraw, yMaxDraw); 
        fFitHist->SetMarkerSize(0.1);
        fFitHist->SetMarkerStyle(20);
        fFitHist->SetMarkerColor(kBlack);
        fFitHist->SetLineColor(kBlack);
        fFitHist->SetLineWidth(3);
        fFitHist->Draw("same pe");
        pad->Update();

        legend->SetBorderSize(0);
        legend->SetTextSize(0.045);
        legend->Draw("same");
        pad->Update();
        cout << "Drawn!" << endl;

    }

    void Debug() {
        cout << endl;
        cout << endl;
        cout << endl;
        cout << "########## Debugging fit components ##########" << endl;
        cout << "Total spline terms: " << this->fFitSplines.size() << endl;
        cout << "Total func terms: " << this->fFitFunc.size() << endl;
        cout << "--------------------------" << endl;    
        int nTerms = this->fFitFunc.size() + this->fFitSplines.size();
        for(int iTerm=0; iTerm<nTerms; iTerm++) {
            cout << "Term name: " << this->fFitFuncComps[iTerm] << endl;
            cout << "Add Mode: " << this->fAddModes[iTerm] << endl;
            cout << "Number of parameters of term: " << this->fNPars[iTerm+1] << endl;
            cout << "Normalization factor: " << this->fNorms[iTerm] << endl;
            int startPar = accumulate(fNPars.begin(), std::next(fNPars.begin(), iTerm+1), 0);
            for(int iPar=1; iPar<this->fNPars[iTerm+1]; iPar++) {
                cout << "Parameter " << iPar-1 << " value: " << this->fFit->GetParameter(startPar + iPar) << endl;
            }
            cout << "--------------------------" << endl;
        }
        cout << endl;
        cout << endl;
        cout << endl;
    }

   private:

    void EvaluateComponents(int basIdx, std::vector<bool> onBaseline, std::vector<double> shifts, std::vector<TString> addComps = {""}) {
        
        double globNorm;
        if(this->fGlobNorm) {
            globNorm = this->fFit->GetParameter(this->fFit->GetNpar()-1);
        } else {
            globNorm = 1;
        }
        cout << "GLOBAL NORM: " << globNorm << endl;
        int normParNumber = 0;
        for(int iNorm=0; iNorm<fNPars.size()-1; iNorm++) {
            fNorms.push_back(this->fFit->GetParameter(normParNumber));
            // cout << "Norm: " << this->fFit->GetParameter(normParNumber) << endl;
            normParNumber += fNPars[iNorm+1];
        }
        for(int iShift=0; iShift<shifts.size()-1; iShift++) {
            cout << "Shift: " << shifts[iShift] << endl;
            normParNumber += fNPars[iShift+1];
        }
        std::vector<TF1*> components;
        // cout << "bas " << basIdx << endl;
        int nTerms = this->fFitFunc.size() + this->fFitSplines.size() - this->fGlobNorm;
        int setPars = 0;
        int iSpline = 0;
        int iFunc = 0;
        cout << "Fit Function terms: " << nTerms << endl;
        cout << "Fit Function number of parameters: " << this->fFit->GetNpar() << endl;
        for(int iTerm=0; iTerm<nTerms; iTerm++) {
            std::string compName = static_cast<std::string>(this->fFitFuncComps[iTerm]);
            cout << "Name of the component: " << compName << endl;
            if(this->fFitFuncComps[iTerm].Contains("spline") && !this->fFitFuncComps[iTerm].Contains("spline3")) {
                cout << "SPLINE" << endl;
                cout << "Set par " << this->fFit->GetParameter(setPars) << endl;
                cout << "Set par +1 " << this->fFit->GetParameter(setPars+1) << endl;
                cout << endl;
                components.push_back(new TF1(Form("iComp_%.0f", iTerm),
                        [&, this, iSpline, iTerm, setPars](double *x, double *pars) {
                            if(this->fNPars[iTerm+1] == 2) {
                                return fFitSplines[iSpline]->Eval(x[0] - this->fFit->GetParameter(setPars+1));
                            } else {
                                return fFitSplines[iSpline]->Eval(x[0]);
                            }
                        }, fFitRangeMin, fFitRangeMax, 0));
                iSpline++;
            } else {
                cout << "Evaluating function " << iFunc << " nPars: " << fNPars[iTerm+1] << endl;
                components.push_back(new TF1(Form("iComp_%.0f", iTerm), fFitFunc[iFunc], fFitRangeMin, fFitRangeMax, fNPars[iTerm+1]));
                iFunc++;                
                for(int iPar=1; iPar<fNPars[iTerm+1]; iPar++) {
                    // cout << "Setting iPar " << iPar << ": " << this->fFit->GetParameter(setPars + iPar) << endl;
                    components.back()->FixParameter(iPar-1, this->fFit->GetParameter(setPars + iPar));
                }
                // cout << endl;
            }
            setPars += fNPars[iTerm+1];
        }

        for(int iFunc=0; iFunc<nTerms; iFunc++) {
            bool toBeSummed = false;
            for(int iAddComp=0; iAddComp<addComps.size(); iAddComp++) {
                
                if(addComps[iAddComp].Contains(std::to_string(iFunc))) {
                    toBeSummed = true;
                    continue;
                }
            }
            
            if(!toBeSummed) {
                // cout << "Component: " << this->fFitFuncComps[iFunc] << endl;
                // cout << iFunc << endl;
                // cout << onBaseline[iFunc] << endl;
                // cout << basIdx << endl;
                // cout << this->fFitFuncComps[iFunc] << endl;
                this->fFitFuncEval.push_back(new TF1(this->fFitFuncComps[iFunc],
                    [&, this, iFunc, components, onBaseline, basIdx, globNorm, shifts](double *x, double *pars) {
                        if(onBaseline[iFunc]) {
                            return globNorm * (this->fNorms[iFunc]*components[iFunc]->Eval(x[0]) + 
                                   this->fNorms[basIdx]*components[basIdx]->Eval(x[0])) + shifts[iFunc];
                        } else if(this->fFitFuncComps[iFunc].Contains("Lednicky")) {
                            return components[iFunc]->Eval(x[0]);
                        } else {
                            return globNorm * this->fNorms[iFunc]*components[iFunc]->Eval(x[0]) + shifts[iFunc];
                        }},
                    fFitRangeMin, fFitRangeMax, 0));            
            }
        }        
            
        if(addComps[0] != "") {
            for(int iAddComp=0; iAddComp<addComps.size(); iAddComp++) {
                int compNumber = this->fFitFuncEval.size();
                this->fFitFuncEval.push_back(new TF1(Form("Comp_%i", compNumber),
                        [&, this, nTerms, addComps, iAddComp, components, onBaseline, basIdx, globNorm](double *x, double *pars) {
                        double sum = 0.;
                        bool addBaseline = true;
                        for(int iComp=0; iComp<nTerms; iComp++) {
                            if(addComps[iAddComp].Contains(std::to_string(iComp))) {
                                sum += this->fNorms[iComp]*components[iComp]->Eval(x[0]);
                                if(!onBaseline[iComp]) addBaseline=false;
                            }
                        }
                        if(addBaseline) {
                            sum += this->fNorms[basIdx]*components[basIdx]->Eval(x[0]);
                        }
                        return globNorm * sum;}, fFitRangeMin, fFitRangeMax, 0));
            }
        }
        // cout << "Evaluated" << endl;
    }

    TH1 *fFitHist = nullptr;
    TF1 *fFit = nullptr;

    std::vector<double (*)(double *x, double *par)> fFitFunc;       // List of function describing each term of the CF model
    std::vector<TString> fFitFuncComps;                             // Function names of fit components
    std::vector<TSpline3 *> fFitSplines;                            // Spline fit components
    std::vector<TF1*> fFitFuncEval;                                 // Fit components evaluated after the fitting
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
