#ifndef FEMPY_CORRELATIONFITTER_HXX_
#define FEMPY_CORRELATIONFITTER_HXX_

#include <map>
#include <string>
#include <tuple>

#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TFitResultPtr.h"
#include "FitFunctions.cxx"

class CorrelationFitter {
   public:
    CorrelationFitter(TH1 *datahist, TH1 *mchist, double fitRangeMin, double fitRangeMax) {
        this->fDataHist = reinterpret_cast<TH1 *>(datahist->Clone());
        this->fMCHist = reinterpret_cast<TH1 *>(mchist->Clone());
        this->fFitRangeMin = fitRangeMin;
        this->fFitRangeMax = fitRangeMax;
        this->fNPars = {0};  // The first parameter has index zero
        this->fAncestorIdx = -2;

        this->fFuncMap.insert({"pol0", Pol0});
        this->fFuncMap.insert({"pol1", Pol1});
        this->fFuncMap.insert({"pol2", Pol2});
        this->fFuncMap.insert({"pol3", Pol3});
        this->fFuncMap.insert({"pol4", Pol4});
        this->fFuncMap.insert({"pol5", Pol5});
        this->fFuncMap.insert({"gaus", Gaus});
        this->fFuncMap.insert({"bw", BreitWigner});
        this->fFuncMap.insert({"voigt", Voigt});
        this->fFuncMap.insert({"ComplexLednicky_Singlet_doublegaussian_lambda", ComplexLednicky_Singlet_doublegaussian_lambda});
        this->fFuncMap.insert({"sillkstar", BreitWignerKStar});
        this->fFuncMap.insert({"spline3", Spline3});
        this->fFuncMap.insert({"spline3range", Spline3Range});
        this->fFuncMap.insert({"powerlaw", PowerLaw});
        this->fFuncMap.insert({"flatpol3", FlatPol3});
        this->fFuncMap.insert({"sillkstar", SillKStar});
    }

    void SetBaselineIdx(double basIdx) {
        this->fBaselineIdx = basIdx; 
    }

    void DrawLegend(TVirtualPad *pad, double lowX, double lowY, double highX, double highY, 
                    std::vector<TString> labels = {}) {
        TLegend *legend = new TLegend(lowX, lowY, highX, highY);
        legend->AddEntry(this->fDataHist, labels[0].Data(), "lp");
        legend->AddEntry(this->fFit, labels[1].Data(), "l");
        for(int iLabel=2; iLabel<labels.size(); iLabel++) {
            if(labels[iLabel].Contains("lambda_flat")) continue;
            legend->AddEntry(this->fFitFuncEval[iLabel-2], labels[iLabel].Data(), "l");
        }
        legend->SetBorderSize(0);
        legend->SetTextSize(0.045);
        legend->Draw("same");
        pad->Update();
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

    void Add(TString name, std::vector<std::tuple<std::string, double, double, double>> pars, std::string addmode, bool onbaseline) {
        if(this->fFuncMap.find(name)!=this->fFuncMap.end()){
            this->fFitFunc.push_back(this->fFuncMap[name]);
            this->fFitFuncComps.push_back(name);
        } else {
            printf("Error: function '%s' is not implemented. Exit!", name.Data());
            exit(1);
        }

        this->fAddModes.push_back(addmode);
        this->fDrawOnBaseline.push_back(onbaseline);
        this->fNPars.push_back(pars.size());
        // Save fit settings
        for (const auto &par : pars) {
            this->fFitPars.insert({this->fFitPars.size(), par});
        }
    }

    void AddSplineHisto(TString name, TH1* splinedhisto, std::vector<std::tuple<std::string, double, double, double>> pars, std::string addmode, bool onbaseline,
                        TString legend="ciaospl") {

        TH1D *splineHisto = static_cast<TH1D*>(splinedhisto);
        TSpline3* sp3 = new TSpline3(splinedhisto);
        this->fFitSplines.push_back(sp3);
        this->fFitFuncComps.push_back(name);
        
        this->fAddModes.push_back(addmode);
        this->fDrawOnBaseline.push_back(onbaseline);
        cout << "SPLINE NUMBER OF PARAMETERS: " << pars.size() << endl;
        this->fNPars.push_back(pars.size());
        // Save fit settings
        for (const auto &par : pars) {
            this->fFitPars.insert({this->fFitPars.size(), par});
        }
    }

    void PrefitComponent(TVirtualPad *pad, TH1 *histo, TH1 *histoSpline, TString name, int startTotPar, int nPars, 
                         double lowFitRange, double uppFitRange, double lowReject=0, double uppReject=0) {
        if(lowReject < 0.0001) {
            lowReject = uppFitRange;
        }
        if(uppReject < 0.0001) {
            uppReject = lowFitRange;
        }
        
        TH1D *splineHisto = static_cast<TH1D*>(histoSpline);
        splineHisto->Rebin(2);
        TSpline3* spSigma1385 = new TSpline3(histoSpline);

        TFile *sillFit = new TFile("~/an/LPi/SillFit.root", "recreate");
        // histoSpline->Scale(1/histoSpline->Integral());
        // cout << "Component name: " << name << endl;
        // TF1 *sillSigma1385 = new TF1("sillSigma1385", this->fFuncMap["sillkstar"], 0, 6000, 3);
        // cout << "Parameter initialization for fit" << endl;
        // sillSigma1385->SetParName(0, "normsill");
        // sillSigma1385->FixParameter(0, 3);
        // //sillSigma1385->SetParLimits(0, 0, 100000000);
        // cout << "0: " << sillSigma1385->GetParName(0) << " " << sillSigma1385->GetParameter(0) << endl;
        // sillSigma1385->SetParName(1, "widthsill");
        // sillSigma1385->FixParameter(1, 5);
        // //sillSigma1385->SetParLimits(1, 0, 600000);
        // cout << "1: " << sillSigma1385->GetParName(1) << " " << sillSigma1385->GetParameter(1) << endl;
        // sillSigma1385->SetParName(2, "masssill");
        // sillSigma1385->FixParameter(2, 1385);
        // //sillSigma1385->SetParLimits(2, 0, 10000);
        // cout << "2: " << sillSigma1385->GetParName(2) << " " << sillSigma1385->GetParameter(2) << endl;
        // cout << "Eval: " << sillSigma1385->Eval(100) << endl;
        
        // int statusSillFit = splineHisto->Fit(sillSigma1385, "SMR+0", "")->Status();
        // sillFit->cd();
        TCanvas *cSillFit = new TCanvas("cSillFit", "cSillFit", 600, 600);
        splineHisto->Draw();
        spSigma1385->Draw("same");
        //sillSigma1385->Draw("same");
        cSillFit->Write();
        sillFit->Close();

        cout << "Component name: " << name << endl;
        TF1 *prefitCompFunc = new TF1(name,
            [&, this](double *x, double *pars) -> double {
                if (x[0] >= lowReject && x[0] <= uppReject) {
                    TF1::RejectPoint();
                    //return this->fFuncMap[name](x, pars) + pars[nPars]*spSigma1385->Eval(x[0]);                    
                    return this->fFuncMap[name](x, pars) + pars[nPars]*this->fFuncMap["sillkstar"](x, &pars[nPars+1]);
                } else if (x[0] < lowFitRange) {
                    TF1::RejectPoint();
                    //return this->fFuncMap[name](x, pars) + pars[nPars]*spSigma1385->Eval(x[0]);                    
                    return this->fFuncMap[name](x, pars) + pars[nPars]*this->fFuncMap["sillkstar"](x, &pars[nPars+1]);
                } else {
                    //return this->fFuncMap[name](x, pars) + pars[nPars]*spSigma1385->Eval(x[0]);                    
                    return this->fFuncMap[name](x, pars) + pars[nPars]*this->fFuncMap["sillkstar"](x, &pars[nPars+1]);
                }
            //}, 0, uppFitRange, nPars+4);
            }, 0, uppFitRange, nPars+1);

        for (size_t iPar = startTotPar; iPar < startTotPar + nPars; iPar++) {
            double lowParLimit;
            double uppParLimit;
            this->fFit->GetParLimits(iPar, lowParLimit, uppParLimit);
            cout << "Setting parameters for prefit component" << endl;
            cout << "iPar" << iPar << ": " << this->fFit->GetParName(iPar) << " " << this->fFit->GetParameter(iPar) 
                 << " " << lowParLimit << " " << uppParLimit << endl;
            prefitCompFunc->SetParName(iPar, this->fFit->GetParName(iPar));
            prefitCompFunc->SetParameter(iPar, this->fFit->GetParameter(iPar));
            prefitCompFunc->SetParLimits(iPar, lowParLimit, uppParLimit);
        }

        //prefitCompFunc->SetParName(nPars, "normsill");
        //prefitCompFunc->SetParameter(nPars, 0.000001);
        //prefitCompFunc->SetParLimits(nPars, 0.00000001, 0.0001);
        //cout << "Setting parameters for prefit component" << endl;
        //cout << "iPar" << nPars << ": " << prefitCompFunc->GetParName(nPars) << " " << prefitCompFunc->GetParameter(nPars) << endl;

        //for (size_t iPar = 1; iPar < prefitCompFunc->GetNpar() - nPars; iPar++) {
        //    double lowParLimit;
        //    double uppParLimit;
        //    sillSigma1385->GetParLimits(iPar-1, lowParLimit, uppParLimit);
        //    cout << "Setting parameters for prefit component" << endl;
        //    cout << "iPar" << iPar << ": " << sillSigma1385->GetParName(iPar-1) << " " << sillSigma1385->GetParameter(iPar-1) << endl;
        //    prefitCompFunc->SetParName(iPar+nPars, sillSigma1385->GetParName(iPar-1));
        //    prefitCompFunc->SetParameter(iPar+nPars, sillSigma1385->GetParameter(iPar-1));
        //    prefitCompFunc->SetParLimits(iPar+nPars, lowParLimit, uppParLimit);
        //}
    
        cout << endl;

        for (size_t iPar = 0; iPar < prefitCompFunc->GetNpar(); iPar++) {
            double lowParLimit;
            double uppParLimit;
            prefitCompFunc->GetParLimits(iPar, lowParLimit, uppParLimit);
            cout << "Parameters for prefit" << endl;
            cout << "iPar" << iPar << ": " << prefitCompFunc->GetParName(iPar) << " " << prefitCompFunc->GetParameter(iPar) 
                 << " " << lowParLimit << " " << uppParLimit << endl;
        }

        cout << "Setting parameters for prefit component" << endl;
        cout << "iPar" << nPars+1 << ": " << "splinenorm" << " " << 0.000001 
             << " " << 0.00000001 << " " << 0.00001 << endl;
        prefitCompFunc->SetParName(nPars+1, "splinenorm");
        prefitCompFunc->SetParLimits(nPars+1, 0.00000001, 0.00001);

        //fMCHist->Rebin(2);
        int status = fMCHist->Fit(prefitCompFunc, "SMR+0", "")->Status();

        for (size_t iPar = 0; iPar < nPars; iPar++) {
            double lowParLimit;
            double uppParLimit;
            prefitCompFunc->GetParLimits(iPar, lowParLimit, uppParLimit);
            cout << "Parameter initialization for fit" << endl;
            cout << "iPar" << iPar + startTotPar << ": " << prefitCompFunc->GetParName(iPar) << " " << prefitCompFunc->GetParameter(iPar) << endl;
            this->fFit->SetParName(iPar + startTotPar, prefitCompFunc->GetParName(iPar));
            this->fFit->FixParameter(iPar + startTotPar, prefitCompFunc->GetParameter(iPar));
        }

        TF1 *prefitPol3 = new TF1("pol3",
            [&, this](double *x, double *pars) -> double {
                return this->fFuncMap["pol3"](x, pars);
            }, 0, uppFitRange, nPars+1);

        prefitPol3->SetParameter(0, prefitCompFunc->GetParameter(0));
        prefitPol3->SetParameter(1, prefitCompFunc->GetParameter(1));
        prefitPol3->SetParameter(2, prefitCompFunc->GetParameter(2));
        prefitPol3->SetParameter(3, prefitCompFunc->GetParameter(3));

        TF1 *prefitSplinedHisto = new TF1("splinedhisto",
            [&, this](double *x, double *pars) -> double {
                //return prefitPol3->Eval(x[0]) + prefitCompFunc->GetParameter(nPars)*spSigma1385->Eval(x[0]);
                return 0.9 + prefitCompFunc->GetParameter(nPars)*spSigma1385->Eval(x[0]);
            }, 0, uppFitRange, 0);

        pad->cd();
        double yMinDraw = 0;
        double yMaxDraw = 1.3 * histo->GetMaximum();
        gPad->DrawFrame(0, yMinDraw, uppFitRange, yMaxDraw, ";k* (GeV/c); C(k*)");
    
        prefitCompFunc->SetNpx(300);
        prefitCompFunc->SetLineColor(kRed);
        prefitCompFunc->SetLineWidth(3);
        prefitCompFunc->Clone()->Draw("same");
        pad->Update();

        prefitSplinedHisto->SetNpx(300);
        prefitSplinedHisto->SetLineColor(kGreen);
        prefitSplinedHisto->SetLineWidth(3);
        prefitSplinedHisto->Clone()->Draw("same");
        pad->Update();

        //sillSigma1385->SetNpx(300);
        //sillSigma1385->SetLineColor(kMagenta);
        //sillSigma1385->SetLineWidth(3);
        //sillSigma1385->Clone()->Draw("same");
        //pad->Update();

        prefitPol3->SetNpx(300);
        prefitPol3->SetLineColor(kBlue);
        prefitPol3->SetLineWidth(3);
        prefitPol3->Clone()->Draw("same");
        pad->Update();

        histo->GetYaxis()->SetRangeUser(yMinDraw, yMaxDraw); 
        histo->SetMarkerSize(0.3);
        histo->SetMarkerStyle(20);
        histo->SetMarkerColor(kBlack);
        histo->SetLineColor(kBlack);
        histo->SetLineWidth(3);
        histo->DrawCopy("same pe");
        pad->Update();
    }

    void BuildFitFunction() {
        cout << "----------- Building the fit function -----------" << endl;
        int totTerms = this->fFitFunc.size() + this->fFitSplines.size() + 1;
        // Build the fit function
        this->fFit = new TF1(
            "fFit",
            [&, this](double *x, double *pars) {
                double result = 0;
                int nTerms = this->fFitFunc.size() + this->fFitSplines.size();
                int nPar = 0;
                int nSplineComp = 0;
                int nFuncComp = 0;
                for (int iTerm = 0; iTerm < nTerms; iTerm++) {
                    if(fFitFuncComps[iTerm].Contains("splinehisto")) {
                        auto func = this->fFitSplines[nSplineComp];
                        if(fAddModes[iTerm] == "mult") {
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
                        if(fAddModes[iTerm] == "mult") {
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

        int baselinePar = accumulate(fNPars.begin(), std::next(fNPars.begin(), fBaselineIdx), 0);
        //cout << "Total number of parameters " << this->fFitPars.size() << endl;
        //cout << "Baseline par " << baselinePar << endl;
        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            auto pars = this->fFitPars[iPar];
            cout << "iPar" << iPar << ": " << std::get<0>(pars) << " " << std::get<1>(pars) << " " << std::get<2>(pars) << " " << std::get<3>(pars) << endl;        
            //if(iPar > baselinePar+1) {
            //    cout << "Prefit function iPar " << iPar-baselinePar << ": " << fPrefitMC->GetParameter(iPar-baselinePar-2) << endl;
            //    cout << "Function iPar " << iPar << endl;
            //    this->fFit->SetParName(iPar, std::get<0>(pars).data());
            //    this->fFit->FixParameter(iPar, fPrefitMC->GetParameter(iPar-baselinePar-2));
            //} else {
                //cout << "Parameter name: " << std::get<0>(pars).data() << endl; 
            this->fFit->SetParName(iPar, std::get<0>(pars).data());
            if (std::get<2>(pars) >= std::get<3>(pars)) {
                //cout << "iPar" << iPar << ": " << std::get<1>(pars) << endl;
                this->fFit->FixParameter(iPar, std::get<1>(pars));
            } else {
                this->fFit->SetParameter(iPar, std::get<1>(pars));
                //cout << "iPar" << iPar << ": " << std::get<2>(pars) << " " << std::get<3>(pars) << endl;
                this->fFit->SetParLimits(iPar, std::get<2>(pars), std::get<3>(pars));
            }
        }
    }

    // Perform the fit
    TFitResultPtr Fit() {
        cout << "----------- Fit Parameter initialization -----------" << endl;
        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            double lowParLimit;
            double uppParLimit;
            this->fFit->GetParLimits(iPar, lowParLimit, uppParLimit);
            //cout << "Setting parameters for prefit component" << endl;
            cout << "iPar" << iPar << ": " << this->fFit->GetParName(iPar) << " " << this->fFit->GetParameter(iPar) 
                 << " " << lowParLimit << " " << uppParLimit << endl;
        }

        TFitResultPtr fitResults = fDataHist->Fit(this->fFit, "SMR+0", "");
        return fitResults;
    }

    TF1 *GetFunction() {
        return this->fFit;
    } 

    void DrawPrefit(TVirtualPad *pad, double lowRangeUser=0.0, double uppRangeUserMult=1.05, 
              std::string title=";k* (MeV/c);C(k*)") {

        pad->cd();
        double yMinDraw = lowRangeUser;
        double yMaxDraw = uppRangeUserMult * fDataHist->GetMaximum();
        
        gPad->DrawFrame(fFitRangeMin, yMinDraw, 1000, yMaxDraw, title.data());
        this->fPrefitMC->SetNpx(300);
        this->fPrefitMC->SetLineColor(kRed);
        this->fPrefitMC->SetLineWidth(3);
        this->fPrefitMC->DrawF1(fFitRangeMin+1,1000,"same");
        pad->Update();

        fMCHist->GetYaxis()->SetRangeUser(yMinDraw, yMaxDraw); 
        fMCHist->SetMarkerSize(0.3);
        fMCHist->SetMarkerStyle(20);
        fMCHist->SetMarkerColor(kBlack);
        fMCHist->SetLineColor(kBlack);
        fMCHist->SetLineWidth(3);
        fMCHist->Draw("same pe");
        pad->Update();
    }

    /*
    Define a canvas before calling this function and pass gPad as TVirtualPad
    */
    void Draw(TVirtualPad *pad, std::vector<TString> addComps = {""}, double lowRangeUser=0.0, double uppRangeUserMult=1.05, 
              std::string title=";k* (MeV/c);C(k*)") {

        EvaluateComponents(addComps); 

        pad->cd();
        double yMinDraw = lowRangeUser;
        double yMaxDraw = uppRangeUserMult * fDataHist->GetMaximum();
        
        //gPad->DrawFrame(fFitRangeMin, yMinDraw, fFitRangeMax, yMaxDraw, title.data());
        gPad->DrawFrame(fFitRangeMin, yMinDraw, 2000, yMaxDraw, title.data());
        this->fFit->SetNpx(300);
        this->fFit->SetLineColor(kRed);
        this->fFit->SetLineWidth(3);
        this->fFit->DrawF1(fFitRangeMin+1,1000,"same");
        pad->Update();

        std::vector<Color_t> colors = {kBlue, kAzure + 2, kGreen, kBlue + 2, kOrange, kCyan, kBlack, kMagenta, kYellow};
        for(int iFuncEval=0; iFuncEval<fFitFuncEval.size(); iFuncEval++) {
            this->fFitFuncEval[iFuncEval]->SetNpx(300);
            this->fFitFuncEval[iFuncEval]->SetLineColor(colors[iFuncEval]);
            this->fFitFuncEval[iFuncEval]->SetLineWidth(3);
            this->fFitFuncEval[iFuncEval]->DrawF1(fFitRangeMin+1,1000,"same");
            pad->Update();
        }
        this->fFitFuncEval[3]->Eval(2);

        fDataHist->GetYaxis()->SetRangeUser(yMinDraw, yMaxDraw); 
        fDataHist->SetMarkerSize(0.3);
        fDataHist->SetMarkerStyle(20);
        fDataHist->SetMarkerColor(kBlack);
        fDataHist->SetLineColor(kBlack);
        fDataHist->SetLineWidth(3);
        fDataHist->Draw("same pe");
        pad->Update();
    }

    void Debug() {
        cout << endl;
        cout << endl;
        cout << endl;
        cout << "########## Debugging fit components ##########" << endl;
        cout << "Total spline terms: " << this->fFitSplines.size() << endl;
        cout << "Total func terms: " << this->fFitFunc.size() << endl;
        cout << "--------------------------" << endl;    
        for(int iTerm=0; iTerm<this->fFitSplines.size() + this->fFitFunc.size(); iTerm++) {
            if(iTerm==fBaselineIdx) cout << "BASELINE" << endl;
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

    void EvaluateComponents(std::vector<TString> addComps = {""}) {
        
        int normParNumber = 0;
        for(int iNorm=0; iNorm<fNPars.size()-1; iNorm++) {
            fNorms.push_back(this->fFit->GetParameter(normParNumber));
            normParNumber += fNPars[iNorm+1];
        }

        std::vector<TF1*> components;
        int nTerms = this->fFitFunc.size() + this->fFitSplines.size();
        int setPars = 0;
        int iSpline = 0;
        int iFunc = 0;
        cout << "Fit Function number of parameters: " << this->fFit->GetNpar() << endl;
        for(int iTerm=0; iTerm<nTerms; iTerm++) {
            std::string compName = static_cast<std::string>(this->fFitFuncComps[iTerm]);
            cout << "Set par " << this->fFit->GetParameter(setPars) << endl;
            cout << "Set par +1 " << this->fFit->GetParameter(setPars+1) << endl;
            if(this->fFitFuncComps[iTerm].Contains("spline")) {
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
                components.push_back(new TF1(Form("iComp_%.0f", iTerm), fFitFunc[iFunc], fFitRangeMin, fFitRangeMax, fNPars[iTerm+1]));
                iFunc++;                
                for(int iPar=1; iPar<fNPars[iTerm+1]; iPar++) {
                    components.back()->FixParameter(iPar-1, this->fFit->GetParameter(setPars + iPar));
                }
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
                this->fFitFuncEval.push_back(new TF1(Form("Comp_%i", iFunc),
                    [&, this, iFunc, components](double *x, double *pars) {
                        if(this->fDrawOnBaseline[iFunc]) {
                            return this->fNorms[iFunc]*components[iFunc]->Eval(x[0]) + 
                                   this->fNorms[this->fBaselineIdx]*components[this->fBaselineIdx]->Eval(x[0]);
                        } else if(this->fFitFuncComps[iFunc].Contains("Lednicky")) {
                            return components[iFunc]->Eval(x[0]);
                        } else {
                            return this->fNorms[iFunc]*components[iFunc]->Eval(x[0]);
                        }},
                    fFitRangeMin, fFitRangeMax, 0));            
            }
        }        
            
        if(addComps[0] != "") {
            for(int iAddComp=0; iAddComp<addComps.size(); iAddComp++) {
                int compNumber = this->fFitFuncEval.size();
                this->fFitFuncEval.push_back(new TF1(Form("Comp_%i", compNumber),
                        [&, this, nTerms, addComps, iAddComp, components](double *x, double *pars) {
                        double sum = 0.;
                        bool addBaseline = true;
                        for(int iComp=0; iComp<nTerms; iComp++) {
                            if(addComps[iAddComp].Contains(std::to_string(iComp))) {
                                sum += this->fNorms[iComp]*components[iComp]->Eval(x[0]);
                                if(!this->fDrawOnBaseline[iComp]) addBaseline=false;
                            }
                        }
                        if(addBaseline) {
                            sum += this->fNorms[this->fBaselineIdx]*components[this->fBaselineIdx]->Eval(x[0]);
                        }
                        return sum;}, fFitRangeMin, fFitRangeMax, 0));
            }
        }
    }

    TH1 *fDataHist = nullptr;
    TH1 *fMCHist = nullptr;
    TF1 *fFit = nullptr;
    TF1 *fPrefitMC = nullptr;

    std::map<TString, double (*)(double *x, double *par)> fFuncMap; // Map containing the implemented functions
    std::vector<double (*)(double *x, double *par)> fFitFunc;       // List of function describing each term of the CF model
    double (*fBaseline)(double *x, double *par);                    // Baseline function (for drawing purposes)
    double fBaselineIdx;                                            // Index of baseline function element
    double fAncestorIdx;                                            // Index of baseline function element
    std::vector<TString> fFitFuncComps;                             // Function names of fit components
    std::vector<bool> fDrawOnBaseline;                              // Options to draw fit components on the baseline function
    std::vector<TSpline3 *> fFitSplines;                            // Spline fit components
    std::vector<TF1*> fFitFuncEval;                                 // Fit components evaluated after the fitting
    std::vector<int> fNPars;                                        // Keeps track of how many parameters each function has
    std::vector<std::string> fAddModes;                             // Select mode of adding the contributions to the model
    std::vector<double> fNorms;                                     // Vector saving the norm factor of each term

    double fFitRangeMin;
    double fFitRangeMax;
    std::map<int, std::tuple<std::string, double, double, double>> fFitPars;    // List of fit parameters
    std::map<int, std::tuple<std::string, double, double, double>> fFitNorms;   // List of fit parameters
};

#endif  // FEMPY_CORRELATIONFITTER_HXX_
