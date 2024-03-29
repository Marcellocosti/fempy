#ifndef FEMPY_MODELFITTER_HXX_
#define FEMPY_MODELFITTER_HXX_

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

class ModelFitter {
   public:
    ModelFitter(TH1 *fithist, double fitRangeMin, double fitRangeMax, 
                         double rejectMin=0, double rejectMax=-1) {
        this->fFitHist = reinterpret_cast<TH1 *>(fithist->Clone());
        this->fFitRangeMin = fitRangeMin;
        this->fFitRangeMax = fitRangeMax;
        this->fRejectMin = rejectMin;
        this->fRejectMax = rejectMax;
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
        
        cout << "Adding " << name << endl;
        if(functions.find(name)!=functions.end()){
            this->fFitFunc.push_back(std::get<0>(functions[name]));
            this->fFitFuncComps.push_back(name);
            // -1 needed because pars includes the norm of the term
            if(pars.size()-1 != std::get<1>(functions[name])) {
                printf("Error: wrong number of parameters for function '%s'!", name.Data());
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

    TF1 *GetComponent(int icomp, int ibaseline=-1) {
        if(!this->fFit) {
            throw std::invalid_argument("Fit not performed, component cannot be evaluated!");
        }
        
        if(ibaseline != -1) {
            TF1 *compWithoutNormAndBaseline = new TF1(this->fFitFuncComps[icomp],
                [&, this, icomp, ibaseline](double *x, double *pars) {
                   return ( this->fFitFuncEval[icomp]->Eval(x[0]) - 
                          ( this->fFitFuncEval[ibaseline]->Eval(x[0])/this->fNorms[ibaseline] ) )
                          / this->fNorms[icomp];  
            }, this->fFitRangeMin, this->fFitRangeMax, 0);
            return compWithoutNormAndBaseline;
        } else {
            TF1 *compWithoutNorm = new TF1(this->fFitFuncComps[icomp],
                [&, this, icomp](double *x, double *pars) {
                   return this->fFitFuncEval[icomp]->Eval(x[0]) / this->fNorms[icomp];  
            }, this->fFitRangeMin, this->fFitRangeMax, 0);
            return compWithoutNorm;
        }
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
        cout << "----------- Fit Parameter initialization -----------" << endl;
        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            double lowParLimit;
            double uppParLimit;
            this->fFit->GetParLimits(iPar, lowParLimit, uppParLimit);
            cout << "iPar" << iPar << ": " << this->fFit->GetParName(iPar) << " " << this->fFit->GetParameter(iPar) 
                 << " " << lowParLimit << " " << uppParLimit << endl;
        }

        cout << "Bin content: " << fFitHist->GetBinContent(1) << endl;
        TFitResultPtr fitResults = fFitHist->Fit(this->fFit, "SMR+0", "");
        return fitResults;
    }

    TF1 *GetFitFunction() {
        return this->fFit;
    } 

    /*
    Define a canvas before calling this function and pass gPad as TVirtualPad
    */
    void Draw(TVirtualPad *pad, std::vector<TString> legLabels, std::vector<double> legCoords,
              std::vector<bool> onBaseline, int linesThickness, int basIdx=-1, std::vector<TString> addComps = {""},
              double lowRangeUser=0.0, double uppRangeUser=1.05, std::string title=";k* (MeV/c);C(k*)") {

        EvaluateComponents(basIdx, onBaseline, addComps); 
        cout << "Drawing" << endl;
        pad->cd();
        double yMinDraw = lowRangeUser;
        double yMaxDraw = uppRangeUser + fFitHist->GetMaximum();
        
        cout << "Drawing" << endl;
        TLegend *legend = new TLegend(legCoords[0], legCoords[1], legCoords[2], legCoords[3]);
        legend->AddEntry(this->fFitHist, legLabels[0].Data(), "lp");
        legend->AddEntry(this->fFit, legLabels[1].Data(), "l");

        cout << "Drawing" << endl;
        //gPad->DrawFrame(fFitRangeMin, yMinDraw, fFitRangeMax, yMaxDraw, title.data());
        gPad->DrawFrame(fFitRangeMin, yMinDraw, 2000, yMaxDraw, title.data());
        this->fFit->SetNpx(300);
        this->fFit->SetLineColor(kRed);
        this->fFit->SetLineWidth(linesThickness);
        this->fFit->DrawF1(fFitRangeMin+1,1000,"same");
        pad->Update();
        
        cout << "Drawingciao" << endl;
        std::vector<Color_t> colors = {kMagenta + 3, kAzure + 2, kGreen, kBlue + 2, kOrange, kCyan, kBlack, kGreen+2};
        for(int iFuncEval=0; iFuncEval<fFitFuncEval.size(); iFuncEval++) {
            cout << "Ciao1" << endl;
            this->fFitFuncEval[iFuncEval]->SetNpx(300);
            cout << "Ciao2" << endl;
            this->fFitFuncEval[iFuncEval]->SetLineColor(colors[iFuncEval]); //.data());
            cout << "Ciao3" << endl;
            this->fFitFuncEval[iFuncEval]->SetLineWidth(linesThickness);
            cout << "Ciao4" << endl;
            this->fFitFuncEval[iFuncEval]->DrawF1(fFitRangeMin+1,1000,"same");
            this->fFitFuncEval[iFuncEval]->Eval(100);
            cout << "Ciao5" << endl;
            pad->Update();
            cout << "Ciao6" << endl;
            if(legLabels[iFuncEval+2].Contains("lambda_flat")) continue;
            cout << "Ciao7" << endl;
            legend->AddEntry(this->fFitFuncEval[iFuncEval], legLabels[iFuncEval+2].Data(), "l");
            cout << "Ciao8" << endl;
            cout << endl;
        }
    
        cout << "Drawing" << endl;
        fFitHist->GetYaxis()->SetRangeUser(yMinDraw, yMaxDraw); 
        fFitHist->SetMarkerSize(0.3);
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

    void EvaluateComponents(int basIdx, std::vector<bool> onBaseline, std::vector<TString> addComps = {""}) {
        
        int normParNumber = 0;
        for(int iNorm=0; iNorm<fNPars.size()-1; iNorm++) {
            fNorms.push_back(this->fFit->GetParameter(normParNumber));
            normParNumber += fNPars[iNorm+1];
        }

        std::vector<TF1*> components;
        cout << 'bas' << basIdx << endl;
        int nTerms = this->fFitFunc.size() + this->fFitSplines.size();
        int setPars = 0;
        int iSpline = 0;
        int iFunc = 0;
        cout << "Fit Function terms: " << nTerms << endl;
        cout << "Fit Function number of parameters: " << this->fFit->GetNpar() << endl;
        for(int iTerm=0; iTerm<nTerms; iTerm++) {
            std::string compName = static_cast<std::string>(this->fFitFuncComps[iTerm]);
            cout << "Name of the component: " << compName << endl;
            if(this->fFitFuncComps[iTerm].Contains("spline") && !this->fFitFuncComps[iTerm].Contains("spline3")) {
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
                    cout << "Setting iPar " << iPar << ": " << this->fFit->GetParameter(setPars + iPar) << endl;
                    components.back()->FixParameter(iPar-1, this->fFit->GetParameter(setPars + iPar));
                }
                cout << endl;
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
                cout << "Component: " << this->fFitFuncComps[iFunc] << endl;
                this->fFitFuncEval.push_back(new TF1(this->fFitFuncComps[iFunc],
                    [&, this, iFunc, components, onBaseline, basIdx](double *x, double *pars) {
                        if(onBaseline[iFunc]) {
                            return this->fNorms[iFunc]*components[iFunc]->Eval(x[0]) + 
                                   this->fNorms[basIdx]*components[basIdx]->Eval(x[0]);
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
                        [&, this, nTerms, addComps, iAddComp, components, onBaseline, basIdx](double *x, double *pars) {
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
                        return sum;}, fFitRangeMin, fFitRangeMax, 0));
            }
        }
        cout << "Evaluated" << endl;
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

    double fFitRangeMin;
    double fFitRangeMax;
    double fRejectMin;
    double fRejectMax;
    std::map<int, std::tuple<std::string, double, double, double>> fFitPars;    // List of fit parameters
};

#endif  // FEMPY_MODELFITTER_HXX_
