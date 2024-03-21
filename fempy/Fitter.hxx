#ifndef FEMPY_FITTER_HXX_
#define FEMPY_FITTER_HXX_

#include <map>
#include <string>
#include <tuple>

#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TMath.h"
//#include "FitFunctions.cxx"

class Fitter {
   public:
    Fitter(TH1 *hist, double fitRangeMin, double fitRangeMax) {
        this->fHist = reinterpret_cast<TH1 *>(hist->Clone());
        this->fFitRangeMin = fitRangeMin;
        this->fFitRangeMax = fitRangeMax;
        this->fNPars = {0};  // The first parameter has index zero

        this->fFuncMap.insert({"pol0", Pol0});
        this->fFuncMap.insert({"pol1", Pol1});
        this->fFuncMap.insert({"pol2", Pol2});
        this->fFuncMap.insert({"pol3", Pol3});
        this->fFuncMap.insert({"pol4", Pol4});
        this->fFuncMap.insert({"pol5", Pol5});
        this->fFuncMap.insert({"gaus", Gaus});
        this->fFuncMap.insert({"bw", BreitWigner});
        this->fFuncMap.insert({"voigt", Voigt});
        this->fFuncMap.insert({"landau", Landau});
        this->fFuncMap.insert({"sillkstar", BreitWignerKStar});
        this->fFuncMap.insert({"spline3", Spline3});
        this->fFuncMap.insert({"spline3range", Spline3Range});
        this->fFuncMap.insert({"powerlaw", PowerLaw});
        this->fFuncMap.insert({"flatpol3", FlatPol3});
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
    void Add(TString name, std::vector<std::tuple<std::string, double, double, double>> pars) {
        cout << "ADDING FUNCTION!" << endl;
        if(this->fFuncMap.find(name)!=this->fFuncMap.end()){
            this->fFitFunc.push_back(this->fFuncMap[name]);
            this->fFitFuncComps.push_back(name);
        } else {
            printf("Error: function '%s' is not implemented. Exit!", name.Data());
            exit(1);
        }

        this->fNPars.push_back(pars.size());
        // Save fit settings
        for (const auto &par : pars) {
            this->fFitPars.insert({this->fFitPars.size(), par});
        }
    }

    // Perform the fit
    int Fit(double rejectLow=0, double rejectUpp=0) {
    
        if(rejectLow < 0.0001) {
            rejectLow = 1;
        }
        if(rejectUpp < 0.0001) {
            rejectUpp = 0;
        }
        // Build the fit function
        this->fFit = new TF1(
            "fFit",
            [&, this](double *x, double *pars) {
                double result = 0;
                for (int iTerm = 0; iTerm < this->fFitFunc.size(); iTerm++) {
                    if (x[0] >= rejectLow && x[0] <= rejectUpp) {
                        TF1::RejectPoint();
                    }
                    auto func = this->fFitFunc[iTerm];
                    result += func(x, &pars[this->fNPars[iTerm]]);  // Shift the index of the parameters
                }
                return result;
            },
            fFitRangeMin, fFitRangeMax, this->fFitPars.size());

        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            auto pars = this->fFitPars[iPar];

            this->fFit->SetParName(iPar, std::get<0>(pars).data());
            if (std::get<2>(pars) > std::get<3>(pars)) {
                cout << "iPar: " << iPar << " " << std::get<1>(pars) << endl;
                this->fFit->FixParameter(iPar, std::get<1>(pars));
            } else {
                this->fFit->SetParameter(iPar, std::get<1>(pars));
                cout << "iPar: " << iPar << " " << std::get<2>(pars) << " " << std::get<3>(pars) << endl;
                this->fFit->SetParLimits(iPar, std::get<2>(pars), std::get<3>(pars));
            }
        }

        //return 1;
        //int status = fHist->Fit(this->fFit, "LSMR+0", "")->Status();
        int status = fHist->Fit(this->fFit, "SMR+0", "")->Status();

        return status;
    }

    TF1 *GetFunction() {
        return this->fFit;
    } 

    TF1 *GetComponent(int icomp) {
        //return this->fFitFuncEval[icomp];
        TF1 *funct = new TF1(Form("Comp %.0f", 0), Pol0, fFitRangeMin, fFitRangeMax, fNPars[0]);
        funct->SetParameter(0, this->fFit->GetParameter(0));
        return funct;
    } 

    /*
    Define a canvas before calling this function and pass gPad as TVirtualPad
    */
    void Draw(TVirtualPad *pad) {
        //return;
        pad->cd();
        fHist->GetYaxis()->SetRangeUser(0, 1.3 * fHist->GetMaximum());
        gPad->DrawFrame(fFitRangeMin, 0, fFitRangeMax, 1.3 * fHist->GetMaximum(),
                        Form("%s;%s;%s", this->fHist->GetTitle(), this->fHist->GetXaxis()->GetTitle(),
                             this->fHist->GetYaxis()->GetTitle()));
        
        this->fFit->SetNpx(300);
        this->fFit->SetLineColor(kRed);
        this->fFit->Draw("same");
        pad->Modified();
        pad->Update();

        std::vector<Color_t> colors = {kMagenta + 3, kAzure + 2, kGreen, kBlue + 2, kOrange};
        int setPars = 0;
        std::vector<TF1*> components;
        for(int iFunc=0; iFunc<fFitFunc.size(); iFunc++) {
            cout << "Evaluating component ... " << iFunc << endl;
            components.push_back(new TF1(Form("iComp_%.0f", iFunc), fFitFunc[iFunc], fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]));
            for(int iPar=0; iPar<fNPars[iFunc+1]; iPar++) {
                cout << "N pars: " << fNPars[iFunc+1] << endl;
                cout << "Setting parameter: " << iPar << " at " << this->fFit->GetParameter(setPars + iPar) << endl;
                components.back()->FixParameter(iPar, this->fFit->GetParameter(setPars + iPar));
            }
            cout << "Value of raw component at 1 MeV: " << components.back()->Eval(1) << endl;
            cout << "Multiplicative factor: " << this->fFit->GetParameter(setPars) << endl;
            TString funcName = this->fFitFuncComps[iFunc];
            this->fFitFuncEval.push_back(new TF1(funcName.Data(),
                [&, this, iFunc, setPars, components](double *x, double *pars) {
                    return components[iFunc]->Eval(x[0]);},
                fFitRangeMin, fFitRangeMax, 0));
            cout << "Value of component after normalization at 1 MeV: " << this->fFitFuncEval.back()->Eval(1) << endl;
            setPars += fNPars[iFunc+1];
            cout << endl;
            this->fFitFuncEval[iFunc]->SetNpx(300);
            this->fFitFuncEval[iFunc]->SetLineColor(colors[iFunc]);
            //this->fFitFuncEval.back()->DrawF1(fFitRangeMin+1,fFitRangeMax,"same");
            this->fFitFuncEval.back()->DrawF1(fFitRangeMin,fFitRangeMax,"same");
            pad->Update();
        }
    
        fHist->SetMarkerSize(0.5);
        fHist->SetMarkerStyle(20);
        fHist->SetMarkerColor(kBlack);
        fHist->SetLineColor(kBlack);
        fHist->SetLineWidth(2);
        fHist->Draw("same pe");
        pad->Update();
    }

    /*
    Define a canvas before calling this function and pass gPad as TVirtualPad
    */
    void DrawSpline(TVirtualPad *pad) {
        //return;
        pad->cd();
        fHist->GetYaxis()->SetRangeUser(0, 1.3 * fHist->GetMaximum());
        gPad->DrawFrame(fFitRangeMin, 0, fFitRangeMax, 1.3 * fHist->GetMaximum(),
                        Form("%s;%s;%s", this->fHist->GetTitle(), this->fHist->GetXaxis()->GetTitle(),
                             this->fHist->GetYaxis()->GetTitle()));
        
        fHist->SetMarkerSize(0.5);
        fHist->SetMarkerStyle(20);
        fHist->SetMarkerColor(kBlack);
        fHist->SetLineColor(kBlack);
        fHist->SetLineWidth(2);
        fHist->Draw("same pe");
        printf("---> %d\n", fHist->GetNbinsX());
        pad->Update();

        TSpline3* sp3 = new TSpline3(fHist);
        sp3->SetNpx(300);
        sp3->SetLineColor(kRed);
        sp3->Draw("same");
        pad->Modified();
        pad->Update();
        
    }

    /*
    Define a canvas before calling this function and pass gPad as TVirtualPad
    */
    void DrawComponents(TVirtualPad *pad) {
        pad->cd();
        //fHist->GetYaxis()->SetRangeUser(0, 1.3 * fHist->GetMaximum());
        //gPad->DrawFrame(fFitRangeMin, 0, fFitRangeMax, 1.3 * fHist->GetMaximum(),
        //                Form("%s;%s;%s", this->fHist->GetTitle(), this->fHist->GetXaxis()->GetTitle(),
        //                     this->fHist->GetYaxis()->GetTitle()));
        

        for(int i=0; i<fNPars.size(); i++) {
            // cout << fNPars[i] << endl;
        }
        // cout << "Drawing components... " << endl;
        if(fFitFuncEval.size() > 0) {
            for(int iFunc=0; iFunc<fFitFuncEval.size(); iFunc++) {
                fFitFuncEval[iFunc]->SetNpx(300);
                fFitFuncEval[iFunc]->SetLineColor(kBlack);
                pad->Modified();
                pad->Update();
                //fFitFuncEval[iFunc]->DrawClone("same");
                fFitFuncEval[iFunc]->Draw("same");
                // cout << "Drawn component!" << endl;
            }
        }

        printf("----------------- | %f \n", fFitFuncEval[0]->Eval(0.1));
        fFitFuncEval[0]->Draw("same r");

        printf("---> %d\n", fHist->GetNbinsX());
        pad->Update();
    }

   private:
    TH1 *fHist = nullptr;
    TF1 *fFit = nullptr;

    std::map<TString, double (*)(double *x, double *par)> fFuncMap; // Map containing the implemented functions
    std::vector<double (*)(double *x, double *par)> fFitFunc;       // List of function describing each term of the CF model
    std::vector<TString> fFitFuncComps;       
    std::vector<TF1*> fFitFuncEval;                                 // List of function describing each term of the CF model
    std::vector<int> fNPars;

    double fFitRangeMin;
    double fFitRangeMax;
    std::map<int, std::tuple<std::string, double, double, double>> fFitPars;  // List of fit parameters
};

#endif  // FEMPY_FITTER_HXX_


//    void FitFunctsEval() {
//        int nPar = 0;
//        for(int iFunc=0; iFunc<fFitFunc.size(); iFunc++) {
//            TF1 *iComp = new TF1(Form("Comp %.0f", 0), fFitFunc[iFunc], fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]);
//            for(int iPar=0; iPar<fNPars[iFunc+1]; iPar++) {
//                // cout << this->fFit->GetParameter(nPar + iPar) << endl;
//                iComp->SetParameter(iPar, this->fFit->GetParameter(nPar + iPar));
//            }
//            nPar += fNPars[iFunc+1];
//            // cout << "Evaluated component!" << endl;
//            fFitFuncEval.push_back(iComp);
//        }
//        // cout << fFitFuncEval.size() << endl;
//    }