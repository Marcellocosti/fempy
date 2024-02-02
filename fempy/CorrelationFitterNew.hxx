#ifndef FEMPY_CORRELATIONFITTERNEW_HXX_
#define FEMPY_CORRELATIONFITTERNEW_HXX_

#include <map>
#include <string>
#include <tuple>

#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TMath.h"
#include "LednickyLambdaPion.cxx"
#include "FitFunctions.cxx"

class CorrelationFitterNew {
   public:
    CorrelationFitterNew(TH1 *datahist, TH1 *mchist, double fitRangeMin, double fitRangeMax) {
        this->fDataHist = reinterpret_cast<TH1 *>(datahist->Clone());
        this->fMCHist = reinterpret_cast<TH1 *>(mchist->Clone());
        this->fFitRangeMin = fitRangeMin;
        this->fFitRangeMax = fitRangeMax;
        this->fNPars = {0};  // The first parameter has index zero
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

    void Add(std::string name, std::vector<std::tuple<std::string, double, double, double>> pars, std::string addmode="") {
        //cout << "Number of parameters: ";
        //cout << pars.size() << endl;

        fAddModes.push_back(addmode);

        if (name == "pol0") {  // Constant function
            this->fFitFunc.push_back(Pol0);
            this->fFitFuncComps.push_back(name);
        } else if (name == "pol1") {  // pol1
            this->fFitFunc.push_back(Pol1);
            this->fFitFuncComps.push_back(name);
        } else if (name == "pol2") {  // pol2
            this->fFitFunc.push_back(Pol2);
            this->fFitFuncComps.push_back(name);
        } else if (name == "pol3") {  // pol3
            this->fFitFunc.push_back(Pol3);
            this->fFitFuncComps.push_back(name);
        } else if (name == "gaus") {  // Breit-Wigner function
            this->fFitFunc.push_back(Gaus);
            this->fFitFuncComps.push_back(name);
        } else if (name == "bw") {    // Breit-Wigner function
            this->fFitFunc.push_back(BreitWigner);
            this->fFitFuncComps.push_back(name);
        } else if (name == "voigt") { // Breit-Wigner function
            this->fFitFunc.push_back(Voigt);
            this->fFitFuncComps.push_back(name);
        } else if (name == "ComplexLednicky_Singlet_doublegaussian_lambda") { // Complex Lednicky
            this->fFitFunc.push_back(ComplexLednicky_Singlet_doublegaussian_lambda);
            this->fFitFuncComps.push_back(name);
        } else if (name == "spline3") { // Spline3
            this->fFitFunc.push_back(Spline3);
            this->fFitFuncComps.push_back(name);
        } else {
            //printf("Error: function '%s' is not implemented. Exit!", name.data());
            exit(1);
        }

        this->fNPars.push_back(pars.size());

        // Save fit settings
        for (const auto &par : pars) {
            this->fFitPars.insert({this->fFitPars.size(), par});
        }
    }

    // Perform the fit
    int Fit() {
        
        // Build the fit function
        this->fFit = new TF1(
            "fFit",
            [&, this](double *x, double *pars) {
                double result = 0;
                int nPar = 0;
                for (int iTerm = 0; iTerm < this->fFitFunc.size(); iTerm++) {
                    auto func = this->fFitFunc[iTerm];
                    if(fAddModes[iTerm] == "mult") {
                        double partResult = result; 
                        result = pars[nPar]*func(x, &pars[nPar+1])*partResult;
                    } else {
                        result += pars[nPar]*func(x, &pars[nPar+1]);
                    }
                    nPar += this->fNPars[iTerm+1];
                }
                return result;
            },
            fFitRangeMin, fFitRangeMax, this->fFitPars.size());

        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            auto pars = this->fFitPars[iPar];
            this->fFit->SetParName(iPar, std::get<0>(pars).data());
            if (std::get<2>(pars) > std::get<3>(pars)) {
                //cout << "iPar: " << iPar << " " << std::get<1>(pars) << endl;
                this->fFit->FixParameter(iPar, std::get<1>(pars));
            } else {
                this->fFit->SetParameter(iPar, std::get<1>(pars));
                //cout << "iPar: " << iPar << " " << std::get<2>(pars) << " " << std::get<3>(pars) << endl;
                this->fFit->SetParLimits(iPar, std::get<2>(pars), std::get<3>(pars));
            }
        }

        //cout << "Evaluation of the fit function prefit: " << fFit->Eval(0.1) << endl;
        int status = fDataHist->Fit(this->fFit, "SMR+0", "")->Status();
        cout << "Evaluation of the fit function after fit: " << fFit->Eval(0.1) << endl;
        
        return status;
    }

    TF1 *GetFunction() {
        return this->fFit;
    } 

    /*
    Define a canvas before calling this function and pass gPad as TVirtualPad
    */
    void Draw(TVirtualPad *pad) {
        pad->cd();
        fDataHist->GetYaxis()->SetRangeUser(0, 1.3 * fDataHist->GetMaximum());
        gPad->DrawFrame(fFitRangeMin, 0, fFitRangeMax, 1.3 * fDataHist->GetMaximum(),
                        Form("%s;%s;%s", this->fDataHist->GetTitle(), this->fDataHist->GetXaxis()->GetTitle(),
                             this->fDataHist->GetYaxis()->GetTitle()));
        
        this->fFit->SetNpx(300);
        this->fFit->SetLineColor(kRed);
        this->fFit->DrawF1(fFitRangeMin,fFitRangeMax,"same");
        pad->Update();

        cout << endl << endl;

        std::vector<Color_t> colors = {kMagenta + 3, kAzure + 2, kGreen, kBlue + 2};
        int setPars = 0;
        std::vector<TF1*> components;
        for(int iFunc=0; iFunc<fFitFunc.size(); iFunc++) {
            cout << "Evaluating component ... " << iFunc << endl;
            components.push_back(new TF1(Form("iComp_%.0f", iFunc), fFitFunc[iFunc], fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]-1));
            for(int iPar=0; iPar<fNPars[iFunc+1]-1; iPar++) {
                cout << "N pars: " << fNPars[iFunc+1] << endl;
                cout << "Setting parameter: " << iPar << " at " << this->fFit->GetParameter(setPars + iPar + 1) << endl;
                components.back()->FixParameter(iPar, this->fFit->GetParameter(setPars + iPar + 1));
            }
            cout << "Value of raw component at 1 MeV: " << components.back()->Eval(1) << endl;
            cout << "Multiplicative factor: " << this->fFit->GetParameter(setPars) << endl;
            std::string funcName = this->fFitFuncComps[iFunc];
            this->fFitFuncEval.push_back(new TF1(funcName.data(),
                [&, this, iFunc, setPars, components](double *x, double *pars) {
                    return this->fFit->GetParameter(setPars)*components[iFunc]->Eval(x[0]);},
                fFitRangeMin, fFitRangeMax, 0));
            cout << "Value of component after normalization at 1 MeV: " << this->fFitFuncEval.back()->Eval(1) << endl;
            setPars += fNPars[iFunc+1];
            cout << endl;
            this->fFitFuncEval[iFunc]->SetNpx(300);
            this->fFitFuncEval[iFunc]->SetLineColor(colors[iFunc]);
            this->fFitFuncEval.back()->Draw("same");
            pad->Update();
        }
        
        printf("size %d\n", fFitFuncEval.size());
        
        for(int iEvalFunc=0; iEvalFunc<fFitFuncEval.size(); iEvalFunc++) {
            cout << "Normalized component " << iEvalFunc << " at 1 MeV: " << this->fFitFuncEval[iEvalFunc]->Eval(1) << endl;
        }
        cout << endl;
 
        fDataHist->SetMarkerSize(1);
        fDataHist->SetMarkerStyle(20);
        fDataHist->SetMarkerColor(kBlack);
        fDataHist->SetLineColor(kBlack);
        fDataHist->SetLineWidth(2);
        fDataHist->Draw("same pe");
        pad->Update();
    }

   private:
    TH1 *fDataHist = nullptr;
    TH1 *fMCHist = nullptr;
    TF1 *fFit = nullptr;

    std::vector<double (*)(double *x, double *par)> fFitFunc;       // List of function describing each term of the CF model
    std::vector<std::string> fFitFuncComps;                         // List of function describing each term of the CF model
    std::vector<TF1*> fFitFuncEval;                                 // List of function describing each term of the CF model
    std::vector<int> fNPars;                                        // Keeps track of how many parameters each function has
    std::vector<std::string> fAddModes;                             // select mode of adding the contributions

    double fFitRangeMin;
    double fFitRangeMax;
    std::map<int, std::tuple<std::string, double, double, double>> fFitPars;        // List of fit parameters
};

#endif  // FEMPY_CORRELATIONFITTERNEW_HXX_



//    TF1 *GetComponent(int icomp) {
//        //return this->fFitFuncEval[icomp];
//        TF1 *funct = new TF1(Form("Comp %.0f", 0), Pol0, fFitRangeMin, fFitRangeMax, fNPars[0]);
//        funct->SetParameter(0, this->fFit->GetParameter(0));
//        return funct;
//    } 

        //for(int iFunc=0; iFunc<fFitFuncEval.size(); iFunc++) {
        //    cout << "Drawing function" << endl;
        //    this->fFitFuncEval[iFunc]->SetNpx(300);
        //    this->fFitFuncEval[iFunc]->SetLineColor(colors[iFunc]);
        //    this->fFitFuncEval[iFunc]->Draw("same");
        //    pad->Update();
        //}

        //cout << "NUMBER OF PARAMETERS: " << endl;
        //cout << fNPars[1] << " " << fNPars[2] << " " << fNPars[3] << " " << fNPars[4] << endl;
        //cout << endl;
        //int setPars = 0;
        //for(int iFunc=0; iFunc<fFitFunc.size(); iFunc++) {
        //    cout << "Evaluating component ..." << endl;
        //    for(int iPar=0; iPar<fNPars[iFunc+1]; iPar++) {
        //        cout << "iPar: " << iPar << " Value: " << this->fFit->GetParameter(setPars + iPar) << endl;
        //    }
        //    setPars += fNPars[iFunc+1];
        //    cout << endl;
        //}


        //cout << "Ciao6" << endl;
        //for(int iFuncEval=0; iFuncEval<fFitFuncEval.size(); iFuncEval++) {
        //    for(int iPar=0; iPar<fFitFuncEval[iFuncEval]->GetNpar(); iPar++) {
        //        cout << "Setted parameter: " << iPar << " at " << this->fFitFuncEval[iFuncEval]->GetParameter(iPar) << endl;
        //    }
        //}
        
//        int setPars = 0;
//        for(int iFunc=0; iFunc<fFitFunc.size(); iFunc++) {
//            cout << "Evaluating component ..." << endl;
//            this->fFitFuncEval.push_back(new TF1(Form("fComp_%.0f", iFunc), fFitFunc[iFunc], 
//                                                      fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]-1));
//            //TF1 *comp = new TF1(Form("fComp_%.0f", iFunc), fFitFunc[iFunc], fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]);
//            cout << "N pars: " << fNPars[iFunc+1]-1 << endl;
//            for(int iPar=0; iPar<fNPars[iFunc+1]-1; iPar++) {
//                cout << "Setting parameter: " << iPar << " at " << this->fFit->GetParameter(setPars + iPar + 1) << endl;
//                this->fFitFuncEval.back()->SetParameter(iPar, this->fFit->GetParameter(setPars + iPar + 1));
//            }
//            //this->fFitFuncEval.push_back(new TF1(Form("fComp_%.0f", iFunc), this->fFit->GetParameter(setPars)*comp->Eval(x[0]), 
//            //                                          fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]));
//            setPars += fNPars[iFunc+1];
//            //cout << "Evaluation of the fit function after fit: " << fFitFuncEval.back()->Eval(0.1) << endl;
//            cout << "Evaluated component!" << endl;
//            //delete comp;
//        }

            //auto fff = new TF1(
            //    Form("fComp_%.0f", iFunc),
            //    [&, this, iFunc, setPars, components](double *x, double *pars) {
            //        return this->fFit->GetParameter(setPars)*components[iFunc]->Eval(x[0]);
            //    },
            //    fFitRangeMin, fFitRangeMax, 0);
            //
            //this->fFitFuncEval.push_back(fff);