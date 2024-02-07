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

    void Add(TString name, std::vector<std::tuple<std::string, double, double, double>> pars, std::string addmode="") {
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
        } else if (name == "landau") { // Complex Lednicky
            this->fFitFunc.push_back(Landau);
            this->fFitFuncComps.push_back(name);
        } else if (name == "sillkstar") { // // breit wigner kstar dependent
            this->fFitFunc.push_back(BreitWignerKStar);
            this->fFitFuncComps.push_back(name);
        } else if (name == "spline3") { // Spline3
            this->fFitFunc.push_back(Spline3);
            this->fFitFuncComps.push_back(name);
        } else if (name == "spline3range") { // breit wigner kstar dependent
            this->fFitFunc.push_back(Spline3Range);
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

    void AddSplineHisto(TString name, TH1* splinedhisto, std::vector<std::tuple<std::string, double, double, double>> pars, std::string addmode="") {

        this->fFitFuncComps.push_back(name);
        this->fAddModes.push_back(addmode);
        
        this->fNPars.push_back(pars.size());
        // Save fit settings
        for (const auto &par : pars) {
            this->fFitPars.insert({this->fFitPars.size(), par});
        }
        
        TH1D *splineHisto = static_cast<TH1D*>(splinedhisto);
        TSpline3* sp3 = new TSpline3(splinedhisto);
        this->fFitSplines.push_back(sp3);
        cout << "Added spline" << endl;
        cout << endl;
    }

    void DrawSpline(TVirtualPad *pad, TH1* hist) {
        pad->cd();
        //double yMaxDraw = 1.3 * fDataHist->GetMaximum();

        TH1D* histo = static_cast<TH1D*>(hist);
        double yMaxDraw = histo->GetBinContent(histo->GetMaximumBin())*1.2;
        fDataHist->GetYaxis()->SetRangeUser(0, yMaxDraw);
        gPad->DrawFrame(fFitRangeMin, 0, fFitRangeMax, yMaxDraw,
                        Form("%s;%s;%s", this->fDataHist->GetTitle(), this->fDataHist->GetXaxis()->GetTitle(),
                             this->fDataHist->GetYaxis()->GetTitle()));

        TSpline3* sp3 = new TSpline3(histo);
        sp3->SetNpx(300);
        sp3->SetLineColor(kRed);
        sp3->Draw("same");

        //this->fFitSplines.back()->SetNpx(300);
        //this->fFitSplines.back()->SetLineColor(kRed);
        //this->fFitSplines.back()->Draw("same");

        histo->SetMarkerSize(0.5);
        histo->SetMarkerStyle(20);
        histo->SetMarkerColor(kBlack);
        histo->SetLineColor(kBlack);
        histo->SetLineWidth(2);
        histo->Draw("same pe");
        pad->Update();

    }

    // Perform the fit
    int Fit() {
        cout << "Fitting" << endl;
        int totTerms = this->fFitFunc.size() + this->fFitSplines.size();
        cout << "Total terms: " << this->fFitFunc.size() << " + " << this->fFitSplines.size() << endl;
        // Build the fit function
        this->fFit = new TF1(
            "fFit",
            [&, this](double *x, double *pars) {
                double result = 0;
                int nTerms = this->fFitFunc.size() + this->fFitSplines.size();
                int nPar = 0;
                int nSplineComp = 0;
                int nFuncComp = 0;
                //cout << "nSplineComp: " << nSplineComp << "  " << "nFuncComp: " << nFuncComp << endl;
                //cout << "nTerms: " << nTerms << endl;
                //cout << "nSplineComp: " << nSplineComp << "  " << "nFuncComp: " << nFuncComp << endl;
                for (int iTerm = 0; iTerm < nTerms; iTerm++) {

                    if(fFitFuncComps[iTerm].Contains("splinehisto")) {
                        auto func = this->fFitSplines[nSplineComp];
                        if(fAddModes[iTerm] == "mult") {
                            double partResult = result; 
                            result = pars[nPar]*this->fFitSplines[nSplineComp]->Eval(x[0])*partResult;
                        }
                        else {
                            result += pars[nPar]*this->fFitSplines[nSplineComp]->Eval(x[0]);
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
                return result;
            },
            fFitRangeMin, fFitRangeMax, this->fFitPars.size());

        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            auto pars = this->fFitPars[iPar];
            this->fFit->SetParName(iPar, std::get<0>(pars).data());
            if (std::get<2>(pars) >= std::get<3>(pars)) {
                cout << "iPar" << iPar << ": " << std::get<1>(pars) << endl;
                this->fFit->FixParameter(iPar, std::get<1>(pars));
            } else {
                this->fFit->SetParameter(iPar, std::get<1>(pars));
                cout << "iPar" << iPar << ": " << std::get<2>(pars) << " " << std::get<3>(pars) << endl;
                this->fFit->SetParLimits(iPar, std::get<2>(pars), std::get<3>(pars));
            }
        }
        
        cout << "Evaluation of the fit function prefit: " << fFit->Eval(0.1) << endl;
        //int status = fDataHist->Fit(this->fFit, "LSMR+0", "")->Status();
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
        //double yMaxDraw = 1.3 * fDataHist->GetMaximum();
        double yMaxDraw = 1.15;

        fDataHist->GetYaxis()->SetRangeUser(0, yMaxDraw);
        gPad->DrawFrame(fFitRangeMin, 0, fFitRangeMax, yMaxDraw,
                        Form("%s;%s;%s", this->fDataHist->GetTitle(), this->fDataHist->GetXaxis()->GetTitle(),
                             this->fDataHist->GetYaxis()->GetTitle()));
        
        this->fFit->SetNpx(300);
        this->fFit->SetLineColor(kRed);
        this->fFit->DrawF1(fFitRangeMin+1,fFitRangeMax,"same");
        pad->Update();

        cout << endl << endl;

        std::vector<Color_t> colors = {kMagenta + 3, kAzure + 2, kGreen, kBlue + 2, kOrange};
        std::vector<TF1*> components;
        
        int nTerms = this->fFitFunc.size() + this->fFitSplines.size();
        int nSplineComp = 0;
        int nFuncComp = 0;
        int setPars = 0;
        cout << "N pars: ";
        for(int iNPars=0; iNPars<nTerms+1; iNPars++) {
            cout << fNPars[iNPars] << " " << endl;
        }
        for(int iFunc=0; iFunc<nTerms; iFunc++) {
            //if(iFunc==1) {setPars += fNPars[iFunc+1]; continue;}
            cout << "Evaluating component ... " << iFunc << endl;
            TString funcName = this->fFitFuncComps[iFunc];
            std::string compName = static_cast<std::string>(this->fFitFuncComps[iFunc]);
            cout << "Name of component: " << funcName << endl;

            if(funcName.Contains("splinehisto")) {
                this->fFitFuncEval.push_back(new TF1(compName.data(),
                    [&, this, setPars, nSplineComp](double *x, double *pars) {
                        return this->fFit->GetParameter(setPars)*this->fFitSplines[nSplineComp]->Eval(x[0]);
                    },
                    fFitRangeMin, fFitRangeMax, 0));
                
                
                this->fFitFuncEval.back()->SetNpx(300);
                this->fFitFuncEval.back()->SetLineColor(colors[iFunc]);
                this->fFitFuncEval.back()->Draw("same");
                pad->Update();
                cout << "Value of component after normalization at 2 MeV: " << this->fFitFuncEval.back()->Eval(2) << endl;
                cout << endl;
                setPars += fNPars[iFunc+1];
                nSplineComp++;
            } else {
                cout << "Ciao" << endl;
                cout << fNPars[iFunc+1] << endl;
                components.push_back(new TF1(Form("iComp_%.0f", iFunc), fFitFunc[nFuncComp], fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]-1));
                for(int iPar=1; iPar<fNPars[iFunc+1]; iPar++) {
                    cout << "N pars: " << fNPars[iFunc+1] << endl;
                    cout << "Setting parameter: " << iPar << " at " << this->fFit->GetParameter(setPars + iPar) << endl;
                    components.back()->FixParameter(iPar-1, this->fFit->GetParameter(setPars + iPar));
                }
                cout << "Value of raw component at 2 MeV: " << components.back()->Eval(2) << endl;
                cout << "Multiplicative factor: " << this->fFit->GetParameter(setPars) << endl;
                this->fFitFuncEval.push_back(new TF1(compName.data(),
                    [&, this, nFuncComp, setPars, components](double *x, double *pars) {
                        return this->fFit->GetParameter(setPars)*components[nFuncComp]->Eval(x[0]);},
                    fFitRangeMin, fFitRangeMax, 0));
                cout << "Value of component after normalization at 2 MeV: " << this->fFitFuncEval[nFuncComp]->Eval(2) << endl;
                setPars += fNPars[iFunc+1];
                cout << endl;
                this->fFitFuncEval.back()->SetNpx(300);
                this->fFitFuncEval.back()->SetLineColor(colors[iFunc]);
                this->fFitFuncEval.back()->DrawF1(fFitRangeMin+1,fFitRangeMax,"same");
                pad->Update();
                nFuncComp++;
            }
        }
        
        printf("size %d\n", fFitFuncEval.size());
        
        for(int iEvalFunc=0; iEvalFunc<fFitFuncEval.size(); iEvalFunc++) {
            cout << "Normalized component " << iEvalFunc << " at 2 MeV: " << this->fFitFuncEval[iEvalFunc]->Eval(200) << endl;
        }
        cout << endl;
 
        fDataHist->SetMarkerSize(0.5);
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
    std::vector<TString> fFitFuncComps;                             // List of function describing each term of the CF model
    std::vector<TSpline3 *> fFitSplines;                            // List of function describing each term of the CF model
    std::vector<TF1*> fFitFuncEval;                                 // List of function describing each term of the CF model
    std::vector<int> fNPars;                                        // Keeps track of how many parameters each function has
    std::vector<std::string> fAddModes;                             // select mode of adding the contributions

    double fFitRangeMin;
    double fFitRangeMax;
    std::map<int, std::tuple<std::string, double, double, double>> fFitPars;        // List of fit parameters
};

#endif  // FEMPY_CORRELATIONFITTERNEW_HXX_


                        //if(iFunc == 0) {
                        //    return this->fFit->GetParameter(fNPars[1])*components[0]->Eval(x[0]); // components[1]->Eval(x[0]);
                        //} 
                        //if(iFunc > 2) {
                        //    return this->fFit->GetParameter(fNPars[1] + fNPars[2])*components[2]->Eval(x[0]) + this->fFit->GetParameter(setPars)*components[iFunc]->Eval(x[0]);
                        //} else{
                        //    return this->fFit->GetParameter(setPars)*components[iFunc]->Eval(x[0]);
                        //}
                    //},

//                    if(fAddModes[iTerm] == "mult") {
//                        double partResult = result; 
//                        if(fFitFuncComps[iTerm]=="splinehisto") {
//                            result = pars[nPar]*sp3->Eval(x[0])*partResult;
//                        }
//                        else {
//                            result = pars[nPar]*func(x, &pars[nPar+1])*partResult;
//                        }
//                    } else {
//                        if(fFitFuncComps[iTerm]=="splinehisto") {
//                            result = pars[nPar]*sp3->Eval(x[0]);
//                        } else {
//                            result += pars[nPar]*func(x, &pars[nPar+1]);
//                        }
//                    }


//    TF1 *GetComponent(int icomp) {
//        //return this->fFitFuncEval[icomp];
//        TF1 *funct = new TF1(Form("Comp %.0f", 0), Pol0, fFitRangeMin, fFitRangeMax, fNPars[0]);
//        funct->SetParameter(0, this->fFit->GetParameter(0));
//        return funct;
//    } 

        //for(int iFunc=0; iFunc<fFitFuncEval.size(); iFunc++) {
        //    //cout << "Drawing function" << endl;
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
        //    //cout << "Evaluating component ..." << endl;
        //    for(int iPar=0; iPar<fNPars[iFunc+1]; iPar++) {
        //        //cout << "iPar: " << iPar << " Value: " << this->fFit->GetParameter(setPars + iPar) << endl;
        //    }
        //    setPars += fNPars[iFunc+1];
        //    //cout << endl;
        //}


        //cout << "Ciao6" << endl;
        //for(int iFuncEval=0; iFuncEval<fFitFuncEval.size(); iFuncEval++) {
        //    for(int iPar=0; iPar<fFitFuncEval[iFuncEval]->GetNpar(); iPar++) {
        //        //cout << "Setted parameter: " << iPar << " at " << this->fFitFuncEval[iFuncEval]->GetParameter(iPar) << endl;
        //    }
        //}
        
//        int setPars = 0;
//        for(int iFunc=0; iFunc<fFitFunc.size(); iFunc++) {
//            //cout << "Evaluating component ..." << endl;
//            this->fFitFuncEval.push_back(new TF1(Form("fComp_%.0f", iFunc), fFitFunc[iFunc], 
//                                                      fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]-1));
//            //TF1 *comp = new TF1(Form("fComp_%.0f", iFunc), fFitFunc[iFunc], fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]);
//            //cout << "N pars: " << fNPars[iFunc+1]-1 << endl;
//            for(int iPar=0; iPar<fNPars[iFunc+1]-1; iPar++) {
//                //cout << "Setting parameter: " << iPar << " at " << this->fFit->GetParameter(setPars + iPar + 1) << endl;
//                this->fFitFuncEval.back()->SetParameter(iPar, this->fFit->GetParameter(setPars + iPar + 1));
//            }
//            //this->fFitFuncEval.push_back(new TF1(Form("fComp_%.0f", iFunc), this->fFit->GetParameter(setPars)*comp->Eval(x[0]), 
//            //                                          fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]));
//            setPars += fNPars[iFunc+1];
//            //cout << "Evaluation of the fit function after fit: " << fFitFuncEval.back()->Eval(0.1) << endl;
//            //cout << "Evaluated component!" << endl;
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