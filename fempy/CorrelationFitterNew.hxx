#ifndef FEMPY_CORRELATIONFITTERNEW_HXX_
#define FEMPY_CORRELATIONFITTERNEW_HXX_

#include <map>
#include <string>
#include <tuple>

#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TMath.h"
#include "TLegend.h"
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

        this->fFuncMap.insert({"pol0", Pol0});
        this->fFuncMap.insert({"pol1", Pol1});
        this->fFuncMap.insert({"pol2", Pol2});
        this->fFuncMap.insert({"pol3", Pol3});
        this->fFuncMap.insert({"gaus", Gaus});
        this->fFuncMap.insert({"bw", BreitWigner});
        this->fFuncMap.insert({"voigt", Voigt});
        this->fFuncMap.insert({"ComplexLednicky_Singlet_doublegaussian_lambda", ComplexLednicky_Singlet_doublegaussian_lambda});
        this->fFuncMap.insert({"landau", Landau});
        this->fFuncMap.insert({"sillkstar", BreitWignerKStar});
        this->fFuncMap.insert({"spline3", Spline3});
        this->fFuncMap.insert({"spline3range", Spline3Range});
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
        this->fNPars.push_back(pars.size());
        // Save fit settings
        for (const auto &par : pars) {
            this->fFitPars.insert({this->fFitPars.size(), par});
        }
    }

    // Perform the fit
    int Fit() {
        cout << "Fitting" << endl;
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
                // cout << "nSplineComp: " << nSplineComp << "  " << "nFuncComp: " << nFuncComp << endl;
                // cout << "nTerms: " << nTerms << endl;
                for (int iTerm = 0; iTerm < nTerms; iTerm++) {
                    // cout << "nSplineComp: " << nSplineComp << "  " << "nFuncComp: " << nFuncComp << endl;
                    if(fFitFuncComps[iTerm].Contains("splinehisto")) {
                        auto func = this->fFitSplines[nSplineComp];
                        if(fAddModes[iTerm] == "mult") {
                            double partResult = result; 
                            result = pars[nPar]*this->fFitSplines[nSplineComp]->Eval(x[0])*partResult;
                        }
                        else {
                            result += pars[nPar]*this->fFitSplines[nSplineComp]->Eval(x[0]);
                        // cout << "Term number " << iTerm << " added" << endl;             
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
                        // cout << "Term number " << iTerm << " added" << endl;             
                    }
                    nPar += this->fNPars[iTerm+1];   
                }
                return result;}, 
                fFitRangeMin, fFitRangeMax, this->fFitPars.size());

        cout << "Fit function created!" << endl;
        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            auto pars = this->fFitPars[iPar];
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
        
        //cout << "Evaluation of the fit function prefit: " << fFit->Eval(0.1) << endl;
        //int status = fDataHist->Fit(this->fFit, "LSMR+0", "")->Status();
        int status = fDataHist->Fit(this->fFit, "SMR+0", "")->Status();
        //cout << "Evaluation of the fit function after fit: " << fFit->Eval(0.1) << endl;
        
        return status;
    }

    TF1 *GetFunction() {
        return this->fFit;
    } 

    /*
    Define a canvas before calling this function and pass gPad as TVirtualPad
    */
    void Draw(TVirtualPad *pad, std::vector<TString> addComps = {""}, double lowRangeUser=0.0, double uppRangeUserMult=1.05, 
              std::string title=";k* (MeV/c);C(k*)") {

        EvaluateComponentsNew(addComps); 

        cout << "Components evaluated" << endl;

        pad->cd();
        double yMinDraw = lowRangeUser;
        double yMaxDraw = uppRangeUserMult * fDataHist->GetMaximum();
        cout << "yMin, yMax set" << endl;
        
        //fDataHist->GetYaxis()->SetRangeUser(0.98, yMaxDraw);
        // gPad->DrawFrame(fFitRangeMin, 0, fFitRangeMax, yMaxDraw,
        //                 Form("%s;%s;%s", this->fDataHist->GetTitle(), this->fDataHist->GetXaxis()->GetTitle(),
        //                      this->fDataHist->GetYaxis()->GetTitle()));
        gPad->DrawFrame(fFitRangeMin, yMinDraw, fFitRangeMax, yMaxDraw, title.data());
        this->fFit->SetNpx(300);
        this->fFit->SetLineColor(kRed);
        this->fFit->SetLineWidth(3);
        this->fFit->DrawF1(fFitRangeMin+1,fFitRangeMax,"same");
        pad->Update();
        cout << "Fit function drawn" << endl;

        std::vector<Color_t> colors = {kBlue, kAzure + 2, kGreen, kBlue + 2, kOrange, kCyan, kBlack, kMagenta, kYellow};
        cout << "Number of components to be drawn: " << fFitFuncEval.size() << endl;
        //for(int iFuncEval=0; iFuncEval<2; iFuncEval++) { working
        //for(int iFuncEval=0; iFuncEval<3; iFuncEval++) {
        for(int iFuncEval=0; iFuncEval<fFitFuncEval.size(); iFuncEval++) {
            this->fFitFuncEval[iFuncEval]->SetNpx(300);
            this->fFitFuncEval[iFuncEval]->SetLineColor(colors[iFuncEval]);
            this->fFitFuncEval[iFuncEval]->SetLineWidth(3);
            this->fFitFuncEval[iFuncEval]->DrawF1(fFitRangeMin+1,fFitRangeMax,"same");
            cout << "Fit component number " << iFuncEval << " drawn" << endl;
            pad->Update();
        }
        this->fFitFuncEval[3]->Eval(2);

        //exit();


        cout << "CIAO1" << endl;
        fDataHist->GetYaxis()->SetRangeUser(yMinDraw, yMaxDraw); 
        cout << "CIAO2" << endl;
        fDataHist->SetMarkerSize(0.3);
        cout << "CIAO3" << endl;
        fDataHist->SetMarkerStyle(20);
        cout << "CIAO4" << endl;
        fDataHist->SetMarkerColor(kBlack);
        cout << "CIAO5" << endl;
        fDataHist->SetLineColor(kBlack);
        cout << "CIAO6" << endl;
        fDataHist->SetLineWidth(3);
        cout << "CIAO7" << endl;
        fDataHist->Draw("same pe");
        cout << "CIAO8" << endl;
        pad->Update();
        cout << "Histogram drawn" << endl;
    }

    void Debug() {
        cout << "########## Debugging fit components ##########" << endl;
        cout << "Total spline terms: " << this->fFitSplines.size() << endl;
        cout << "Total func terms: " << this->fFitFunc.size() << endl;
        cout << "--------------------------" << endl;    
        for(int iTerm=0; iTerm<this->fFitSplines.size() + this->fFitFunc.size(); iTerm++) {
            if(iTerm==fBaselineIdx) cout << "BASELINE" << endl;
            cout << "Add Mode: " << this->fAddModes[iTerm] << endl;
            cout << "Number of parameters of term " << iTerm << ": " << this->fNPars[iTerm+1] << endl;
            int startPar = accumulate(fNPars.begin(), std::next(fNPars.begin(), iTerm+1), 0);
            for(int iPar=0; iPar<this->fFitSplines.size() + this->fFitFunc.size(); iPar++) {
                cout << "Parameter " << iPar << "value: " << this->fNPars[startPar + iPar] << endl;
            }
            cout << "--------------------------" << endl;
        }
    }

   private:

    void EvaluateComponentsNew(std::vector<TString> addComps = {""}) {
        
        cout << "fNPars size: " << fNPars.size() << endl;
        for(int iNPar = 0; iNPar<fNPars.size(); iNPar++) {
            cout << fNPars[iNPar] << " ";
        }
        cout << endl;
        int normParNumber = 0;
        for(int iNorm=0; iNorm<fNPars.size()-1; iNorm++) {
            cout << "Norm par number: " << normParNumber << endl;
            cout << "Getting parameter number " << normParNumber << endl;
            fNorms.push_back(this->fFit->GetParameter(normParNumber));
            cout << "Norm of " << iNorm << "-th component: " << this->fFit->GetParameter(normParNumber) << endl;
            normParNumber += fNPars[iNorm+1];
            cout << endl;
        }

        cout << endl;
        for(int iPar=0; iPar<this->fFit->GetNpar(); iPar++) {
            cout << "Fit function par " << iPar << "-th: " << this->fFit->GetParameter(iPar) << endl;
        }
        cout << endl;

        std::vector<TF1*> components;
        int nTerms = this->fFitFunc.size() + this->fFitSplines.size();
        cout << "Number of fit components: " << nTerms << endl;
        int setPars = 0;
        int iSpline = 0;
        int iFunc = 0;
        for(int iTerm=0; iTerm<nTerms; iTerm++) {
            std::string compName = static_cast<std::string>(this->fFitFuncComps[iTerm]);
            if(this->fFitFuncComps[iTerm].Contains("spline")) {
                components.push_back(new TF1(Form("iComp_%.0f", iTerm),
                        [&, this, iSpline](double *x, double *pars) {
                            return fFitSplines[iSpline]->Eval(x[0]);
                        }, fFitRangeMin, fFitRangeMax, 0));
                iSpline++;
            } else {
                components.push_back(new TF1(Form("iComp_%.0f", iTerm), fFitFunc[iFunc], fFitRangeMin, fFitRangeMax, fNPars[iTerm+1]));
                iFunc++;                
                for(int iPar=1; iPar<fNPars[iTerm+1]; iPar++) {
                    cout << "Setting parameter: (" << iPar-1 << "," << this->fFit->GetParameter(setPars + iPar) << "), " << endl;
                    components.back()->FixParameter(iPar-1, this->fFit->GetParameter(setPars + iPar));
                }
            }
            setPars += fNPars[iTerm+1];    
            cout << endl;
        }
        
        for(int iFunc=0; iFunc<nTerms; iFunc++) {
            bool toBeSummed = false;
            for(int iAddComp=0; iAddComp<addComps.size(); iAddComp++) {
                
                if(addComps[iAddComp].Contains(std::to_string(iFunc))) {
                    cout << "Component to be summed" << endl;
                    toBeSummed = true;
                    continue;
                }
            }
            
            if(!toBeSummed) {
                cout << "Adding component" << endl;
                cout << "Norm: " << this->fNorms[iFunc] << endl;
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
            
        cout << "FitFuncEval size: " << fFitFuncEval.size() << endl;
        //if(addComps[0] != "") {
        //if(addComps[0] == "3456") {
        cout << addComps[0].Atoi() << endl;                
        if(addComps[0].EqualTo("3456")) {
            for(int iAddComp=0; iAddComp<addComps.size(); iAddComp++) {
                cout << "Combined term " << addComps[iAddComp] << endl;
                int compNumber = this->fFitFuncEval.size();
                cout << "Component number " << compNumber << endl;
                //this->fFitFuncEval.push_back(new TF1("Ciao",
                this->fFitFuncEval.push_back(new TF1(Form("Comp_%i", compNumber),
                        [&, this, nTerms, addComps, iAddComp, components](double *x, double *pars) {
                        //cout << "Ciaoooooooooooooo1" << endl;                
                        double sum = 0.;
                        bool addBaseline = true;
                        //cout << "Ciaoooooooooooooo2" << endl;                
                        for(int iComp=0; iComp<nTerms; iComp++) {
                            //cout << "Ciaoooooooooooooo3" << endl;                
                            //cout << std::to_string(iComp) << endl;                
                            //cout << addComps[iAddComp].Atoi() << endl;                
                            //cout << addComps[iAddComp].Contains("ciao") << endl;                
                            //cout << "Ciaoooooooooooooo4" << endl;                
                            if(addComps[iAddComp].Contains(std::to_string(iComp))) {
                                cout << iComp;
                                sum += this->fNorms[iComp]*components[iComp]->Eval(x[0]);
                                if(!this->fDrawOnBaseline[iComp]) addBaseline=false;
                            }
                            //cout << "Ciaoooooooooooooo4" << endl;                
                        }
                        //cout << "Ciaoooooooooooooo5" << endl;                
                        if(addBaseline) {
                            sum +=  this->fNorms[this->fBaselineIdx]*components[this->fBaselineIdx]->Eval(x[0]);
                        }
                        return sum;}, fFitRangeMin, fFitRangeMax, 0));
            }
        }
        cout << "FitFuncEval size after joint components: " << fFitFuncEval.size() << endl;
    }

    TH1 *fDataHist = nullptr;
    TH1 *fMCHist = nullptr;
    TF1 *fFit = nullptr;

    std::map<TString, double (*)(double *x, double *par)> fFuncMap; // Map containing the implemented functions
    std::vector<double (*)(double *x, double *par)> fFitFunc;       // List of function describing each term of the CF model
    double (*fBaseline)(double *x, double *par);                    // Baseline function (for drawing purposes)
    double fBaselineIdx;                                            // Index of baseline function element
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

#endif  // FEMPY_CORRELATIONFITTERNEW_HXX_



//    EvaluateComponents(std::vector<TString> addComps) {
//        
//        std::vector<TF1*> components;
//        int nTerms = this->fFitFunc.size() + this->fFitSplines.size() + 1;
//        int setPars = 0;
//        for(int iFunc=0; iFunc<nTerms; iFunc++) {
//            if(iFunc==fBaselineIdx) {
//                setPars += fNPars[iFunc+1];
//                continue;
//            }
//            //cout << "Evaluating component ... " << iFunc << endl;
//            TString funcName = this->fFitFuncComps[iFunc];
//            std::string compName = static_cast<std::string>(this->fFitFuncComps[iFunc]);
//            components.push_back(new TF1(Form("iComp_%.0f", iFunc), fFitFunc[nFuncComp], fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]-1));
//            for(int iPar=1; iPar<fNPars[iFunc+1]; iPar++) {
//                cout << "Setting parameter: (" << iPar << "," << this->fFit->GetParameter(setPars + iPar) << "), " << endl;
//                components.back()->FixParameter(iPar-1, this->fFit->GetParameter(setPars + iPar));
//            }
//            //cout << "Name of component: " << funcName << endl;
//            if(fDrawOnBaseline[iFunc]) {
//                //cout << "Drawing spline" << endl;
//                this->fFitFuncEval.push_back(new TF1(compName.data(),
//                    [&, this, setPars, iTerm](double *x, double *pars) {
//                        return this->fFitFuncEval[fBaselineIdx]->Eval(x[0]) + this->fFit->GetParameter(setPars)*this->fFitFunc[iTerm]->Eval(x[0]);
//                    },
//                    fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]-1));
//            } else {
//                this->fFitFuncEval.push_back(new TF1(compName.data(),
//                    [&, this, setPars, iTerm](double *x, double *pars) {
//                        return this->fFit->GetParameter(setPars)*this->fFitFunc[iTerm]->Eval(x[0]);
//                    },
//                    fFitRangeMin, fFitRangeMax, fNPars[iFunc+1]-1));
//            }
//            //cout << "Value of component after normalization at 2 MeV: " << this->fFitFuncEval.back()->Eval(2) << endl;
//            //cout << endl;
//            setPars += fNPars[iFunc+1];    
//        } else {
//                //cout << fNPars[iFunc+1] << endl;
//                
//                cout << "Value of raw component at 2 MeV: " << components.back()->Eval(2) << ", " << "Multiplicative factor: " << this->fFit->GetParameter(setPars) << endl;
//                if(funcName.Contains("Lednicky")) {
//                    this->fFitFuncEval.push_back(new TF1("model",
//                            [&, this, setPars, components](double *x, double *pars) {
//                                return components.back()->Eval(x[0]);},
//                            fFitRangeMin, fFitRangeMax, 0));
//                    
//                    cout << "Value of baseline at 2 MeV: " << this->fFitFuncEval[0]->Eval(2) << endl;
//                    cout << "Value of component after normalization at 2 MeV: " << this->fFitFuncEval.back()->Eval(2) << endl;
//                    //cout << endl;
//                    setPars += fNPars[iFunc+1];
//                    nFuncComp++;
//                } else {
//                    if(fDrawOnBaseline[iFunc]) { 
//                        this->fFitFuncEval.push_back(new TF1(compName.data(),
//                            [&, this, setPars, components](double *x, double *pars) {
//                                return this->fFitFuncEval[0]->Eval(x[0]) + this->fFit->GetParameter(setPars)*components.back()->Eval(x[0]);},
//                            fFitRangeMin, fFitRangeMax, 0));
//                    } else { 
//                        this->fFitFuncEval.push_back(new TF1(compName.data(),
//                            [&, this, setPars, components](double *x, double *pars) {
//                                return this->fFit->GetParameter(setPars)*components.back()->Eval(x[0]);},
//                            fFitRangeMin, fFitRangeMax, 0));
//                    }
//                    //cout << "Value of baseline at 2 MeV: " << this->fFitFuncEval[0]->Eval(2) << endl;
//                    //cout << "Value of component after normalization at 2 MeV: " << this->fFitFuncEval.back()->Eval(2) << endl;
//                    //cout << endl;
//                    setPars += fNPars[iFunc+1];
//                    nFuncComp++;
//                }
//            }
//        }
//