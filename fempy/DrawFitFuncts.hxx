#ifndef FEMPY_DRAWFITFUNCTS_HXX_
#define FEMPY_DRAWFITFUNCTS_HXX_

#include <map>
#include <string>
#include <tuple>
#include <stdexcept>
#include <numeric>

#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TVirtualPad.h"
#include "THashList.h"
#include "TFitResultPtr.h"

class DrawFitFuncts {
   public:
    DrawFitFuncts(TH1 *fithist, double drawRangeMin, double drawRangeMax, 
                  bool globNorm=0, int basIdx=-1) {
    
        this->fFitHist = fithist;
        this->fBasIdx = basIdx; 
        this->fDrawRangeMin = drawRangeMin;
        this->fDrawRangeMax = drawRangeMax;
        this->fGlobNorm = globNorm; 
    }

    void SetParHist(TH1 *parHist) {
        this->fParHist = parHist;
    }

    void SetTotalFitFunc(TF1 *totalFitFunc) {
        this->fFit = totalFitFunc;
    }

    void SetBasIdx(int basIdx) {
        this->fBasIdx = basIdx;
    }

    void SetGlobNorm(bool setGlobNorm) {
        this->fGlobNorm = setGlobNorm;
    }

    void AddFitCompName(TString fitFuncComp) {
        cout << "Adding " << fitFuncComp << endl;
        this->fFitFuncComps.push_back(fitFuncComp);
    }

    void AddSplineHisto(TH1 *splineHisto) {
        this->fSplineHistos.push_back(splineHisto);
    }

    void EvaluateToBeDrawnComponents(std::vector<bool> onBaseline, std::vector<bool> multNorm, std::vector<bool> multGlobNorm, 
                            std::vector<double> shifts, int basIdx=-1, std::vector<TString> addComps = {""}, bool debug=1) {

        // Warnings that prevent the evaluation from being successful
        if(basIdx == -1){
            std::cerr << "Warning: Baseline is not fixed!" << std::endl;
        }
        if(onBaseline.size() != this->fFitFuncComps.size()){
            std::cerr << "Warning: onbaseline status not defined for all components!" << std::endl;
        }
        if(multNorm.size() != this->fFitFuncComps.size()){
            std::cerr << "Warning: multnorm status not defined for all components!" << std::endl;
        }
        if(multGlobNorm.size() != this->fFitFuncComps.size()){
            std::cerr << "Warning: multglobnorm status not defined for all components!" << std::endl;
        }

        cout << "Number of parameters in the histogram: " << this->fParHist->GetNbinsX() << endl;
        cout << "Number of components to be drawn: " << this->fFitFuncComps.size() << endl;
        for(int iFunc=0; iFunc<this->fFitFuncComps.size(); iFunc++) {
            cout << "Component " << this->fFitFuncComps[iFunc] << endl; 
        }
        // Evaluate the single fit components alone
        int startPar=0;
        int iSpline=0;
        std::vector<TF1 *> rawComps;
        std::vector<TSpline3 *> rawSplineComps;
        std::vector<int> nParsComps;
        for(int iFunc=0; iFunc<this->fFitFuncComps.size(); iFunc++) {
            if(this->fFitFuncComps[iFunc].Contains("spline")) {
                nParsComps.push_back(1);
                rawSplineComps.push_back(new TSpline3(this->fSplineHistos[iSpline]));
                rawComps.push_back(new TF1(this->fFitFuncComps[iFunc], 
                            [&, this, rawSplineComps]
                            (double *x, double *pars) {
                            return rawSplineComps.back()->Eval(x[0] - pars[0]);}, 
                            this->fDrawRangeMin, this->fDrawRangeMax, 1));
                cout << "Setting spline shift of " << this->fFitFuncComps[iFunc] << " to " << 
                        this->fParHist->GetXaxis()->GetLabels()->At(startPar+iFunc+1)->GetName() << endl; 
                cout << "Bin content: " << this->fParHist->GetBinContent(startPar+iFunc+2) << endl;
                rawComps.back()->FixParameter(0, this->fParHist->GetBinContent(startPar+iFunc+2));
                startPar += 1; // 1 + 1
                iSpline++;
            } else {
                nParsComps.push_back(std::get<1>(functions[this->fFitFuncComps[iFunc]]));
                rawComps.push_back(new TF1(this->fFitFuncComps[iFunc], std::get<0>(functions[this->fFitFuncComps[iFunc]]), 
                                           fDrawRangeMin, fDrawRangeMax, std::get<1>(functions[this->fFitFuncComps[iFunc]])));
                for(int iPar=0; iPar<rawComps.back()->GetNpar(); iPar++) {
                    // cout << hist->GetXaxis()->GetLabels()->At(bin-1)->GetName() << endl;
                    cout << "Setting par n. " << iPar << " of " << this->fFitFuncComps[iFunc] << " to " << 
                            this->fParHist->GetXaxis()->GetLabels()->At(startPar+iFunc+iPar+1+this->fGlobNorm)->GetName() << endl; 
                    // cout << "Setting par n. " << iPar << "of " << this->fFitFuncComps[iFunc] << " to " << this->fParHist->GetBinContent(startPar+iFunc+iPar) << endl; 
                    cout << "StartPar: " << startPar << endl;
                    cout << "iFunc: " << iFunc << endl;
                    cout << "iPar: " << iPar << endl;
                    cout << "Bin number: " << startPar+iFunc+iPar+1 << endl;
                    cout << "Bin content: " << this->fParHist->GetBinContent(startPar+iFunc+iPar+2+this->fGlobNorm) << endl;
                    cout << endl;
                    rawComps.back()->FixParameter(iPar, this->fParHist->GetBinContent(startPar+iFunc+iPar+2+this->fGlobNorm));
                }
                cout << "startPar: " << startPar << endl;
                startPar += std::get<1>(functions[this->fFitFuncComps[iFunc]]);
            }
            cout << endl;
        }

        // save the normalization constant for which each component has to be multiplied when drawing
        std::vector<double> norms;
        for(int iFunc=0; iFunc<this->fFitFuncComps.size(); iFunc++) {
            if(multNorm[iFunc]) {
                int normIdx = accumulate(nParsComps.begin(), std::next(nParsComps.begin(), iFunc), 0) + iFunc + this->fGlobNorm;
                cout << "Norm index " << normIdx << endl;
                cout << "Setting component norm to: " << this->fParHist->GetXaxis()->GetLabels()->At(normIdx)->GetName();
                cout << " , value: " << this->fParHist->GetBinContent(normIdx+1) << endl;
                norms.push_back(this->fParHist->GetBinContent(normIdx+1));
            } else {
                cout << "Setting component norm to 1" << endl;
                norms.push_back(1.);
            }
        } 

        // append to the raw components vector the functions that are sum of more than one component
        if(addComps[0] != "") {
            for(int iAddComp=0; iAddComp<addComps.size(); iAddComp++) {

                // push back the components multiplied by their norm
                rawComps.push_back(new TF1("SumComp_" + addComps[iAddComp], 
                        [&, this, rawComps, addComps]
                        (double *x, double *pars) {
                        double sum=0.;
                        for(int iFunc=0; iFunc<onBaseline.size(); iFunc++) {
                            if(addComps[iAddComp].Contains(std::to_string(iFunc))) {
                                sum += norms[iFunc] * rawComps[iFunc]->Eval(x[0]);
                            }
                        }
                        return sum;}, this->fDrawRangeMin, this->fDrawRangeMax, 0));

                // determine whether, when drawing, the newly added component has to be drawn on 
                // the baseline and has to be multiplied for the global normalization constant
                bool addBaseline = true;
                bool addMultGlobNorm = true;
                for(int iFunc=0; iFunc<onBaseline.size(); iFunc++) {
                    if(addComps[iAddComp].Contains(std::to_string(iFunc))) {

                        // all the components should have the property set to true for the sum component
                        // to also have the same feature
                        if(!onBaseline[iFunc]) {
                            cout << "Not all components are to be drawn on the baseline, their sum will not be ";
                            cout << "drawn on the baseline!" << endl;
                            addBaseline = false;
                        }
                        if(!multGlobNorm[iFunc]) {
                            cout << "Not all components are to be multiplied for the global normalization constant, their sum will not be ";
                            cout << "multiplied!" << endl;
                            addMultGlobNorm = false;
                        }

                    }
                }
                if(addBaseline) {
                    onBaseline.push_back(1);
                } else {
                    onBaseline.push_back(0);
                }
                if(addMultGlobNorm) {
                    multGlobNorm.push_back(1);
                } else {
                    multGlobNorm.push_back(0);
                }
     
                // kill the components to be summed by setting the normalization constants to zero 
                for(int iFunc=0; iFunc<this->fFitFuncComps.size(); iFunc++) {
                    if(addComps[iAddComp].Contains(std::to_string(iFunc))) {
                        onBaseline[iFunc] = false;
                        multNorm[iFunc] = 0.0000;
                        multGlobNorm[iFunc] = 0.0000;
                    }
                }
            }
        }

        // Define the baseline with its norm, if not indicated it is set to 1
        double basNorm;
        TF1 *bas = nullptr;
        if(basIdx != -1) {
            int previousCompsPars = accumulate(nParsComps.begin(), std::next(nParsComps.begin(), basIdx), 0) + basIdx;
            cout << "Setting baseline norm to: " << this->fParHist->GetBinContent(previousCompsPars + this->fGlobNorm + 1) << endl;
            basNorm = this->fParHist->GetBinContent(previousCompsPars + this->fGlobNorm + 1);
            bas = new TF1(this->fFitFuncComps[basIdx],
                [&, this, rawComps, multNorm, multGlobNorm, basIdx]
                    (double *x, double *pars) {
                       return rawComps[basIdx]->Eval(x[0]);
                    }, this->fDrawRangeMin, this->fDrawRangeMax, 0);
        } else {
            basNorm = 1.000000;
            bas = new TF1(this->fFitFuncComps[basIdx],
                [&, this]
                    (double *x, double *pars) {
                       return 1.0000;
                    }, this->fDrawRangeMin, this->fDrawRangeMax, 0);
        }

        // Define the weight of the baseline for each component, if it specified not to draw the component
        // on the baseline it will be set to zero 
        std::vector<double> onBasNorms;
        for(int iFunc=0; iFunc<this->fFitFuncComps.size(); iFunc++) {
            if(onBaseline[iFunc]) {
                cout << "Setting norm of the baseline to: " << basNorm << endl;
                onBasNorms.push_back(basNorm);
            } else {
                cout << "Setting norm of the baseline to 0" << endl;
                onBasNorms.push_back(0.00000);
            }
        } 

        // Define the global normalization constant for which every component will be multiplied, 1 if we want 
        // to draw the component as not multiplied
        std::vector<double> globNorms;
        cout << "Number of global norms " << multGlobNorm.size()  << endl;
        for(int iFunc=0; iFunc<this->fFitFuncComps.size(); iFunc++) {
            cout << "Component number " << iFunc << endl;
            if(multGlobNorm[iFunc]) {
                cout << "Setting global norm to: " << this->fGlobNorm << endl;
                globNorms.push_back(this->fGlobNorm);
            } else {
                cout << "Setting global norm to 1" << endl;
                globNorms.push_back(1.);
            }
        } 

        // Define the final functions that will be drawn on the canvas
        for(int iRawComp=0; iRawComp<rawComps.size(); iRawComp++) {
            this->fDrawFuncs.push_back(new TF1(this->fFitFuncComps[iRawComp],
                [&, this, globNorms, iRawComp, norms, rawComps, onBasNorms, bas]
                (double *x, double *pars) {
                   return globNorms[iRawComp] * (norms[iRawComp]*rawComps[iRawComp]->Eval(x[0]) + onBasNorms[iRawComp]*bas->Eval(x[0]));  
                }, this->fDrawRangeMin, this->fDrawRangeMax, 0));
        }
        cout << "Raw components defined!" << endl;

        // Debug couts
        if(debug) {
            cout << "Norms: " << endl;
            cout << "[";
            for(int iNorm=0; iNorm<norms.size(); iNorm++) {
                if(iNorm<norms.size()-1) cout << norms[iNorm] << ", ";
                else cout << norms[iNorm];
            }
            cout << "]" << endl;
            cout << endl;
            cout << "onBasNorms: " << endl;
            cout << "[";
            for(int iOnBasNorm=0; iOnBasNorm<onBasNorms.size(); iOnBasNorm++) {
                if(iOnBasNorm<onBasNorms.size()-1) cout << onBasNorms[iOnBasNorm] << ", ";
                else cout << onBasNorms[iOnBasNorm];
            }
            cout << "]" << endl;
            cout << endl;
            cout << "globNorms: " << endl;
            cout << "[";
            for(int iGlobNorm=0; iGlobNorm<globNorms.size(); iGlobNorm++) {
                if(iGlobNorm<globNorms.size()-1) cout << globNorms[iGlobNorm] << ", ";
                else cout << globNorms[iGlobNorm];
            }
            cout << "]" << endl;
            cout << endl;
        }
    }

    /*
    Define a canvas before calling this function and pass gPad as TVirtualPad
    */
    void Draw(TVirtualPad *pad, std::vector<TString> legLabels, std::vector<double> legCoords, int linesThickness, 
              double lowRangeUser=0.0, double uppRangeUser=1.05, std::string title=";k* (MeV/c);C(k*)") {

        cout << "Start drawing!" << endl;
        pad->cd();
        double yMinDraw = lowRangeUser;
        double yMaxDraw = uppRangeUser + fFitHist->GetMaximum();
        
        TLegend *legend = new TLegend(legCoords[0], legCoords[1], legCoords[2], legCoords[3]);
        legend->AddEntry(this->fFitHist, legLabels[0].Data(), "lp");
        legend->AddEntry(this->fFit, legLabels[1].Data(), "l");

        gPad->DrawFrame(fDrawRangeMin, yMinDraw, fDrawRangeMax, yMaxDraw, title.data());
                
        std::vector<Color_t> colors = {kMagenta + 3, kAzure + 2, kGreen, kOrange, kBlue + 2, kCyan, kBlack, kGreen+2};
        cout << "Number of components to be drawn: " << fDrawFuncs.size() << endl;
        for(int iFuncEval=0; iFuncEval<fDrawFuncs.size(); iFuncEval++) {
            cout << "Drawing the component " << iFuncEval << endl;
            this->fDrawFuncs[iFuncEval]->SetNpx(300);
            this->fDrawFuncs[iFuncEval]->SetLineColor(colors[iFuncEval]); //.data());
            this->fDrawFuncs[iFuncEval]->SetLineWidth(linesThickness);
            this->fDrawFuncs[iFuncEval]->DrawF1(fDrawRangeMin+1,fDrawRangeMax,"same");
            this->fDrawFuncs[iFuncEval]->Draw("same");
            cout << "Legend label: " << legLabels[iFuncEval+2] << endl;
            if(legLabels[iFuncEval+2].Contains("lambda_flat")) continue;
            legend->AddEntry(this->fDrawFuncs[iFuncEval], legLabels[iFuncEval+2].Data(), "l");
        }

        this->fFit->SetNpx(300);
        this->fFit->SetLineColor(kRed);
        this->fFit->SetLineWidth(linesThickness);
        this->fFit->DrawF1(fDrawRangeMin+1,fDrawRangeMax,"same");
        pad->Update();

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

        cout << "Finish drawing!" << endl;
    }

   private:

    TH1 *fFitHist = nullptr;
    TH1 *fParHist = nullptr;
    TF1 *fFit = nullptr;
    bool fGlobNorm; 
    int fBasIdx; 
    double fDrawRangeMin;
    double fDrawRangeMax;

    std::vector<TString> fFitFuncComps;                             // Function names of fit components
    std::vector<TF1*> fDrawFuncs;                                   // Fit components evaluated after the fitting
    std::vector<TF1*> fRawFuncs;                                   // Fit components evaluated after the fitting
    std::vector<TH1*> fSplineHistos;                                // Fit components evaluated after the fitting

};

#endif  // FEMPY_DRAWFITFUNCTS_HXX_
