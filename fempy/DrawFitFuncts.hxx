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

#if LOG_LEVEL
#define DEBUG(msg) std::cout << msg << std::endl
#else
#define DEBUG(msg)
#endif

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

        DEBUG("Total number of parameters with subcomponents: " << this->fParHist->GetNbinsX());
        DEBUG("Number of components to be drawn: " << this->fFitFuncComps.size());
        for(int iFunc=0; iFunc<this->fFitFuncComps.size(); iFunc++) {
            DEBUG("Component " << this->fFitFuncComps[iFunc]); 
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
                DEBUG("Set spline shift to " << 
                      this->fParHist->GetBinContent(startPar+iFunc+2); 
                      cout << ", content of bin no. " << startPar+iFunc+2;
                      cout << ", label " << this->fParHist->GetXaxis()->GetLabels()->At(startPar+iFunc+1)->GetName());
                rawComps.back()->FixParameter(0, this->fParHist->GetBinContent(startPar+iFunc+2));
                startPar += 1; 
                iSpline++;
            } else {
                nParsComps.push_back(std::get<1>(functions[this->fFitFuncComps[iFunc]]));
                rawComps.push_back(new TF1(this->fFitFuncComps[iFunc], std::get<0>(functions[this->fFitFuncComps[iFunc]]), 
                                           fDrawRangeMin, fDrawRangeMax, std::get<1>(functions[this->fFitFuncComps[iFunc]])));
                DEBUG("--------------------------------");
                DEBUG("Set pars of comp " << iFunc << ", named " << this->fFitFuncComps[iFunc] << ", having " << rawComps.back()->GetNpar() << " parameters" << endl; 
                      cout << "StartPar: " << startPar);
                for(int iPar=0; iPar<rawComps.back()->GetNpar(); iPar++) {
                    DEBUG("Set par n. " << iPar << " to " << this->fParHist->GetXaxis()->GetLabels()->At(startPar+iFunc+iPar+1+this->fGlobNorm)->GetName(); 
                          cout << ", bin no. " << startPar+iFunc+iPar+1 << " bin content: " << this->fParHist->GetBinContent(startPar+iFunc+iPar+2+this->fGlobNorm));
                    rawComps.back()->FixParameter(iPar, this->fParHist->GetBinContent(startPar+iFunc+iPar+2+this->fGlobNorm));
                }
                startPar += std::get<1>(functions[this->fFitFuncComps[iFunc]]);
                DEBUG("--------------------------------");
            }
        }

        // save the normalization constant for which each component has to be multiplied when drawing
        std::vector<double> norms;
        DEBUG("--------------------------------");
        std::cout << std::showpos;
        cout.precision(4);
        std::cout << std::scientific;
        for(int iFunc=0; iFunc<this->fFitFuncComps.size(); iFunc++) {
            if(multNorm[iFunc]) {
                int normIdx = accumulate(nParsComps.begin(), std::next(nParsComps.begin(), iFunc), 0) + iFunc + this->fGlobNorm;
                DEBUG("Set component " + std::to_string(iFunc) + " norm to: ";
                cout << this->fParHist->GetXaxis()->GetLabels()->At(normIdx)->GetName() << ", val: ";
                cout << this->fParHist->GetBinContent(normIdx+1) << ", norm idx: ";
                cout << std::to_string(normIdx));
                norms.push_back(this->fParHist->GetBinContent(normIdx+1));
            } else {
                DEBUG("Set component norm to 1");
                norms.push_back(1.);
            }
        } 
        DEBUG("--------------------------------");
        std::cout << std::noshowpos;

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
        DEBUG("--------------------------------");
        if(basIdx != -1) {
            int previousCompsPars = accumulate(nParsComps.begin(), std::next(nParsComps.begin(), basIdx), 0) + basIdx;
            DEBUG("Set baseline norm for function to: " << this->fParHist->GetBinContent(previousCompsPars + this->fGlobNorm + 1));
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
        DEBUG("--------------------------------");

        // Define the weight of the baseline for each component, if it specified not to draw the component
        // on the baseline it will be set to zero 
        std::vector<double> onBasNorms;
        DEBUG("--------------------------------");
        for(int iFunc=0; iFunc<this->fFitFuncComps.size(); iFunc++) {
            if(onBaseline[iFunc]) {
                DEBUG("Set norm of the baseline for function " << iFunc << " to: " << basNorm);
                onBasNorms.push_back(basNorm);
            } else {
                DEBUG("Set norm of the baseline for function " << iFunc << " to: " << static_cast<double>(int(0)));
                onBasNorms.push_back(0.00000);
            }
        } 
        DEBUG("--------------------------------");

        // Define the global normalization constant for which every component will be multiplied, 1 if we want 
        // to draw the component as not multiplied
        std::vector<double> globNorms;
        DEBUG("--------------------------------");
        DEBUG("Number of global norms " << multGlobNorm.size());
        for(int iFunc=0; iFunc<this->fFitFuncComps.size(); iFunc++) {
            if(multGlobNorm[iFunc]) {
                DEBUG("Set global norm for function " << iFunc << " to: " << this->fGlobNorm);
                globNorms.push_back(this->fGlobNorm);
            } else {
                DEBUG("Set global norm for function " << iFunc << " to 1");
                globNorms.push_back(1.);
            }
        } 
        DEBUG("--------------------------------");

        // Define the final functions that will be drawn on the canvas
        for(int iRawComp=0; iRawComp<rawComps.size(); iRawComp++) {
            this->fDrawFuncs.push_back(new TF1(this->fFitFuncComps[iRawComp],
                [&, this, globNorms, iRawComp, norms, rawComps, onBasNorms, bas]
                (double *x, double *pars) {
                   return globNorms[iRawComp] * (norms[iRawComp]*rawComps[iRawComp]->Eval(x[0]) + onBasNorms[iRawComp]*bas->Eval(x[0]));  
                }, this->fDrawRangeMin, this->fDrawRangeMax, 0));
        }
        DEBUG("Raw components defined!");

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
        DEBUG("--------------------------------");
        DEBUG("Number of components to be drawn: " << fDrawFuncs.size());
        for(int iFuncEval=0; iFuncEval<fDrawFuncs.size(); iFuncEval++) {
            this->fDrawFuncs[iFuncEval]->SetNpx(300);
            this->fDrawFuncs[iFuncEval]->SetLineColor(colors[iFuncEval]); //.data());
            this->fDrawFuncs[iFuncEval]->SetLineWidth(linesThickness);
            this->fDrawFuncs[iFuncEval]->DrawF1(fDrawRangeMin+1,fDrawRangeMax,"same");
            this->fDrawFuncs[iFuncEval]->Draw("same");
            DEBUG("Drawing the component " << iFuncEval << " with legend label: " << legLabels[iFuncEval+2]);
            if(legLabels[iFuncEval+2].Contains("lambda_flat")) continue;
            legend->AddEntry(this->fDrawFuncs[iFuncEval], legLabels[iFuncEval+2].Data(), "l");
        }
        DEBUG("--------------------------------");

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
    std::vector<TF1*> fRawFuncs;                                    // Fit components evaluated after the fitting
    std::vector<TH1*> fSplineHistos;                                // Fit components evaluated after the fitting

};

#endif  // FEMPY_DRAWFITFUNCTS_HXX_
