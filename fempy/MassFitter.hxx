#ifndef FEMPY_MASSFITTER_HXX_
#define FEMPY_MASSFITTER_HXX_

#include <map>
#include <tuple>
#include <string>
#include "yaml-cpp/yaml.h"


double Gaus(double *x, double *par) {
    double norm = 1. / TMath::Sqrt((2. * TMath::Pi())) / par[2];
    return norm * par[0] * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / 2. / par[2] / par[2]);
}

double DoubleGaus(double *x, double *par) {
    double norm1 = 1. / TMath::Sqrt((2. * TMath::Pi())) / par[2];
    double norm2 = 1. / TMath::Sqrt((2. * TMath::Pi())) / par[5];
    return norm1 * par[0] * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / 2. / par[2] / par[2]) +
           norm2 * par[3] * TMath::Exp(-(x[0] - par[4]) * (x[0] - par[4]) / 2. / par[5] / par[5]);
}

double Hat(double *x, double *par) {
    // p0: total yield
    // p1: mean
    // p2: sigma of thin gaussian
    // p3: fraction of narrow gaussian yield
    // p4: wide/narrow gaussian width

    double normThin = 1. / TMath::Sqrt((2. * TMath::Pi())) / par[2];
    double gThin = normThin * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / 2. / par[2] / par[2]);

    double wideSigma = par[2] * par[4];
    double normWide = 1. / TMath::Sqrt((2. * TMath::Pi())) / wideSigma;
    double gWide = normWide * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / 2. / wideSigma / wideSigma);

    return par[0] * (par[3] * gThin + (1 - par[3]) * gWide);
}

double PowEx(double *x, double *par) {
    // p0: total yield
    // p1: slope
    Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    if (x[0] < mpi) return 0;
    return par[0] * TMath::Sqrt(x[0] - mpi) * TMath::Exp(-1. * par[1] * (x[0] - mpi));
}

double Pol1(double *x, double *par) { return par[0] + par[1] * x[0]; }

double Pol5(double *x, double *par) { return par[0] + 
                                             par[1] * x[0] + 
                                             par[2] * x[0] * x[0] +
                                             par[3] * x[0] * x[0] * x[0] +
                                             par[4] * x[0] * x[0] * x[0] * x[0] + 
                                             par[5] * x[0] * x[0] * x[0] * x[0] * x[0]; }

double Exp(double *x, double *par) {
    // p0: total yield
    // p1: slope
    return par[0] * TMath::Exp(par[1] * x[0]);
}

double Spline3(double *x, double *par){
    int numKnots = 6;
    Double_t xKnots[numKnots];
    Double_t yKnots[numKnots];
    for(int iKnot=0; iKnot<numKnots; iKnot++){
        xKnots[iKnot] = par[iKnot];
        yKnots[iKnot] = par[numKnots+iKnot];
    }
    TSpline3* sp3 = new TSpline3("sp3", xKnots, yKnots, numKnots, "");
    return sp3->Eval(x[0]);
}

double Spline5(double *x, double *par){
    int numKnots = 6;
    Double_t xKnots[numKnots];
    Double_t yKnots[numKnots];
    for(int iKnot=0; iKnot<numKnots; iKnot++){
        xKnots[iKnot] = par[iKnot];
        yKnots[iKnot] = par[numKnots+iKnot];
    }
    TSpline5* sp5 = new TSpline5("sp5", xKnots, yKnots, numKnots, "");
    return sp5->Eval(x[0]);
}

class MassFitter {
 public:
    std::map<std::string, int> nPars = {{"ngaus", 3}, {"gaus", 3}, {"doublegaus", 6}, {"hat", 5},
                                        {"pol1", 2}, {"pol5", 7}, {"powex", 2}, {"exp", 2}, 
                                        {"spline3", 12}, {"spline5", 12}}; 

    enum SgnFuncs { kGaus = 0 }; // not used, delete?
    enum BkgFuncs { kPol1 = 0 }; // not used, delete?
    MassFitter(TH1 *hist, std::string sgnFuncName, std::string bkgFuncName, double fitRangeMin, double fitRangeMax) {
        this->fHist = reinterpret_cast<TH1 *>(hist->Clone());

        this->fFitRangeMin = fitRangeMin;
        this->fFitRangeMax = fitRangeMax;

        this->fIntLowEdge = fitRangeMin;
        this->fIntUppEdge = fitRangeMax;
        
        this->fSgnFuncName = sgnFuncName;
        this->fBkgFuncName = bkgFuncName;        

        this->fNSgnPars = nPars[this->fSgnFuncName];
        this->fNBkgPars = nPars[this->fBkgFuncName];
        int nTotPars = this->fNSgnPars + this->fNBkgPars; 

        if (fSgnFuncName == "gaus") {
            this->fSgnFunc = Gaus;
        } else if (fSgnFuncName == "hat") {
            this->fSgnFunc = Hat;
        } else if (fSgnFuncName == "doublegaus") {
            this->fSgnFunc = DoubleGaus;
        } else {
            printf("Function not implemented\n");
            exit(1);
        }

        if (fBkgFuncName == "powex") {
            this->fBkgFunc = PowEx;
        } else if (fBkgFuncName == "exp") {
            this->fBkgFunc = Exp;
        } else if (fBkgFuncName == "pol1") {
            this->fBkgFunc = Pol1;
        } else if (fBkgFuncName == "pol5") {
            this->fBkgFunc = Pol5;
        } else if (fBkgFuncName == "spline3") {
            this->fBkgFunc = Spline3;    
        } else if (fBkgFuncName == "spline5") {
            this->fBkgFunc = Spline5;    
        } else {
            printf("Function not implemented\n");
            exit(1);
        }

        this->fFit = new TF1(
            "fTot",
            [&, this](double *x, double *pars) {
                return this->fSgnFunc(x, pars) + this->fBkgFunc(x, &pars[this->fNSgnPars]);
            },
            fFitRangeMin, fFitRangeMax, nTotPars);
        fFit->SetNpx(300);

        this->fPrefit = new TF1(
            "fPrefit",
            [&, this, bkgFuncName](double *x, double *pars) -> double {
                if (fBkgFuncName == "powex" && std::abs(x[0] - 0.1455) < 0.003) {
                    TF1::RejectPoint();
                    return this->fBkgFunc(x, pars);
                } else if (fBkgFuncName == "exp" && std::abs(x[0] - 1.87) < 0.05) {
                    TF1::RejectPoint();
                    return this->fBkgFunc(x, pars);
                } else if (fBkgFuncName == "spline3" && std::abs(x[0] - 1.115683) < 0.012) {
                    fPrefitRangeMin = 1.115683 - 0.012;
                    fPrefitRangeMax = 1.115683 + 0.012;
                    TF1::RejectPoint();
                    return this->fBkgFunc(x, pars);
                } else if (fBkgFuncName == "spline5" && std::abs(x[0] - 1.115683) < 0.012) {
                    fPrefitRangeMin = 1.115683 - 0.012;
                    fPrefitRangeMax = 1.115683 + 0.012;
                    TF1::RejectPoint();
                    return this->fBkgFunc(x, pars);
                } else {
                    return this->fBkgFunc(x, pars);
                }
            },
            this->fFitRangeMin, this->fFitRangeMax, this->fNBkgPars);
    }

    void SetIntegrationEdges(double lowEdge, double uppEdge) {
        cout << this->fIntLowEdge << endl; 
        cout << this->fIntUppEdge << endl; 
        this->fIntLowEdge = lowEdge;
        this->fIntUppEdge = uppEdge;
        cout << this->fIntLowEdge << endl; 
        cout << this->fIntUppEdge << endl;
    }

    void SetIntegrationEdges(double nSigma) {
        if (!fFit) return;
        double probs[2] = {TMath::Freq(-nSigma), TMath::Freq(nSigma)};
        double quantiles[2];
        fSgn->GetQuantiles(2, quantiles, probs);
        this->fIntLowEdge = quantiles[0];
        this->fIntUppEdge = quantiles[1];
    }

    /*
    It is your responsibility to give the correct number of parameters. All fit parameters must be initialized.
    Set the fit parameters. Signal parameters first, then background.
    format: {index_of_parameter, {init_value, lower_lim, upper_lim}}
    */
    void SetFitSettings(std::string name) {
        if (name == "421_gaus_exp_alice3_siblings") {
            fFitPars = {
                // gaussian
                {0, {"norm", 1000, 0, 1.e6}},
                {1, {"mean", 1.86484, 1.86484*0.98, 1.86484*1.02}},
                {2, {"sigma", 10, 0.005, 0.05}},

                // exp
                {3, {"norm", 10, 0.01, 1.e4}},
                {4, {"slope", -0.5, -10, 10}},
            };
        } else if (name == "413_gaus_powex_DstarFemto_SE") {
            fFitPars = {
                // gaussian
                {0, {"norm", 1, 0.05, 5}},
                {1, {"mean", 0.1455, 0.1443, 0.1458}},
                {2, {"sigma", 0.0006, 0.0005, 0.0007}},
                
                // powex
                {3, {"norm", 2000, 100, 6000}},
                {4, {"slope", 20, 10, 30}},
            };
        } else if (name == "413_gaus_powex_DstarKFemto_SE") {
            fFitPars = {
                // gaussian
                {0, {"norm", 1, 0.01, 5}},
                {1, {"mean", 0.1455, 0.1443, 0.1458}},
                {2, {"sigma", 0.0006, 0.0005, 0.0007}},
                
                // powex
                {3, {"norm", 2000, 100, 6000}},
                {4, {"slope", 20, 2, 30}},
            };
        } else if (name == "413_gaus_powex_DstarFemto_ME") {
            fFitPars = {
                // gaussian
                {0, {"norm", 200, 10, 500}},
                {1, {"mean", 0.1455, 0.1443, 0.1458}},
                {2, {"sigma", 0.0006, 0.0005, 0.0007}},
                
                // powex
                {3, {"norm", 200000, 10000, 500000}},
                {4, {"slope", 20, 15, 25}},
            };
        } else if (name == "411_gaus_exp_DstarFemto") {
            fFitPars = { {0, {"norm", 0.5, 0, 50}},
                {0, {"mean", 1.87, 1.84, 1.9}},
                {2, {"sigma", 0.006, 0.0001, 0.01}},

                // exp
                {3, {"norm", 10, 0, 1.e6}},
                {4, {"slope", -0.5, -10, 0}},
            };
        } else if (name == "411_hat_exp_DstarFemto") {
            fFitPars = {    // gauss1
                {0, {"norm", 0.1, 0, 50}},
                {1, {"mean", 0.145, 0.144, 0.146}},
                {2, {"sigma", 0.001, 0.0002, 0.002}},

                // gauss2
                {3, {"Yfrac", 0.5, 0.1, 1}},
                {4, {"sigmaFrac", 1.5, 1, 5}},

                // exp
                {5, {"norm", 1, 0, 1e6}},
                {6, {"slope", 0.1, -100, 100}},
            };
        } else {
            std::cout << "The set of parameters '" << name << "' is not valid. Exit!" << std::endl;
            exit(1);
        }
    }

    void Add(std::vector<std::tuple<std::string, double, double, double>> pars) {
        for (int iPar = 0; iPar < this->fNSgnPars + this->fNBkgPars; iPar++){
            fFitPars.insert({iPar, pars[iPar]});
        }
    }

    void AddPrefit(std::vector<std::tuple<std::string, double, double, double>> pars) {
        for (int iPar = 0; iPar < this->fNBkgPars; iPar++){
            fBkgPreFitPars.insert({iPar, pars[iPar]});
        }
    }

    int Prefit() {
        for (int iPar = 0; iPar < this->fNBkgPars; iPar++) {
            auto pars = fBkgPreFitPars[iPar];
            //std::cout << iPar << " " 
            //          << std::get<0>(pars)
            //          << " " << std::get<1>(pars)
            //          << " " << std::get<2>(pars)
            //          << " " << std::get<3>(pars)
            //          << std::endl;
            this->fPrefit->SetParName(iPar, std::get<0>(pars).data());
            this->fPrefit->SetParameter(iPar, std::get<1>(pars));
            this->fPrefit->SetParLimits(iPar, std::get<2>(pars), std::get<3>(pars));
        }
        
        int statusPrefit = fHist->Fit(this->fPrefit, "SMRL+0", "")->Status();

        for(int iPar = this->fNSgnPars; iPar < this->fNSgnPars + this->fNBkgPars; iPar++) {
            double lowParLimit, uppParLimit;
            this->fPrefit->GetParLimits(iPar - this->fNSgnPars, lowParLimit, uppParLimit);
            fFitPars[iPar] = {this->fPrefit->GetParName(iPar - this->fNSgnPars), 
                              this->fPrefit->GetParameter(iPar - this->fNSgnPars),
                              lowParLimit, uppParLimit};
        }
        //cout << "Prefit done!" << endl;
        return statusPrefit;
    }

    int Fit() {
        
        for (int iPar = 0; iPar < this->fNSgnPars + this->fNBkgPars; iPar++) {
            auto pars = fFitPars[iPar];

            this->fFit->SetParName(iPar, std::get<0>(pars).data());
            this->fFit->SetParameter(iPar, std::get<1>(pars));
            this->fFit->SetParLimits(iPar, std::get<2>(pars), std::get<3>(pars));
        }

        printf("\n\n\nPerforming the full fit:\n");
        int status = fHist->Fit(this->fFit, "SMRL+0", "")->Status();

        // decompose the fit function in its contributions
        this->fBkg = new TF1(fBkgFuncName.data(), fBkgFunc, fFitRangeMin, fFitRangeMax, fNBkgPars);
        for(int iBkgPar = 0; iBkgPar < fNBkgPars; iBkgPar++){
            this->fBkg->SetParameter(iBkgPar, this->fFit->GetParameter(this->fNSgnPars + iBkgPar));
        }

        if (this->fSgnFuncName == "hat") {
            this->fHatThin = new TF1("fHatThin", Gaus, fFitRangeMin, fFitRangeMax, 3);
            this->fHatThin->SetParameter(0, this->fFit->GetParameter(0) * this->fFit->GetParameter(3));
            this->fHatThin->SetParameter(1, this->fFit->GetParameter(1));
            this->fHatThin->SetParameter(2, this->fFit->GetParameter(2));

            this->fHatWide = new TF1("fHatWide", Gaus, fFitRangeMin, fFitRangeMax, 3);
            this->fHatWide->SetParameter(0, this->fFit->GetParameter(0) * (1 - this->fFit->GetParameter(3)));
            this->fHatWide->SetParameter(1, this->fFit->GetParameter(1));
            this->fHatWide->SetParameter(2, this->fFit->GetParameter(2) * this->fFit->GetParameter(4));

            this->fSgn = new TF1("fSgn", Hat, fFitRangeMin, fFitRangeMax, 5);
            this->fSgn->SetParameter(0, this->fFit->GetParameter(0));
            this->fSgn->SetParameter(1, this->fFit->GetParameter(1));
            this->fSgn->SetParameter(2, this->fFit->GetParameter(2));
            this->fSgn->SetParameter(3, this->fFit->GetParameter(3));
            this->fSgn->SetParameter(4, this->fFit->GetParameter(4));

            this->fSgn = new TF1(fSgnFuncName.data(), fSgnFunc, fFitRangeMin, fFitRangeMax, fNSgnPars);
            for(int iSgnPar = 0; iSgnPar < fNSgnPars; iSgnPar++){
                this->fSgn->SetParameter(iSgnPar, this->fFit->GetParameter(iSgnPar));
            }
        }
        if (this->fSgnFuncName == "doublegaus") {
            this->fHatThin = new TF1("fHatThin", Gaus, fFitRangeMin, fFitRangeMax, 3);
            this->fHatThin->SetParameter(0, this->fFit->GetParameter(0));
            this->fHatThin->SetParameter(1, this->fFit->GetParameter(1));
            this->fHatThin->SetParameter(2, this->fFit->GetParameter(2));

            this->fHatWide = new TF1("fHatWide", Gaus, fFitRangeMin, fFitRangeMax, 3);
            this->fHatWide->SetParameter(0, this->fFit->GetParameter(3));
            this->fHatWide->SetParameter(1, this->fFit->GetParameter(4));
            this->fHatWide->SetParameter(2, this->fFit->GetParameter(5));

            this->fSgn = new TF1("fSgn", Hat, fFitRangeMin, fFitRangeMax, 5);
            this->fSgn->SetParameter(0, this->fFit->GetParameter(0));
            this->fSgn->SetParameter(1, this->fFit->GetParameter(1));
            this->fSgn->SetParameter(2, this->fFit->GetParameter(2));
            this->fSgn->SetParameter(3, this->fFit->GetParameter(3));
            this->fSgn->SetParameter(4, this->fFit->GetParameter(4));
            this->fSgn->SetParameter(5, this->fFit->GetParameter(5));

            this->fSgn = new TF1(fSgnFuncName.data(), fSgnFunc, fFitRangeMin, fFitRangeMax, fNSgnPars);
            for(int iSgnPar = 0; iSgnPar < fNSgnPars; iSgnPar++){
                this->fSgn->SetParameter(iSgnPar, this->fFit->GetParameter(iSgnPar));
            }
        }
        return status;
    }

    /*
    Define a canvas before calling this function and pass gPad as TVirtualPad
    */
    void Draw(TVirtualPad *pad, TString fitDrawOpts, std::string method) {
        pad->cd();
        fHist->GetYaxis()->SetRangeUser(0, 1.3 * fHist->GetMaximum());
        gPad->DrawFrame(fFitRangeMin, 0, fFitRangeMax, 1.3 * fHist->GetMaximum(),
                        ";M(p#pi) (GeV/c^{2});Counts");

        if(fitDrawOpts.Contains('i')) {
            TLatex tl;
            tl.SetTextSize(0.035);
            tl.SetTextFont(42);
            double nSigma = 2;
            double step = 0.05;
            int iStep = 0;

            tl.DrawLatexNDC(.15, .85 - step * iStep++, Form("#chi^{2}/NDF = %.2f", 
            fFit->GetChisquare() / fFit->GetNDF()));
            pad->Update();
            tl.DrawLatexNDC(.15, .85 - step * iStep++, Form("#chi^{2}/NDF SgnWindow = %.2f", 
                            GetSgnWindowChi2Ndf()));
            pad->Update();
            tl.DrawLatexNDC(.15, .85 - step * iStep++,
                            Form("S(%.2f#sigma) = %.2f #pm %.2f", nSigma, 
                                 this->GetSignal(method), this->GetSignalUnc(method)));
            tl.DrawLatexNDC(.15, .85 - step * iStep++,
                Form("B(%.2f#sigma) = %.2f #pm %.2f", nSigma, this->GetBackground(), 
                     this->GetBackgroundUnc()));
            tl.DrawLatexNDC(.15, .85 - step * iStep++, Form("Counts = %.2f", this->GetCounts()));
            pad->Update();

        }

        TLine* canvaLine = new TLine(0, 0, 1, 1);
        double lowMult = 0.;
        double uppMult = 0.;
        int padCoord = 1;
        if(fitDrawOpts.Contains('c')) {
            lowMult = 0.6;
            uppMult = 1.4;
            padCoord = 0;
        }
        if(fitDrawOpts.Contains('s')) {
            canvaLine->SetLineStyle(1);
            canvaLine->SetLineWidth(3);
            canvaLine->SetLineColor(6);

            double lowBinContent = fHist->GetBinContent(fHist->FindBin(fIntLowEdge));
            canvaLine->DrawLine(fIntLowEdge, pad->GetUymin() * padCoord + lowBinContent * lowMult,
                                fIntLowEdge, pad->GetUymax() * 0.5 * padCoord + lowBinContent * uppMult);
            double uppBinContent = fHist->GetBinContent(fHist->FindBin(fIntUppEdge));
            canvaLine->DrawLine(fIntUppEdge, pad->GetUymin() * padCoord + uppBinContent * lowMult,
                                fIntUppEdge, pad->GetUymax() * 0.5 * padCoord + uppBinContent * uppMult);
        }

        if(fitDrawOpts.Contains('p')) {
            canvaLine->SetLineStyle(1);
            canvaLine->SetLineWidth(3);
            canvaLine->SetLineColor(4);

            double lowBinContent = fHist->GetBinContent(fHist->FindBin(fPrefitRangeMin));
            canvaLine->DrawLine(fPrefitRangeMin, pad->GetUymin() * padCoord + lowBinContent * lowMult,
                                fPrefitRangeMin, pad->GetUymax() * 0.5 * padCoord + lowBinContent * uppMult);
            double uppBinContent = fHist->GetBinContent(fHist->FindBin(fPrefitRangeMax));
            cout << uppBinContent << endl;
            canvaLine->DrawLine(fPrefitRangeMax, pad->GetUymin() * padCoord + uppBinContent * lowMult,
                                fPrefitRangeMax, pad->GetUymax() * 0.5 * padCoord + uppBinContent * uppMult);
        }

        if(fitDrawOpts.Contains('k')) {
            canvaLine->SetLineStyle(8);
            canvaLine->SetLineWidth(3);
            canvaLine->SetLineColor(1);

            for(int xKnot = 0; xKnot < this->fNBkgPars/2; xKnot++) {
                double xKnotPos = this->fBkg->GetParameter(xKnot);
                double yKnotPos = fHist->GetBinContent(fHist->FindBin(xKnotPos));
                canvaLine->DrawLine(xKnotPos, pad->GetUymin() * padCoord + yKnotPos * lowMult,
                                    xKnotPos, pad->GetUymax() * 0.5 * padCoord + yKnotPos * uppMult);
            }
        }

        this->fBkg->SetNpx(300);
        this->fBkg->SetLineColor(kGray + 2);
        this->fBkg->Draw("same");

        this->fSgn->SetNpx(300);
        this->fSgn->SetLineColor(kBlue + 2);
        this->fSgn->Draw("same");

        if (this->fSgnFuncName == "hat" || this->fSgnFuncName == "doublegaus") {
            this->fHatThin->SetLineColor(kMagenta + 3);
            this->fHatThin->SetNpx(300);
            this->fHatThin->Draw("same");

            this->fHatWide->SetNpx(300);
            this->fHatWide->SetLineColor(kAzure + 2);
            this->fHatWide->Draw("same");
        } 

        fPrefit->SetLineStyle(9);
        fPrefit->SetLineColor(kGray + 2);
        fPrefit->Draw("same");

        fFit->SetLineColor(kRed);
        fFit->SetLineWidth(3);
        fFit->Draw("same");

        hist->SetMarkerSize(1);
        hist->SetMarkerStyle(20);
        hist->SetMarkerColor(kBlack);
        hist->SetLineColor(kBlack);
        hist->SetLineWidth(2);
        hist->Draw("same pe");

        TLatex tl;
        tl.SetTextSize(0.035);
        tl.SetTextFont(42);
        double nSigma = 2;
        double step = 0.05;
        int iStep = 0;
        tl.DrawLatexNDC(.15, .85 - step * iStep++, Form("#chi^{2}/NDF = %.2f", fFit->GetChisquare() / fFit->GetNDF()));
        tl.DrawLatexNDC(.15, .85 - step * iStep++,
                        Form("S(%.2f#sigma) = %.2f #pm %.2f", nSigma, this->GetSignal(nSigma, method),
                             this->GetSignalUnc(nSigma, method)));

        tl.DrawLatexNDC(
            .15, .85 - step * iStep++,
            Form("B(%.2f#sigma) = %.2f #pm %.2f", nSigma, this->GetBackground(nSigma), this->GetBackgroundUnc(nSigma)));

        tl.DrawLatexNDC(.15, .85 - step * iStep++, Form("Counts = %.2f", this->GetCounts()));

        // Write on the plot the fit parameters
        for (int iPar = 0; iPar < fFit->GetNpar(); iPar++) {
            double par = fFit->GetParameter(iPar);
            double parMin;
            double parMax;
            fFit->GetParLimits(iPar, parMin, parMax);
            double range = parMax - parMin;
            
            // Check if the fit pars are at limit
            if (par - parMin < 1.e-4 * range || parMax - par < 1.e-4 * range) {
                tl.SetTextColor(2);
                tl.DrawLatexNDC(.6, .85 - step * iPar, Form("%s*** = %.2e", fFit->GetParName(iPar), par));
                tl.SetTextColor(1);
            } else {
                tl.DrawLatexNDC(.6, .85 - step * iPar, Form("%s = %.2e", fFit->GetParName(iPar), par));
            }
        }
        if (fFit->GetChisquare() / fFit->GetNDF() > 150) {
            TLatex tlDanger;
            tlDanger.SetTextSize(0.07);
            tlDanger.SetTextFont(42);
            tlDanger.SetTextColor(2);
            tlDanger.DrawLatexNDC(.5, .4, "Danger: #chi^{2}/NDF > 150");
        }
        
        pad->Update();
    }

    double GetMean() {
        if (!fFit) return -1;
        if (this->fSgnFuncName == "gaus" || this->fSgnFuncName == "hat") return fFit->GetParameter(1);
        return -1;
    }

    double GetMeanUnc() {
        if (!fFit) return -1;
        if (this->fSgnFuncName == "gaus" || this->fSgnFuncName == "hat") return fFit->GetParError(1);
        return -1;
    }

    double GetWidth() {
        if (!fFit) return -1;

        double probs[2] = {TMath::Freq(-1), TMath::Freq(1)};
        double quantiles[2];
        fSgn->GetQuantiles(2, quantiles, probs);
        return (quantiles[1] - quantiles[0]) / 2;
    }

    double GetWidthUnc() {
        if (!fFit) return -1;

        if (this->fSgnFuncName == "gaus")
            return fFit->GetParError(2);
        if (this->fSgnFuncName == "doublegaus")
            return TMath::Sqrt(fFit->GetParError(2)*fFit->GetParError(2) + fFit->GetParError(5)*fFit->GetParError(5));
        else if (this->fSgnFuncName == "hat")
            return fFit->GetParError(2) / fFit->GetParameter(2) * GetWidth();  // assume rel error is the same
        return -1;
    }

    double GetCounts() {
        if (!fFit) return -1;
        int firstBin = fHist->GetXaxis()->FindBin(fFitRangeMin * 1.0001);
        int lastBin = fHist->GetXaxis()->FindBin(fFitRangeMax * 0.9999);

        double totCounts = fHist->Integral(firstBin, lastBin);
        double bkgCounts =
            fBkg->Integral(fHist->GetXaxis()->GetBinLowEdge(firstBin), fHist->GetXaxis()->GetBinLowEdge(lastBin + 1)) /
            this->fHist->GetBinWidth(1);
        return totCounts - bkgCounts;
    }

    double GetSignal(std::string method) {
        if (!fFit) return -1;

        if (method == "sgn_int") {
            return fSgn->Integral(fIntLowEdge, fIntUppEdge) / this->fHist->GetBinWidth(1);
        } else if (method == "data_minus_bkg") {
            int firstBin = fHist->GetXaxis()->FindBin(fIntLowEdge * 1.0001);
            int lastBin = fHist->GetXaxis()->FindBin(fIntUppEdge * 0.9999);

            double data = this->fHist->Integral(firstBin, lastBin);
            double bkg = GetBackground();
            return data - bkg;
        } else {
            printf("'%s' is an invalid signal extraction method. Exit!\n", method.data());
            exit(1);
        }
        return -1;
    }
    double GetChi2Ndf() { return fFit->GetChisquare() / fFit->GetNDF(); }

    double GetSgnWindowChi2Ndf() { 

        double lowSgnWindowEdge = this->fFit->GetParameter(1) - 2 * this->GetWidth();
        double uppSgnWindowEdge = this->fFit->GetParameter(1) + 2 * this->GetWidth();

        double chiSquare = 0.;
        int nBins = 0;
        for(int iBin = 0; iBin < this->fHist->GetNbinsX(); iBin++) {
            double iBinCenter = this->fHist->GetBinCenter(iBin); 
            if(iBinCenter >= fFitRangeMin && iBinCenter <= fFitRangeMax) nBins++;
            if(iBinCenter >= lowSgnWindowEdge && iBinCenter <= uppSgnWindowEdge) { 
                double iBinContent = this->fHist->GetBinContent(iBin); 
                chiSquare += (this->fFit->Eval(iBinCenter) - iBinContent) * (this->fFit->Eval(iBinCenter) - iBinContent) / iBinContent;
            }
        }
        return chiSquare / (nBins - this->fNSgnPars - this->fNBkgPars/2);
    }

    double GetBackground() {
        return fBkg->Integral(fIntLowEdge, fIntUppEdge) / this->fHist->GetBinWidth(1);
    }

    double GetBackgroundUnc() {
        Int_t leftBand = this->fHist->FindBin(this->GetMean() - 6 * this->GetWidth());
        Int_t rightBand = this->fHist->FindBin(this->GetMean() + 6 * this->GetWidth());

        int start = this->fHist->FindBin(this->fFitRangeMin * 1.0001);
        int end = this->fHist->FindBin(this->fFitRangeMax * 0.9999);
        double SidebandBkg =
            this->fHist->Integral(1, leftBand) + this->fHist->Integral(rightBand, this->fHist->GetNbinsX());

        double sum2 = 0;
        for (Int_t i = 1; i <= leftBand; i++) {
            sum2 += this->fHist->GetBinError(i) * this->fHist->GetBinError(i);
        }
        for (Int_t i = rightBand; i <= this->fHist->GetNbinsX(); i++) {
            sum2 += this->fHist->GetBinError(i) * this->fHist->GetBinError(i);
        }

        return TMath::Sqrt(sum2) / SidebandBkg * this->GetBackground();
    }

    double GetBackgroundPrefit() {
        return fPrefit->Integral(fIntLowEdge, fIntUppEdge) / this->fHist->GetBinWidth(1);
    }

    double GetBackgroundPrefitUnc() {
        Int_t leftBand = this->fHist->FindBin(this->GetMean() - 6 * this->GetWidth());
        Int_t rightBand = this->fHist->FindBin(this->GetMean() + 6 * this->GetWidth());

        int start = this->fHist->FindBin(this->fFitRangeMin * 1.0001);
        int end = this->fHist->FindBin(this->fFitRangeMax * 0.9999);
        double SidebandBkg =
            this->fHist->Integral(1, leftBand) + this->fHist->Integral(rightBand, this->fHist->GetNbinsX());

        double sum2 = 0;
        for (Int_t i = 1; i <= leftBand; i++) {
            sum2 += this->fHist->GetBinError(i) * this->fHist->GetBinError(i);
        }
        for (Int_t i = rightBand; i <= this->fHist->GetNbinsX(); i++) {
            sum2 += this->fHist->GetBinError(i) * this->fHist->GetBinError(i);
        }

        return TMath::Sqrt(sum2) / SidebandBkg * this->GetBackgroundPrefit();
    }

    double GetSignalUnc(std::string method) {
        if (!fFit) return -1;
        if (this->GetSignal(method) <= 0) return 0;

        return fFit->GetParError(0) / fFit->GetParameter(0) * this->GetSignal(method);
    }

    double GetPurity(std::string method) {
        return this->GetSignal(method) / (this->GetSignal(method) + this->GetBackground()); 
    }

    double GetPurityAllUnc(std::string method) {
        return TMath::Sqrt( this->GetPurityStatUnc(method) * this->GetPurityStatUnc(method) +
                            this->GetPuritySystUnc(method) * this->GetPuritySystUnc(method));
    }
    
    double GetPurityStatUnc(std::string method) {
        
        double sgnInt = this->GetSignal(method);
        double sgnIntErr = this->GetSignalUnc(method);
            
        double bkgInt = this->GetBackground();
        double bkgIntErr = this->GetBackgroundUnc();

        double dSgn = 1 / (sgnInt + bkgInt) - sgnInt / ( (sgnInt + bkgInt)*(sgnInt + bkgInt) );
        double dBkg = sgnInt / ( (sgnInt + bkgInt)*(sgnInt + bkgInt) );

        return TMath::Sqrt( (dSgn * sgnIntErr)*(dSgn * sgnIntErr) + (dBkg * bkgIntErr)*(dBkg * bkgIntErr) );
    }

    double GetPuritySystUnc(std::string method) {
        return TMath::Abs( this->GetPurity(method) - this->GetPurityPrefit(method) ) / 2;
    }

    double GetPurityPrefit(std::string method) {
        return this->GetSignal(method) / (this->GetSignal(method) + this->GetBackgroundPrefit()); 
    }

    double GetPurityPrefitUnc(std::string method) {
        
        double sgnInt = this->GetSignal(method);
        double sgnIntErr = this->GetSignalUnc(method);
            
        double bkgInt = this->GetBackgroundPrefit();
        double bkgIntErr = this->GetBackgroundPrefitUnc();

        double dSgn = 1 / (sgnInt + bkgInt) - sgnInt / ( (sgnInt + bkgInt)*(sgnInt + bkgInt) );
        double dBkg = sgnInt / ( (sgnInt + bkgInt)*(sgnInt + bkgInt) );

        return TMath::Sqrt( (dSgn * sgnIntErr)*(dSgn * sgnIntErr) + (dBkg * bkgIntErr)*(dBkg * bkgIntErr) );
    }

    TF1* GetSgnFunc() {
        //if (!fFit) return -1;
        return fSgn;
    }

    TF1* GetBkgFunc() {
        //if (!fFit) return -1;
        return fBkg;
    }


 private:

    TH1 *fHist = nullptr;
    TF1 *fSgn = nullptr;
    TF1 *fHatThin = nullptr;
    TF1 *fHatWide = nullptr;
    TF1 *fBkg = nullptr;
    TF1 *fFit = nullptr;
    TF1 *fPrefit = nullptr;

    std::string fSgnFuncName;
    std::string fBkgFuncName;

    double (*fSgnFunc)(double *x, double *par);
    double (*fBkgFunc)(double *x, double *par);

    int fNSgnPars;
    int fNBkgPars;
    double fFitRangeMin;
    double fFitRangeMax;
    double fPrefitRangeMin;
    double fPrefitRangeMax;
    double fIntLowEdge;
    double fIntUppEdge;
    std::map<int, std::tuple<std::string, double, double, double>> fFitPars;
    std::map<int, std::tuple<std::string, double, double, double>> fBkgPreFitPars;
};

#endif  // FEMPY_MASSFITTER_HXX_
