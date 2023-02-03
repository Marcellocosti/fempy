#include "TLatex.h"
#include "TObject.h"

double Gaus(double *x, double *par) {
    double norm = 1. / TMath::Sqrt((2. * TMath::Pi())) / par[2];
    return norm * par[0] * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / 2. / par[2] / par[2]);
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
double Pol1(double *x, double *par) { return par[0] + par[1] * x[0]; };

class MassFitter {
   public:
    std::map<std::string, int> nPars = {
        {"ngaus", 3}, {"gaus", 3}, {"hat", 5}, {"pol1", 2}, {"powex", 2},
    };

    enum SgnFuncs { kGaus = 0 };
    enum BkgFuncs { kPol1 = 0 };
    MassFitter(TH1 *hist, std::string sgnFuncName, std::string bkgFuncName, double fitRangeMin, double fitRangeMax) {
        this->hist = (TH1 *)hist->Clone();
        this->fitRangeMin = fitRangeMin;
        this->fitRangeMax = fitRangeMax;
        this->sgnFuncName = sgnFuncName;
        this->bkgFuncName = bkgFuncName;

        this->nSgnPars = nPars[this->sgnFuncName];
        this->nBkgPars = nPars[this->bkgFuncName];
        int nTotPars = this->nSgnPars + this->nBkgPars;

        if (sgnFuncName == "gaus")
            this->sgnFunc = Gaus;
        else if (sgnFuncName == "hat")
            this->sgnFunc = Hat;
        else {
            printf("Function not implemented\n");
            exit(1);
        }

        if (bkgFuncName == "powex")
            this->bkgFunc = PowEx;
        else if (bkgFuncName == "pol1")
            this->bkgFunc = Pol1;
        else {
            printf("Function not implemented\n");
            exit(1);
        }

        this->fFit = new TF1(
            "fTot",
            [&, this](double *x, double *pars) {
                return this->sgnFunc(x, pars) + this->bkgFunc(x, &pars[this->nSgnPars]);
            },
            fitRangeMin, fitRangeMax, nTotPars);
        fFit->SetNpx(300);

        this->fPrefit = new TF1(
            "fPrefit",
            [&, this](double *x, double *pars) -> double {
                if (std::abs(x[0] - 0.1455) < 0.003) {
                    TF1::RejectPoint();
                    return this->bkgFunc(x, pars);

                    // return 0;
                }
                return this->bkgFunc(x, pars);
            },
            this->fitRangeMin, fitRangeMax, nBkgPars);

        if (sgnFuncName == "gaus") {
            this->fFit->SetParName(0, "norm");
            this->fFit->SetParameter(0, 0.1);
            this->fFit->SetParLimits(0, 0, 50);
            this->fFit->SetParName(1, "mean");
            this->fFit->SetParameter(1, 0.145);
            this->fFit->SetParLimits(1, 0.144, 0.146);
            this->fFit->SetParName(2, "sigma");
            this->fFit->SetParameter(2, 0.001);
            this->fFit->SetParLimits(2, 0.0002, 0.002);
        } else if (sgnFuncName == "hat") {
            // g1
            this->fFit->SetParName(0, "norm");
            this->fFit->SetParameter(0, 0.1);
            this->fFit->SetParLimits(0, 0, 50);
            this->fFit->SetParName(1, "mean");
            this->fFit->SetParameter(1, 0.145);
            this->fFit->SetParLimits(1, 0.144, 0.146);
            this->fFit->SetParName(2, "sigma");
            this->fFit->SetParameter(2, 0.001);
            this->fFit->SetParLimits(2, 0.0002, 0.002);

            // g2
            this->fFit->SetParName(3, "Yfrac");
            this->fFit->SetParameter(3, 0.5);
            this->fFit->SetParLimits(3, 0.1, 1);
            this->fFit->SetParName(4, "sigmaFrac");
            this->fFit->SetParameter(4, 1.5);
            this->fFit->SetParLimits(4, 1, 5);
        }

        if (bkgFuncName == "powex") {
            this->fFit->SetParName(this->nSgnPars + 0, "norm");
            this->fFit->SetParameter(this->nSgnPars + 0, 0.5);
            this->fFit->SetParLimits(this->nSgnPars + 0, 200, 3000);
            this->fFit->SetParName(this->nSgnPars + 1, "slope");
            this->fFit->SetParameter(this->nSgnPars + 1, 0.1);
            this->fFit->SetParLimits(this->nSgnPars + 1, 0, 100);
        } else if (bkgFuncName == "pol1")
            bkgFunc = Pol1;
    }

    void Fit() {
        // prefit
        // this->fPrefit = new TF1("fPrefit", [&, this](double *x, double * par) {
        //     if (std::abs(x[0] - 0.145) < 0.001)
        //         TF1::RejectPoint();
        //         return 0;
        //     return this->bkgFunc.EvalPar(x, par);

        // }, this->fitRangeMin, fitRangeMax, nBkgPars);
        printf("\n\n\nPerfomring the prefit to the background:\n");
        hist->Fit(this->fPrefit, "MR0+", "");

        // set the bkg parameters based on prefit
        for (int iPar = 0; iPar < this->nBkgPars; iPar++) {
            fFit->SetParLimits(this->nSgnPars + iPar, fPrefit->GetParameter(iPar) * 0.8,
                               fPrefit->GetParameter(iPar) * 1.2);
        }

        printf("\n\n\nPerfomring the full:\n");
        int status = hist->Fit(this->fFit, "VSMRL+0", "")->Status();

        // decompose the fit function in its contributions
        if (this->bkgFuncName == "pol1") {
            this->fBkg = new TF1("pol1", Pol1, fitRangeMin, fitRangeMax, 2);
            this->fBkg->SetParameter(0, this->fFit->GetParameter(this->nSgnPars + 0));
            this->fBkg->SetParameter(1, this->fFit->GetParameter(this->nSgnPars + 1));
        } else if (this->bkgFuncName == "powex") {
            this->fBkg = new TF1("fPowEx", PowEx, fitRangeMin, fitRangeMax, 2);
            this->fBkg->SetParameter(0, this->fFit->GetParameter(this->nSgnPars + 0));
            this->fBkg->SetParameter(1, this->fFit->GetParameter(this->nSgnPars + 1));
        }

        if (this->sgnFuncName == "hat") {
            this->fHatThin = new TF1("fHatThin", Gaus, fitRangeMin, fitRangeMax, 3);
            this->fHatThin->SetParameter(0, this->fFit->GetParameter(0) * this->fFit->GetParameter(3));
            this->fHatThin->SetParameter(1, this->fFit->GetParameter(1));
            this->fHatThin->SetParameter(2, this->fFit->GetParameter(2));

            this->fHatWide = new TF1("fHatThin", Gaus, fitRangeMin, fitRangeMax, 3);
            this->fHatWide->SetParameter(0, this->fFit->GetParameter(0) * (1 - this->fFit->GetParameter(3)));
            this->fHatWide->SetParameter(1, this->fFit->GetParameter(1));
            this->fHatWide->SetParameter(2, this->fFit->GetParameter(2) * this->fFit->GetParameter(4));
        } else if (this->sgnFuncName == "gaus") {
            this->fSgn = new TF1("fSgn", Gaus, fitRangeMin, fitRangeMax, 3);
            this->fSgn->SetParameter(0, this->fFit->GetParameter(0));
            this->fSgn->SetParameter(1, this->fFit->GetParameter(1));
            this->fSgn->SetParameter(2, this->fFit->GetParameter(2));
        }

        return status;
    }

    void Draw(TVirtualPad *pad) {
        pad->cd();
        hist->GetYaxis()->SetRangeUser(0, 1.3 * hist->GetMaximum());
        gPad->DrawFrame(fitRangeMin, 0, fitRangeMax, 1.3 * hist->GetMaximum(),
                        Form("%s;%s;%s", this->hist->GetTitle(), this->hist->GetXaxis()->GetTitle(),
                             this->hist->GetYaxis()->GetTitle()));

        if (this->bkgFuncName == "pol1") {
            this->fBkg->SetNpx(300);
            this->fBkg->SetLineColor(kGray + 2);
            this->fBkg->Draw("same");
        } else if (this->bkgFuncName == "powex") {
            this->fBkg->SetNpx(300);
            this->fBkg->SetLineColor(kGray + 2);
            this->fBkg->Draw("same");
        }

        if (this->sgnFuncName == "hat") {
            this->fHatThin->SetLineColor(kMagenta + 3);
            this->fHatThin->SetNpx(300);
            this->fHatThin->Draw("same");

            this->fHatWide->SetNpx(300);
            this->fHatWide->SetLineColor(kAzure + 2);
            this->fHatWide->Draw("same");
        } else if (this->sgnFuncName == "gaus") {
            this->fSgn->SetNpx(300);
            this->fSgn->SetLineColor(kAzure + 2);
            this->fSgn->Draw("same");
        }

        fPrefit->SetLineStyle(9);
        fPrefit->SetLineColor(kGray + 2);
        fPrefit->Draw("same");

        fFit->SetLineColor(kRed);
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
        tl.DrawLatexNDC(
            .15, .85 - step * iStep++,
            Form("S(%.2f#sigma) = %.2f #pm %.2f", nSigma, this->GetSignal(nSigma), this->GetSignalUnc(nSigma)));

        tl.DrawLatexNDC(
            .15, .85 - step * iStep++,
            Form("B(%.2f#sigma) = %.2f #pm %.2f", nSigma, this->GetBackground(nSigma), this->GetBackgroundUnc(nSigma)));

        tl.DrawLatexNDC(.15, .85 - step * iStep++, Form("Counts = %.2f", this->GetCounts()));
        pad->Update();
    }

    double GetMean() {
        if (!fFit) return -1;
        if (this->sgnFuncName == "gaus" || this->sgnFuncName == "hat") return fFit->GetParameter(1);
        return -1;
    }

    double GetWidth() {
        if (!fFit) return -1;
        if (this->sgnFuncName == "gaus" || this->sgnFuncName == "hat") return fFit->GetParameter(2);
        return -1;
    }

    double GetWidthUnc() {
        if (!fFit) return -1;
        if (this->sgnFuncName == "gaus" || this->sgnFuncName == "hat") return fFit->GetParError(2);
        return -1;
    }

    double GetCounts() {
        if (!fFit) return -1;
        int firstBin = hist->GetXaxis()->FindBin(fitRangeMin * 1.0001);
        int lastBin = hist->GetXaxis()->FindBin(fitRangeMax * 0.9999);

        double totCounts = hist->Integral(firstBin, lastBin);
        double bkgCounts =
            fBkg->Integral(hist->GetXaxis()->GetBinLowEdge(firstBin), hist->GetXaxis()->GetBinLowEdge(lastBin + 1)) /
            this->hist->GetBinWidth(1);
        return totCounts - bkgCounts;
    }

    double GetSignal(double nSigma) {
        if (!fFit) return -1;

        double start = this->GetMean() - nSigma * this->GetWidth();
        double end = this->GetMean() + nSigma * this->GetWidth();

        if (this->sgnFuncName == "gaus") {
            return fSgn->Integral(start, end) / this->hist->GetBinWidth(1);
        } else if (this->sgnFuncName == "hat") {
            return (fHatThin->Integral(start, end) + fHatWide->Integral(start, end)) / this->hist->GetBinWidth(1);
        }
        return -1;
    }
    double GetChi2Ndf() { return fFit->GetChisquare() / fFit->GetNDF(); }

    double GetBackground(double nSigma) {
        if (!fFit) return -1;

        double start = this->GetMean() - nSigma * this->GetWidth();
        double end = this->GetMean() + nSigma * this->GetWidth();

        return fBkg->Integral(start, end) / this->hist->GetBinWidth(1);
    }

    double GetBackgroundUnc(double nSigma) {
        Int_t leftBand = this->hist->FindBin(this->GetMean() - 6 * this->GetWidth());
        Int_t rightBand = this->hist->FindBin(this->GetMean() + 6 * this->GetWidth());

        int start = this->hist->FindBin(this->fitRangeMin * 1.0001);
        int end = this->hist->FindBin(this->fitRangeMax * 0.9999);
        double SidebandBkg = this->hist->Integral(start, leftBand) + this->hist->Integral(rightBand, end);

        double sum2 = 0;
        for (Int_t i = start; i <= leftBand; i++) {
            sum2 += this->hist->GetBinError(i) * this->hist->GetBinError(i);
        }
        for (Int_t i = rightBand; i <= end; i++) {
            sum2 += this->hist->GetBinError(i) * this->hist->GetBinError(i);
        }

        return TMath::Sqrt(sum2) / SidebandBkg * this->GetBackground(nSigma);
    }

    double GetSignalUnc(double nSigma) {
        if (!fFit) return -1;
        if (this->GetSignal(nSigma) <= 0) return 0;

        return fFit->GetParError(0) / fFit->GetParameter(0) * this->GetSignal(nSigma);
    }

   private:
    TH1 *hist = nullptr;
    TF1 *fSgn = nullptr;
    TF1 *fHatThin = nullptr;
    TF1 *fHatWide = nullptr;
    TF1 *fBkg = nullptr;
    TF1 *fFit = nullptr;
    TF1 *fPrefit = nullptr;

    std::string sgnFuncName;
    std::string bkgFuncName;

    double (*sgnFunc)(double *x, double *par);
    double (*bkgFunc)(double *x, double *par);

    int nSgnPars;
    int nBkgPars;
    double fitRangeMin;
    double fitRangeMax;
};