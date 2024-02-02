#include <map>
#include <string>
#include <tuple>

#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TMath.h"

double Pol0(double *x, double *par) { return par[0]; }

double Pol1(double *x, double *par) { return Pol0(x, par) + par[1] * x[0]; }

double Pol2(double *x, double *par) { return Pol1(x, par) + par[2] * pow(x[0], 2); }

double Pol3(double *x, double *par) { return Pol2(x, par) + par[3] * pow(x[0], 3); }

double Pol4(double *x, double *par) { return Pol3(x, par) + par[4] * pow(x[0], 4); }

double Pol5(double *x, double *par) { return Pol4(x, par) + par[5] * pow(x[0], 5); }

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

double BreitWigner(double *x, double *par) {
    double kstar = x[0];

    double yield = par[0];
    double mean = par[1];
    double gamma = par[2];

    return yield * gamma / TMath::Pi() / (gamma * gamma + (kstar - mean) * (kstar - mean));
}

// Convolution of a Breit-Wigner and a gaussian
double Voigt(double *x, double *par) {
    double kstar = x[0];

    double yield = par[0];
    double mean = par[1];
    double sigma = par[2];
    double gamma = par[3];

    return  yield * TMath::Voigt(kstar - mean, sigma, gamma);
}

double Exp(double *x, double *par) {
    // p0: total yield
    // p1: slope
    return par[0] * TMath::Exp(par[1] * x[0]);
}

double PowEx(double *x, double *par) {
    // p0: total yield
    // p1: slope
    Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    if (x[0] < mpi) return 0;
    return par[0] * TMath::Sqrt(x[0] - mpi) * TMath::Exp(-1. * par[1] * (x[0] - mpi));
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
