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

double BreitWignerKStar(double *x, double *par) {

  // p0: normalisation
  // p1: mass
  // p2: width

  if (x[0] < 0)
    return 0;

  double massPion = 139.57039;
  double massLambda = 1115.683;

  double kStarMJacobian = x[0]/sqrt( x[0]*x[0] + massPion*massPion ) + x[0]/sqrt( x[0]*x[0] + massLambda*massLambda );

  return BreitWigner(x, par) * abs(kStarMJacobian);

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
    int numKnots = 10;
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

double Spline3Range(double *x, double *par){
    int numKnots = 10;
    Double_t xKnots[numKnots];
    Double_t yKnots[numKnots];
    for(int iKnot=0; iKnot<numKnots; iKnot++){
        xKnots[iKnot] = par[iKnot];
        yKnots[iKnot] = par[numKnots+iKnot];
    }
    TSpline3* sp3 = new TSpline3("sp3", xKnots, yKnots, numKnots, "");
    
    //if(x[0]<95 || x[0]>) return 0;
    //else {
    //    return sp3->Eval(x[0]);
    //}
    return sp3->Eval(x[0]);
}

double Spline3Histo(double *x, double *par){

    TFile *histoFile = TFile::Open("/home/mdicostanzo/an/LPi/Analysis/SimAllMothersMerged.root", "r");
    TH1D *splineHisto = static_cast<TH1D*>(histoFile->Get("Pairs/hSE_2113122_NoDirectSigmaXi_smearednew"));
    TSpline3* sp3 = new TSpline3(splineHisto);
    
    return sp3->Eval(x[0]);

}

 double Landau(double *x, double *par)
 {
    // par[0]: norm
    // par[1]: mean
    // par[2]: sigma

    //if (par[2] <= 0) return 0;
    //double den = ::ROOT::Math::landau_pdf( (x[0]-par[1])/par[2] );
    //if (!par[0]) return den;
    //return den/par[2];
 
    if (par[2] <= 0) return 0;
    double den = ::ROOT::Math::landau_pdf( (x[0]-par[1])/par[2] );
    if (!par[0]) return den;
    return par[0]*den;
 
}

//Double_t Sill_WithPS(Double_t *x, Double_t *par) {
//
//  // x[0]: k*
//  // par: [0] "normalisation" constant
//  //      [1] width
//  //      [2] mass
//  //      [3] resonance pT
//  //      [4] kinetic decoupling temperature
//
//  Double_t LambdaFactor = 290;
//
//  if (x[0] < 0)
//    return 0;
//
//  // double t = x[0];
//
//  Double_t kstar = x[0];
//  Double_t Thresh = MassPion + MassProton;
//
//  double MotherMass = sqrt(kstar * kstar + MassPion * MassPion) + sqrt(kstar * kstar + MassProton * MassProton);
//
//  if (MotherMass < Thresh)
//    return 0;
//
//  Double_t Width = par[1] * par[2] / sqrt(par[2] * par[2] - Thresh * Thresh);
//  Double_t arg0 = 2 * MotherMass / TMath::Pi();
//  Double_t arg1 = sqrt(MotherMass * MotherMass - Thresh * Thresh) * Width;
//  Double_t arg2 = pow(MotherMass * MotherMass - par[2] * par[2], 2.);
//  Double_t arg3 = pow(sqrt(MotherMass * MotherMass - Thresh * Thresh) * Width, 2.);
//
//  Double_t ResonancePt = par[3];
//  Double_t Temperature = par[4]; // orignal paper 160
//  Double_t PhaseSpace = (MotherMass / sqrt(MotherMass * MotherMass + ResonancePt * ResonancePt)) * exp(-sqrt(MotherMass * MotherMass + ResonancePt * ResonancePt) / Temperature);
//  Double_t PhaseSpaceNorm = exp(-sqrt(Thresh * Thresh + ResonancePt * ResonancePt) / Temperature) * Temperature;
//
//  // return (par[0]*arg0*arg1/(arg2 + arg3))*PhaseSpace/(PhaseSpaceNorm*ME_Value); // *abs(Derivative_dM_dkStar(kstar));
//
//  return (par[0] * arg0 * arg1 / (arg2 + arg3)) * PhaseSpace / (PhaseSpaceNorm)*abs(Derivative_dM_dkStar(kstar));
//
//}

