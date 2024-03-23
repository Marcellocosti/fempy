#include <map>
#include <string>
#include <tuple>
#include <complex>

#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TMath.h"
#include "gsl/gsl_sf_dawson.h"

double Pol0(double *x, double *par) { return par[0]; }

double Pol1(double *x, double *par) { return Pol0(x, par) + par[1] * x[0]; }

double Pol2(double *x, double *par) { return Pol1(x, par) + par[2] * pow(x[0], 2); }

double Pol3(double *x, double *par) { return Pol2(x, par) + par[3] * pow(x[0], 3); }

double Pol4(double *x, double *par) { return Pol3(x, par) + par[4] * pow(x[0], 4); }

double Pol5(double *x, double *par) { return Pol4(x, par) + par[5] * pow(x[0], 5); }

double FlatPol3(double *x, double *par) { 
    return 1 + par[4] * (Pol3(x, par) - 1); 
}

double PowerLaw(double *x, double *par) { 
    //return par[0] * TMath::Exp(-par[1] * x[0]) * Pol1(x, &par[2]); 
    return par[0] + par[1] * TMath::Exp(-par[2] * x[0]) + Pol3(x, &par[3]); 
}

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
    double mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    if (x[0] < mpi) return 0;
    return par[0] * TMath::Sqrt(x[0] - mpi) * TMath::Exp(-1. * par[1] * (x[0] - mpi));
}

double Spline3(double *x, double *par){
    int numKnots = 10;
    double xKnots[numKnots];
    double yKnots[numKnots];
    for(int iKnot=0; iKnot<numKnots; iKnot++){
        xKnots[iKnot] = par[iKnot];
        yKnots[iKnot] = par[numKnots+iKnot];
    }
    TSpline3* sp3 = new TSpline3("sp3", xKnots, yKnots, numKnots, "");
    return sp3->Eval(x[0]);
}

double Spline5(double *x, double *par){
    int numKnots = 6;
    double xKnots[numKnots];
    double yKnots[numKnots];
    for(int iKnot=0; iKnot<numKnots; iKnot++){
        xKnots[iKnot] = par[iKnot];
        yKnots[iKnot] = par[numKnots+iKnot];
    }
    TSpline5* sp5 = new TSpline5("sp5", xKnots, yKnots, numKnots, "");
    return sp5->Eval(x[0]);
}

double Spline3Range(double *x, double *par){
    int numKnots = 10;
    double xKnots[numKnots];
    double yKnots[numKnots];
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

double SillKStar(double *x, double *par) {

  // x[0]: k*
  // par: [0] "normalisation" constant
  //      [1] width
  //      [2] mass

  if (x[0] < 0)
    return 0;

  double kstar = x[0];
  double massPion = 139.570;
  double massLambda = 1115.683;
  double thresh = massPion + massLambda;

  double motherMass = sqrt(kstar * kstar + massPion * massPion) + sqrt(kstar * kstar + massLambda * massLambda);

  if (motherMass < thresh)
    return 0;

  double width = par[1] * par[2] / sqrt(par[2] * par[2] - thresh * thresh);
  double arg0 = 2 * motherMass / TMath::Pi();
  double arg1 = sqrt(motherMass * motherMass - thresh * thresh) * width;
  double arg2 = pow(motherMass * motherMass - par[2] * par[2], 2.);
  double arg3 = pow(sqrt(motherMass * motherMass - thresh * thresh) * width, 2.);

  double jacobian = kstar / sqrt(kstar * kstar + massPion * massPion) + kstar / sqrt(kstar * kstar + massLambda * massLambda);
  return (par[0] * arg0 * arg1 / (arg2 + arg3)) * abs(jacobian);
}


//double Sill_WithPS(double *x, double *par) {
//
//  // x[0]: k*
//  // par: [0] "normalisation" constant
//  //      [1] width
//  //      [2] mass
//  //      [3] resonance pT
//  //      [4] kinetic decoupling temperature
//
//  double LambdaFactor = 290;
//
//  if (x[0] < 0)
//    return 0;
//
//  // double t = x[0];
//
//  double kstar = x[0];
//  double Thresh = MassPion + MassProton;
//
//  double MotherMass = sqrt(kstar * kstar + MassPion * MassPion) + sqrt(kstar * kstar + MassProton * MassProton);
//
//  if (MotherMass < Thresh)
//    return 0;
//
//  double width = par[1] * par[2] / sqrt(par[2] * par[2] - Thresh * Thresh);
//  double arg0 = 2 * MotherMass / TMath::Pi();
//  double arg1 = sqrt(MotherMass * MotherMass - Thresh * Thresh) * width;
//  double arg2 = pow(MotherMass * MotherMass - par[2] * par[2], 2.);
//  double arg3 = pow(sqrt(MotherMass * MotherMass - Thresh * Thresh) * width, 2.);
//
//  double ResonancePt = par[3];
//  double Temperature = par[4]; // orignal paper 160
//  double PhaseSpace = (MotherMass / sqrt(MotherMass * MotherMass + ResonancePt * ResonancePt)) * exp(-sqrt(MotherMass * MotherMass + ResonancePt * ResonancePt) / Temperature);
//  double PhaseSpaceNorm = exp(-sqrt(Thresh * Thresh + ResonancePt * ResonancePt) / Temperature) * Temperature;
//
//  // return (par[0]*arg0*arg1/(arg2 + arg3))*PhaseSpace/(PhaseSpaceNorm*ME_Value); // *abs(Derivative_dM_dkStar(kstar));
//
//  return (par[0] * arg0 * arg1 / (arg2 + arg3)) * PhaseSpace / (PhaseSpaceNorm)*abs(Derivative_dM_dkStar(kstar));
//
//}

const double FmToNu(5.067731237e-3);
const double Pi(3.141592653589793);
const std::complex<double> i(0,1);

double GeneralLednicky(const double &Momentum, const double &GaussR,
                       const complex<double> &ScattLenSin, const double &EffRangeSin) {

    // Taken from https://github.com/dimihayl/DLM/blob/c40f03eac38006f89eac8e5fa1533c9e48f2b455/CATS_Extentions/DLM_CkModels.cpp#L215
    if (GaussR != GaussR)
    {
        printf("\033[1;33mWARNING:\033[0m GeneralLednicky got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Radius = GaussR * FmToNu;
    const complex<double> IsLen1 = 1. / (ScattLenSin * FmToNu + 1e-64);
    const double eRan1 = EffRangeSin * FmToNu;

    double F1 = gsl_sf_dawson(2. * Momentum * Radius) / (2. * Momentum * Radius);
    double F2 = (1. - exp(-4. * Momentum * Momentum * Radius * Radius)) / (2. * Momentum * Radius);
    complex<double> ScattAmplSin = pow(IsLen1 + 0.5 * eRan1 * Momentum * Momentum - i * Momentum, -1.);

    double CkValue = 0.;
    CkValue += 0.5 * pow(abs(ScattAmplSin) / Radius, 2) * (1. - (eRan1) / (2 * sqrt(Pi) * Radius)) +
               2 * real(ScattAmplSin) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmplSin) * F2 / Radius;
    // double term1 = +2*real(ScattAmplSin)*F1/(sqrt(Pi)*Radius);
    // double term2 = -imag(ScattAmplSin)*F2/Radius;
    // so far this is the eq. for singlet only, w/o QS
    //  std::cout << "------- In General Lednicky analytical -------" << std::endl;
    //  printf("term1 = %.4f \n", term1);
    //  printf("term2 = %.4f \n", term2);
    //  std::cout << "----------------------------------------------" << std::endl;

    CkValue += 1;

    return CkValue;
}

double ComplexLednicky_Singlet_doublegaussian_lambda(double *x, double *par)
{
    // Taken from https://github.com/dimihayl/DLM/blob/c40f03eac38006f89eac8e5fa1533c9e48f2b455/CATS_Extentions/DLM_CkModels.cpp#L504C8-L504C53
    double kStar = x[0];
    double potPar0 = par[0];     // real part of the scattering length
    double potPar1 = par[1];     // imaginary part of the scattering length
    double potPar2 = par[2];     // effective range
    double sourcePar0 = par[3];  // radius of the first gaussian
    double sourcePar1 = par[4];  // radius of the second gaussian
    double sourcePar2 = par[5];  // relative weight of the two gaussians
    double sourcePar3 = par[6];  // normalization of the gaussians
    complex<double> ScatLen(potPar0, potPar1);
    return sourcePar3 * (sourcePar2 * GeneralLednicky(kStar, sourcePar0, ScatLen, potPar2) 
           + (1 - sourcePar2) * GeneralLednicky(kStar, sourcePar1, ScatLen, potPar2)) + 1. - sourcePar3;
}
