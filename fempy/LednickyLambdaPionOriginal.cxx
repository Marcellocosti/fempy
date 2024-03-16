
#include <map>
#include <string>
#include <tuple>
#include <complex>

#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TMath.h"
#include "gsl/gsl_sf_dawson.h"

const double FmToNu(5.067731237e-3);
const double Pi(3.141592653589793);
const std::complex<double> i(0,1);

double GeneralLednicky(const double &Momentum, const double &GaussR,
                       const complex<double> &ScattLenSin, const double &EffRangeSin,
                       const complex<double> &ScattLenTri, const double &EffRangeTri,
                       const bool &SinOnly, const bool &QS, const bool &InverseScatLen)
{

    if (GaussR != GaussR)
    {
        printf("\033[1;33mWARNING:\033[0m GeneralLednicky got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Radius = GaussR * FmToNu;
    const complex<double> IsLen1 = InverseScatLen ? ScattLenSin / FmToNu : 1. / (ScattLenSin * FmToNu + 1e-64);
    const double eRan1 = EffRangeSin * FmToNu;
    const complex<double> IsLen3 = InverseScatLen ? ScattLenTri / FmToNu : 1. / (ScattLenTri * FmToNu + 1e-64);
    const double eRan3 = EffRangeTri * FmToNu;

    double F1 = gsl_sf_dawson(2. * Momentum * Radius) / (2. * Momentum * Radius);
    double F2 = (1. - exp(-4. * Momentum * Momentum * Radius * Radius)) / (2. * Momentum * Radius);
    complex<double> ScattAmplSin = pow(IsLen1 + 0.5 * eRan1 * Momentum * Momentum - i * Momentum, -1.);

    double CkValue = 0.;
    // double term1 = +2*real(ScattAmplSin)*F1/(sqrt(Pi)*Radius);
    // double term2 = -imag(ScattAmplSin)*F2/Radius;
    CkValue += 0.5 * pow(abs(ScattAmplSin) / Radius, 2) * (1. - (eRan1) / (2 * sqrt(Pi) * Radius)) +
               2 * real(ScattAmplSin) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmplSin) * F2 / Radius;
    // so far this is the eq. for singlet only, w/o QS
    //  std::cout << "------- In General Lednicky analytical -------" << std::endl;
    //  printf("term1 = %.4f \n", term1);
    //  printf("term2 = %.4f \n", term2);
    //  std::cout << "----------------------------------------------" << std::endl;

    // if we need to include the triplet, we add the term with a weight factor of 3 more than the singlet.
    // since however the correct norm. coeff. are 0.25 and 0.75 we need to divide by 4 to get the final result
    if (!SinOnly)
    {
        complex<double> ScattAmplTri = pow(IsLen3 + 0.5 * eRan3 * Momentum * Momentum - i * Momentum, -1.);
        CkValue += 3 * (0.5 * pow(abs(ScattAmplTri) / Radius, 2) * (1 - (eRan3) / (2 * sqrt(Pi) * Radius)) +
                        2 * real(ScattAmplTri) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmplTri) * F2 / Radius);
        CkValue *= 0.25;
    }
    // if we have to include QS we need to add a correction factor and normalize by factor of 1/2
    if (QS)
    {
        CkValue -= exp(-Radius * Radius * 4. * Momentum * Momentum);
        CkValue *= 0.5;
    }

    CkValue += 1;

    return CkValue;
}

double ComplexLednicky_Singlet_doublegaussian_lambda(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    complex<double> ScatLen(PotPar[0], PotPar[1]);
    return SourcePar[3] * (SourcePar[2] * GeneralLednicky(Momentum, SourcePar[0], ScatLen, PotPar[2], 0, 0, true, false, false) + (1 - SourcePar[2]) * GeneralLednicky(Momentum, SourcePar[1], ScatLen, PotPar[2], 0, 0, true, false, false)) + 1. - SourcePar[3];
}

double LednickyCover(double *x, double *par)
{
    return ComplexLednicky_Singlet_doublegaussian_lambda(x[0], &par[3], &par[0]);
}